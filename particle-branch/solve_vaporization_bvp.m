function rate_info = solve_vaporization_bvp(bvp_deps)
% SOLVE_VAPORIZATION_BVP - V29 (最终版: 完全无量纲化)
%
% 核心思想:
% 通过选取特征长度(r_p), 特征温度(T_p)等，将整个BVP问题
% (微分方程，边界条件，初始猜测) 转化为在量级为1附近的无量纲空间中求解。
% 这能从根本上解决由于变量尺度差异巨大和物理刚性导致的数值不稳定
% 和雅可比矩阵病态(ill-conditioned)的问题。
%
% 全局依赖注入: 依然是解决深层调用链中上下文丢失问题的最可靠方法。
global GLOBAL_BVP_DEPS;

    % =========================================================================
    % 步骤1: 定义特征量和无量纲参数
    % =========================================================================
    pState = bvp_deps.pState;
    params = bvp_deps.params;
    
    % 特征量
    L_char = pState.r_p;  % 特征长度: 颗粒半径
    T_char = pState.T_p;  % 特征温度: 颗粒沸点
    
    % 存储特征量和物理依赖，供子函数使用
    GLOBAL_BVP_DEPS = bvp_deps;
    GLOBAL_BVP_DEPS.L_char = L_char;
    GLOBAL_BVP_DEPS.T_char = T_char;

    % =========================================================================
    % 步骤2: 构建无量纲初始猜测 (solinit)
    % =========================================================================
    
    % 物理猜测值
    % V29.7 修正: 使用健壮的火焰半径猜测值(位于颗粒和固定外边界中间)
    initial_radius = params.initial_diameter / 2;
    % 根据要求，将远场边界设置为原始颗粒半径的5倍
    r_inf_phys_fixed = 5 * pState.r_p; 
    r_f_phys_guess = pState.r_p * 1.5; % 修改：将火焰半径猜测增大到1.5倍颗粒半径，避免过于接近表面
    
    % 修改：基于能量平衡的蒸发速率猜测
    % 计算辐射热损失
    T_surf = pState.T_p;
    T_amb = params.ambient_temperature;
    A_p = 4 * pi * pState.r_p^2;
    sigma = params.sigma;
    emissivity = params.emissivity;
    L_v = params.materials.Mg.L_evap_Mg;
    
    % 辐射热损失
    Q_rad = emissivity * sigma * A_p * (T_surf^4 - T_amb^4);
    % 估计蒸发速率
    m_dot_phys_guess = Q_rad / L_v;
    
    % 如果计算出的值太小，使用一个合理的最小值
    min_m_dot = 1e-10; % 增大最小值，提高数值稳定性
    if m_dot_phys_guess < min_m_dot
        m_dot_phys_guess = min_m_dot;
    end

%%%%%无量纲初值猜测%%%%%ss
    p_guess_nd(1) = r_f_phys_guess / L_char; % r_f*
    p_guess_nd(2) = m_dot_phys_guess / (4 * pi * L_char * params.rho_D_gas); % Pe_m (佩克莱数)
    
    %  fprintf('  > 猜测质量流率分子: %.3e , 无量纲分母: %.3e \n', m_dot_phys_guess, 4 * pi * L_char * params.rho_D_gas);
    % 使用固定的外部边界
    r_inf_nd = r_inf_phys_fixed / L_char; % 使用当前半径r_p进行无量纲化

    %%%%%%%%网格进行加密和优化，尤其是火焰面附近%%%%%%%
        % 区域1：颗粒表面到火焰面 - 显著增加点数并强化火焰面附近加密
        n_points_r1 = 50; % 大幅增加点数以提高分辨率
        % 使用非线性分布在火焰面附近加密
        r1_ratio = 1.5; % 增强加密程度，使点更集中在火焰面附近
        r1_factor = (r1_ratio^(1/(n_points_r1-1)));
        r1_raw = zeros(1, n_points_r1);
        for i = 1:n_points_r1
            r1_raw(i) = (r1_ratio^((i-1)/(n_points_r1-1)) - 1) / (r1_ratio - 1);
        end
        % 将[0,1]范围映射到[1,p_guess_nd(1)]
        r1_points = 1 + r1_raw * (p_guess_nd(1) - 1);

        % 区域2：火焰面到远场边界 - 同样显著增加点数
        n_points_r2 = 50; % 增加远场区域点数
        % 使用非线性分布在火焰面附近加密
        r2_ratio = 1.8; % 增强加密程度，更好地捕捉火焰面附近梯度
        r2_raw = zeros(1, n_points_r2);
        for i = 1:n_points_r2
            r2_raw(i) = (r2_ratio^((i-1)/(n_points_r2-1)) - 1) / (r2_ratio - 1);
        end
        % 将[0,1]范围映射到[p_guess_nd(1),r_inf_nd]
        r2_points = p_guess_nd(1) + r2_raw * (r_inf_nd - p_guess_nd(1));

        % 确保使用相同的精确火焰点值
        r_f_exact = p_guess_nd(1);  % 使用精确的参数值，而不是可能受精度影响的数组值

        % 关键：明确地在数组中包含两个完全相同的火焰点，确保bvp4c识别为多区域
        x_mesh_nd = [r1_points(1:end-1), r_f_exact, r_f_exact, r2_points(2:end)];

        % 验证网格配置
        fprintf('网格点总数: %d\n', length(x_mesh_nd));
        fprintf('区域1(颗粒表面到火焰面)点数: %d\n', n_points_r1);
        fprintf('区域2(火焰面到远场)点数: %d\n', n_points_r2);
        
        % 验证重复点是否完全相同
        % 更新检查重复点的索引，因为点数增加了
        flame_idx = n_points_r1; % 火焰点索引对应于第一区域的末尾
        fprintf('重复点检查: r_f_nd(1) == r_f_nd(2): %d\n', x_mesh_nd(flame_idx) == x_mesh_nd(flame_idx+1));
        
    % V29.3 的修正保持不变
    GLOBAL_BVP_DEPS.r_f_nd_guess = p_guess_nd(1);
    GLOBAL_BVP_DEPS.r_inf_nd = r_inf_nd;
    solinit = bvpinit(x_mesh_nd, @initial_guess_func_nd, p_guess_nd);
    
    % =========================================================================
    % 步骤3: 实施逐步求解策略
    % =========================================================================
    
    % 保存原始猜测和配置，用于回退
    original_solinit = solinit;
    
    % 打印计算信息
    fprintf('开始边界值问题 (BVP) 求解...\n');
    fprintf('  > 颗粒半径: %.3e m, 外边界: %.3e m\n', L_char, r_inf_nd * L_char);
    fprintf('  > 初始火焰半径猜测: %.3e m (r_f/r_p = %.2f)\n', p_guess_nd(1) * L_char, p_guess_nd(1));
    fprintf('  > 初始蒸发速率猜测: %.3e kg/s (Pe_m = %.3e)\n', m_dot_phys_guess, p_guess_nd(2));
    
    % 生成并分析初始猜测的图形
    plot_initial_guess(solinit);
    
    % 逐步求解的标志
    GLOBAL_BVP_DEPS.simplified_mode = false;
    GLOBAL_BVP_DEPS.simplified_reaction = false;
    GLOBAL_BVP_DEPS.relaxation_level = 0.8; % 默认松弛级别
    
% 三阶段求解策略
    try
        fprintf('正在执行逐步求解:\n');
        
        % 阶段1: 最简化模型 - 高度松弛的极简模型
        fprintf('  > 阶段1/4: 求解高度松弛的极简模型...\n');
        GLOBAL_BVP_DEPS.simplified_mode = true;
        GLOBAL_BVP_DEPS.simplified_reaction = true;
        GLOBAL_BVP_DEPS.relaxation_level = 0.3; % 使用更强的松弛因子提高稳定性
        options_very_simple = bvpset('Stats','on', 'NMax', 6000, 'RelTol', 1e-1, 'AbsTol', 1e-2); 
        % 添加雅可比矩阵分析和监控
        try
            sol_very_simple = bvp4c(@ode_func_nd, @bc_func_nd, solinit, options_very_simple);
        catch ME
            if strcmp(ME.identifier, 'MATLAB:bvp4c:SingJac')
                fprintf('检测到奇异雅可比矩阵，正在分析原因...\n');
                analyze_jacobian(@bc_func_nd, solinit, p_guess_nd);
            end
            rethrow(ME);
        end
        fprintf('    * 阶段1求解完成。火焰位置: r_f/r_p = %.3f, Pe_m = %.3e\n', ...
                sol_very_simple.parameters(1), sol_very_simple.parameters(2));
                
        % 阶段2: 简化模型 - 中等松弛
        fprintf('  > 阶段2/4: 求解中等松弛的简化模型...\n');
        GLOBAL_BVP_DEPS.simplified_mode = true;
        GLOBAL_BVP_DEPS.simplified_reaction = true;
        GLOBAL_BVP_DEPS.relaxation_level = 0.5; % 使用更强的松弛因子
        
        % 使用重新构造的初始化，确保火焰面正确
        % 从上一个解中获取参数，但重新构造初始猜测
        new_p_guess = sol_very_simple.parameters;
        GLOBAL_BVP_DEPS.r_f_nd_guess = new_p_guess(1);
        solinit_simple = bvpinit(x_mesh_nd, @initial_guess_func_nd, new_p_guess);
        
        options_simple = bvpset('Stats','on', 'NMax', 6000, 'RelTol', 5e-2, 'AbsTol', 5e-3); 
        sol_simple = bvp4c(@ode_func_nd, @bc_func_nd, solinit_simple, options_simple);
        fprintf('    * 阶段2求解完成。火焰位置: r_f/r_p = %.3f, Pe_m = %.3e\n', ...
                sol_simple.parameters(1), sol_simple.parameters(2));
        
        % 阶段3: 部分简化模型 - 启用辐射但保持简化反应热
        fprintf('  > 阶段3/4: 求解部分简化模型...\n');
        GLOBAL_BVP_DEPS.simplified_mode = false;
        GLOBAL_BVP_DEPS.simplified_reaction = true;
        GLOBAL_BVP_DEPS.relaxation_level = 0.8; % 轻微松弛
        solinit_medium = bvpinit(sol_simple, p_guess_nd);
        options_medium = bvpset('Stats','on', 'NMax', 4000, 'RelTol', 5e-3, 'AbsTol', 5e-4);
        sol_medium = bvp4c(@ode_func_nd, @bc_func_nd, solinit_medium, options_medium);
        fprintf('    * 阶段3求解完成。火焰位置: r_f/r_p = %.3f, Pe_m = %.3e\n', ...
                sol_medium.parameters(1), sol_medium.parameters(2));
        
        % 阶段4: 完整模型 - 完整的辐射和反应热
        fprintf('  > 阶段4/4: 求解完整模型...\n');
        GLOBAL_BVP_DEPS.simplified_mode = false;
        GLOBAL_BVP_DEPS.simplified_reaction = false;
        GLOBAL_BVP_DEPS.relaxation_level = 0.9; % 最小松弛
        solinit_full = bvpinit(sol_medium, p_guess_nd);
        options_full = bvpset('Stats','on', 'NMax', 5000, 'RelTol', 1e-3, 'AbsTol', 1e-5);
        sol_nd = bvp4c(@ode_func_nd, @bc_func_nd, solinit_full, options_full);
        fprintf('    * 阶段4求解完成。火焰位置: r_f/r_p = %.3f, Pe_m = %.3e\n', ...
                sol_nd.parameters(1), sol_nd.parameters(2));
        
    catch ME
        % 如果逐步求解失败，尝试回退到先前成功的阶段
        fprintf('遇到错误: %s\n', ME.message);
        
        if exist('sol_medium', 'var')
            fprintf('回退到阶段2的解...\n');
            sol_nd = sol_medium;
        elseif exist('sol_simple', 'var')
            fprintf('回退到阶段1的解...\n');
            sol_nd = sol_simple;
        else
            % 所有阶段都失败，尝试调整初始猜测
            fprintf('所有阶段都失败，尝试不同的初始猜测...\n');
            
            % 尝试更保守的火焰位置和Pe_m值
            conservative_p_guess = p_guess_nd;
            conservative_p_guess(1) = 2.0; % 增大火焰位置猜测，使其远离颗粒表面
            conservative_p_guess(2) = p_guess_nd(2) * 0.5; % 使用较小的Pe_m，但不要太小
            
            GLOBAL_BVP_DEPS.simplified_mode = true;
            GLOBAL_BVP_DEPS.simplified_reaction = true;
            GLOBAL_BVP_DEPS.r_f_nd_guess = conservative_p_guess(1);
            
            try
                % 创建新的初始猜测
                conservative_solinit = bvpinit(x_mesh_nd, @initial_guess_func_nd, conservative_p_guess);
                sol_nd = bvp4c(@ode_func_nd, @bc_func_nd, conservative_solinit, options_simple);
                fprintf('使用保守猜测求解成功。\n');
            catch
                % 如果依然失败，清理并重新抛出原始错误
                clear global GLOBAL_BVP_DEPS;
                rethrow(ME);
            end
        end
    end
    
% 清理全局变量
    clear global GLOBAL_BVP_DEPS;

    % =========================================================================
    % 步骤4: 将结果重新物理化并返回
    % =========================================================================
    
    p_sol_nd = sol_nd.parameters;
    Pe_m_sol = p_sol_nd(2);
    
    % 重新物理化关键结果
    r_f_sol = p_sol_nd(1) * L_char;
    m_dot_sol = Pe_m_sol * (4 * pi * L_char * params.rho_D_gas);
    
    % 计算表面的CO通量 (需要物理量)
    y_surf_nd = deval(sol_nd, 1); % r*=1
    theta_surf = y_surf_nd(1);
    Yco_surf = y_surf_nd(4);
    dYco_dr_star_surf = y_surf_nd(9);
    
    dYco_dr_surf = dYco_dr_star_surf / L_char; % 物理梯度
    A_p = 4 * pi * L_char^2;
    m_dot_co_surf = -(A_p * params.rho_D_gas) * dYco_dr_surf + m_dot_sol * Yco_surf;
    
    rate_info = package_rate_info_func(m_dot_sol, m_dot_co_surf, r_f_sol, theta_surf * T_char, params);
    
    fprintf('BVP求解完成。火焰位置: %.3e m (r_f/r_p = %.2f), 蒸发速率: %.3e kg/s\n', ...
            r_f_sol, r_f_sol/L_char, m_dot_sol);
end

% =========================================================================
% 无量纲化的BVP核心函数
% =========================================================================

function dzdr_nd = ode_func_nd(r_nd, z_nd, region, p_nd)
    global GLOBAL_BVP_DEPS;
    deps = GLOBAL_BVP_DEPS;
    
    Pe_m = p_nd(2); % 佩克莱数
    theta = z_nd(1,:); % 无量纲温度
    T_phys = theta * deps.T_char; % 计算物理温度以获取物性
    
    props = deps.physicalModel.get_gas_properties(T_phys);
    k_gas = props.k_gas;
    cp_gas = props.Cp_gas;
    
    % 计算无量纲数 Le (刘易斯数)
    Le = k_gas / (deps.params.rho_D_gas * cp_gas);
    Le = 1 ;
    
    dzdr_nd = zeros(size(z_nd));
    dzdr_nd(1:5,:) = z_nd(6:10,:); % 定义没有改变
    
    % 无量纲化的能量方程
    dzdr_nd(6,:) = (Pe_m ./ (Le .* r_nd.^2)) .* z_nd(6,:) - (2./r_nd) .* z_nd(6,:);
    % 无量纲化的组分方程
    dzdr_nd(7:10,:) = (Pe_m ./ r_nd.^2) .* z_nd(7:10,:) - (2./r_nd) .* z_nd(7:10,:);
end

function res = bc_func_nd(ya, yb, p_nd)
    global GLOBAL_BVP_DEPS;
    deps = GLOBAL_BVP_DEPS;

    % --- 解包无量纲参数和状态 ---
    r_f_nd = p_nd(1);
    Pe_m = p_nd(2);
    
    y_at_rp_nd = ya(:,1);
    y_at_rf_left_nd = yb(:,1);
    y_at_rf_right_nd = ya(:,2);
    y_at_rinf_nd = yb(:,2);

    % 从全局设置读取松弛级别，如果未设置则使用默认值
    if isfield(deps, 'relaxation_level')
            base_relaxation = deps.relaxation_level;
        else
            base_relaxation = 0.8; % 默认松弛因子
    end
    
    % 为不同类型的条件设置不同的松弛级别
    relaxation = base_relaxation; % 基本松弛因子
    
    % 根据阶段调整守恒条件的松弛强度
    % 初始阶段使用更弱的守恒约束，后期阶段逐步加强
    if isfield(deps, 'simplified_mode') && deps.simplified_mode
        conservation_strength = 0.3; % 简化模式下使用较弱的守恒约束
    else
        conservation_strength = base_relaxation * 0.8; % 完整模式下逐步加强守恒约束
    end

    res = zeros(22,1); % 恢复到22个条件，满足bvp4c的期望
    
    % --- 组1: 颗粒表面边界条件 (r* = 1) ---
    res(1) = relaxation * (y_at_rp_nd(1) - 1.0); % BC 1: 表面温度等于沸点 (无量纲温度为1)
    
    % BC 2: 表面Mg质量分数由克-克方程在沸点下的情况决定
    % (严格来说，Y_mg应为1，这里保留完整公式以备未来扩展)
    T_surf_phys = y_at_rp_nd(1) * deps.T_char;
    T_amb_phys = deps.params.ambient_temperature;
    mw = deps.params.materials;
    L_v = mw.Mg.L_evap_Mg;
    R_Mg = deps.params.R_u / mw.Mg.molar_mass;
    M_mix = deps.physicalModel.get_ambient_mixture_molar_mass();
    % fprintf('M_mix: %f\n', M_mix);
    exponent = (L_v / R_Mg) * (1/T_amb_phys - 1/T_surf_phys);
    % fprintf('Y_mg_theory_theory: %f\n', (mw.Mg.molar_mass / M_mix) * exp(exponent));
    Y_mg_theory = min((mw.Mg.molar_mass / M_mix) * exp(exponent), 1.0);
    res(2) = y_at_rp_nd(2) - Y_mg_theory;
    
    res(3) = y_at_rp_nd(4); % BC 3: 表面CO浓度为零 (表面CO以无穷大速率反应)
    res(4) = y_at_rp_nd(5); % BC 4: 表面MgO浓度为零 (表面MgO以无穷大速率沉积)
    
    % BC 5: 表面能量守恒 (修改编号和考虑表面反应)
    % [传导热量] + [表面反应放热] - [辐射损失] - [蒸发耗热] = 0
    theta_surf = y_at_rp_nd(1);
    dtheta_dr_star_surf = y_at_rp_nd(6);
    % 注意：由于表面CO浓度为零，这里直接使用梯度计算通量
    dYco_dr_star_surf = y_at_rp_nd(9);
    dYmgo_dr_star_surf = y_at_rp_nd(10);
    
    % CO和MgO在表面的质量通量 (纯扩散，因为表面浓度为0)
    J_co_surf_nd = -dYco_dr_star_surf; 
    J_mgo_surf_nd = -dYmgo_dr_star_surf;
    
    k_surf = deps.physicalModel.get_gas_properties(T_surf_phys).k_gas;
    k_char = deps.physicalModel.get_gas_properties(deps.T_char).k_gas;
    
    % V32: 使用恒定的反应热以保证数值稳定性
    H_reac_s = deps.params.reaction_heats.surface_reac_H;
    
    % 确保surface_CO_reac_H参数存在，如果不存在则使用默认值
    if isfield(deps.params.reaction_heats, 'surface_CO_reac_H')
        H_reac_co = deps.params.reaction_heats.surface_CO_reac_H;
    else
        % 默认值：假设CO表面反应热与表面主反应热相当
        H_reac_co = H_reac_s;
    end
    
    cp_char = deps.physicalModel.get_gas_properties(deps.T_char).Cp_gas;
    Ja_s = cp_char * deps.T_char / (-H_reac_s / mw.CO.molar_mass);
    Ja_co = cp_char * deps.T_char / (-H_reac_co / mw.CO.molar_mass); % CO表面反应的无量纲数
    Ja_v = cp_char * deps.T_char / L_v;
    
    % 对于简化模式，调整辐射项
    Q_rad_nd = deps.params.emissivity*deps.params.sigma*(T_surf_phys^4 - T_amb_phys^4) * deps.L_char / (k_char * deps.T_char);
    if isfield(deps, 'simplified_mode') && deps.simplified_mode
        % 在简化模式中减小辐射项
        Q_rad_nd = Q_rad_nd * 0.1;
    end
    
    % 修改能量平衡方程：添加CO表面反应放热
    res(5) = -(k_surf/k_char)*dtheta_dr_star_surf + J_co_surf_nd/Ja_s + J_co_surf_nd/Ja_co - Q_rad_nd - Pe_m/Ja_v;

    % --- 组2: 远场边界条件 (r* -> inf) ---
    % (共5个条件)
    res(6) = relaxation * (y_at_rinf_nd(1) - deps.params.ambient_temperature / deps.T_char); % BC 6: 温度等于环境温度
    res(7) = relaxation * (y_at_rinf_nd(3) - 1.0); % BC 7: CO2质量分数为1 (纯CO2环境)
    
    % BC 8: 火焰面左侧所有组分质量分数之和为1（替换为更基本的物理约束，避免与元素守恒冗余）
    res(8) = relaxation * (sum(y_at_rf_left_nd(2:5)) - 1.0);
    res(9) = relaxation * y_at_rinf_nd(4); % BC 9: CO质量分数为0
    res(10) = relaxation * y_at_rinf_nd(5); % BC 10: MgO质量分数为0

    % --- 组3: 火焰面界面条件 (r* = r_f*) ---
    
    % (3a) 状态连续性 (5个条件) - 始终使用松弛因子
    flame_relaxation = relaxation; % 使用全局松弛因子
        for i = 1:5
        res(10+i) = flame_relaxation * (y_at_rf_left_nd(i) - y_at_rf_right_nd(i)); 
    end
    
    % (3b) 反应物耗尽 (2个条件) - 无穷大速率反应的关键条件
    res(16) = flame_relaxation * y_at_rf_right_nd(2); % BC 16: 火焰右侧无Mg
    res(17) = flame_relaxation * y_at_rf_left_nd(3);  % BC 17: 火焰左侧无CO2
    
    % (3c) 能量跳变 (1个条件) - 无穷大速率反应关键条件
    theta_f = y_at_rf_left_nd(1); 
    T_f_phys = theta_f * deps.T_char;
    
    % V32: 使用恒定的反应热以保证数值稳定性
    H_reac_f = deps.params.reaction_heats.flame_reac_H;
    


    Ja_f = cp_char * deps.T_char / (-H_reac_f / mw.Mg.molar_mass);
    k_f = deps.physicalModel.get_gas_properties(T_f_phys).k_gas;
    
    % 使用梯度计算反应物通量，而非完整通量表达式
    J_mg_flame = -y_at_rf_left_nd(7); % 在火焰面，Y_mg=0，所以只有梯度项
    
    res(18) = flame_relaxation * (-(k_f/k_char)*(y_at_rf_right_nd(6) - y_at_rf_left_nd(6)) - J_mg_flame/Ja_f); % BC 18: 温度梯度跳变等于火焰放热
    
    % (3d) 修改：简化的化学计量通量平衡 (只保留两个核心条件)
    calc_flux_nd = @(Y, dYdr) -dYdr + Pe_m * Y;
    
    % 使用梯度计算反应物通量
    J_mg_flame = -y_at_rf_left_nd(7); % 在火焰面，Y_mg=0，所以只有梯度项
    J_co2_flame = -y_at_rf_right_nd(8); % 在火焰面，Y_co2=0，所以只有梯度项
    
    % 产物通量计算
    J_co_left_nd = calc_flux_nd(y_at_rf_left_nd(4), y_at_rf_left_nd(9));
    J_co_right_nd = calc_flux_nd(y_at_rf_right_nd(4), y_at_rf_right_nd(9));
    
    J_mgo_left_nd = calc_flux_nd(y_at_rf_left_nd(5), y_at_rf_left_nd(10));
    J_mgo_right_nd = calc_flux_nd(y_at_rf_right_nd(5), y_at_rf_right_nd(10));

    % 使用强到弱不同程度的松弛因子处理化学计量和守恒条件
    chem_relaxation = relaxation * 0.9; % 主要化学计量条件的松弛因子
    weak_relaxation = conservation_strength; % 次要守恒条件的松弛因子，随求解阶段动态调整
    
    % 计算CO的净通量（修复J_co_net未定义问题）
    J_co_net = J_co_right_nd - J_co_left_nd;

    % 添加验证，删除人为修正
    %if sign(J_co_left_nd) == sign(J_co_right_nd)
     %   fprintf('警告: CO在火焰面两侧通量方向相同，J_left=%.3e, J_right=%.3e\n', J_co_left_nd, J_co_right_nd);
     %   fprintf('  > 这表明物理模型中可能存在通量计算或坐标系定义问题\n');
    %end

    % 替换为元素守恒条件 
    % 镁元素守恒：进入的Mg = 出去的MgO中的Mg
    J_mg_total = J_mg_flame;  % 左侧Mg通量
    % 修正MgO通量计算，考虑方向性（净流出 = 右侧 - 左侧）
    J_mgo_net = J_mgo_right_nd - J_mgo_left_nd;  % 两侧MgO净通量，与CO通量计算方式一致
    res(19) = chem_relaxation * (J_mg_total/mw.Mg.molar_mass - J_mgo_net/(mw.MgO.molar_mass/mw.Mg.molar_mass));
    
    % 碳元素守恒：进入的CO2中的C = 出去的CO中的C
    % 使用更弱的松弛因子来提高数值稳定性
    res(20) = weak_relaxation * 0.8 * (J_co2_flame/mw.CO2.molar_mass - J_co_net/mw.CO.molar_mass);
    
    % 替换为产物两侧扩散比例约束：假设产物向两侧扩散的摩尔比例相等
    % 使用更稳定的条件：CO产物总通量守恒
    % 使用更强的松弛因子确保CO通量方向正确
    res(21) = relaxation * (J_co_right_nd + J_co_left_nd);
    
    % 远场温度梯度趋近于零而非严格为零的条件
    % 使用一个足够小但非零的值，避免数值奇异性
    epsilon = 1e-4;  % 增大epsilon值，提高数值稳定性
    % 从全局变量中获取r_inf_nd和r_f_nd
    r_inf_nd_val = deps.r_inf_nd;  % 从全局依赖中获取
    r_f_nd_val = p_nd(1);  % 火焰位置从参数中获取
    r_ratio = r_inf_nd_val / r_f_nd_val;  % 远场与火焰面的比例
    
    % 远场温度梯度趋近于零而非严格为零的条件
    % 使用一个足够小但非零的值，避免数值奇异性
    epsilon = 1e-4;  % 增大epsilon值，提高数值稳定性
    % 从全局变量中获取r_inf_nd和r_f_nd
    r_inf_nd_val = deps.r_inf_nd;  % 从全局依赖中获取
    r_f_nd_val = p_nd(1);  % 火焰位置从参数中获取
    r_ratio = r_inf_nd_val / r_f_nd_val;  % 远场与火焰面的比例
    
    % 温度梯度条件
    temp_residual = y_at_rinf_nd(6) - epsilon/r_ratio^2;
    % 组分梯度条件 - 加权计入温度梯度条件
    comp_factor = 0.2;  % 组分梯度影响因子
    co2_residual = y_at_rinf_nd(8) + epsilon * 2.0 / r_inf_nd_val^2; % CO2梯度(负值)
    co_residual = y_at_rinf_nd(9) - epsilon * 1.5 / r_inf_nd_val^2;  % CO梯度(正值)
    mgo_residual = y_at_rinf_nd(10) - epsilon * 1.5 / r_inf_nd_val^2; % MgO梯度(正值)
    
    % 组合残差，将组分梯度约束整合到温度梯度约束中
    combined_residual = temp_residual + comp_factor * (abs(co2_residual) + abs(co_residual) + abs(mgo_residual));
    
    res(22) = weak_relaxation * combined_residual;

end

function y_nd = initial_guess_func_nd(x, region)
    global GLOBAL_BVP_DEPS;
    r_f_nd = GLOBAL_BVP_DEPS.r_f_nd_guess;
    r_inf_nd = GLOBAL_BVP_DEPS.r_inf_nd;
    
    % 加载材料参数，用于计算化学计量比
    deps = GLOBAL_BVP_DEPS;
    mw = deps.params.materials;

    y_nd = zeros(10, 1);
    theta_f_guess = 2.0; % T_f ~ 2*T_p
    
    % 自动判断区域，处理可能缺少region参数的情况
    if nargin < 2 || isempty(region)
        % 根据x的值判断所在区域
        if x <= r_f_nd
            region = 1; % 颗粒表面到火焰面
        else
            region = 2; % 火焰面到远场边界
        end
    end
    
    % 平滑过渡函数
    delta1 = (r_f_nd - 1) / 4;
    delta2 = (r_inf_nd - r_f_nd) / 4;
    frac1 = (tanh((x - (1 + r_f_nd)/2) / (delta1+eps)) + 1) / 2;
    frac2 = (tanh((x - (r_f_nd + r_inf_nd)/2) / (delta2+eps)) + 1) / 2;
    
    if region == 1
        % 区域1：颗粒表面到火焰面
        y_nd(1) = 1 + (theta_f_guess - 1) * frac1;  % 温度
        
        % 基于化学计量比计算组分分布
        % 火焰面处: 全是产物CO和MgO
        % 表面处: 主要是Mg
        
        % 表面Mg猜测值 (表面处的主要组分)
        Y_mg_surf = 0.95;
        
        % 火焰面产物
        % 按照反应 Mg + CO2 -> MgO + CO 的化学计量关系
        % 假设火焰面处产物的摩尔比是1:1，计算质量分数
        total_molar_mass = mw.CO.molar_mass + mw.MgO.molar_mass;
        Y_co_flame = mw.CO.molar_mass / total_molar_mass;
        Y_mgo_flame = mw.MgO.molar_mass / total_molar_mass;
        
        % 使用平滑过渡，确保组分总和为1
        y_nd(2) = Y_mg_surf * (1 - frac1);       % Mg 从表面到火焰面递减
        y_nd(3) = 0;                             % 区域1无CO2
        y_nd(4) = Y_co_flame * frac1;            % CO 从表面到火焰面递增
        y_nd(5) = Y_mgo_flame * frac1;           % MgO 从表面到火焰面递增
        
        % 确保质量分数之和为1
        sum_mass = y_nd(2) + y_nd(3) + y_nd(4) + y_nd(5);
        if sum_mass > 0
            scaling = 1.0 / sum_mass;
            y_nd(2:5) = y_nd(2:5) * scaling;
        end
        
    else % region == 2
        % 区域2：火焰面到远场
        y_nd(1) = theta_f_guess - (theta_f_guess - 0.5) * frac2; % 温度，远场约0.5T_p
        
        % 火焰面产物，与区域1相同
        total_molar_mass = mw.CO.molar_mass + mw.MgO.molar_mass;
        Y_co_flame = mw.CO.molar_mass / total_molar_mass;
        Y_mgo_flame = mw.MgO.molar_mass / total_molar_mass;
        
        % 组分分布
        y_nd(2) = 0;                                  % 区域2无Mg
        y_nd(3) = frac2;                              % CO2 从火焰面到远场递增
        y_nd(4) = Y_co_flame * (1 - frac2);           % CO 从火焰面到远场递减
        y_nd(5) = Y_mgo_flame * (1 - frac2);          % MgO 从火焰面到远场递减
        
        % 确保质量分数之和为1
        sum_mass = y_nd(2) + y_nd(3) + y_nd(4) + y_nd(5);
        if sum_mass > 0
            scaling = 1.0 / sum_mass;
            y_nd(2:5) = y_nd(2:5) * scaling;
        end
        
        % 为远场梯度提供非零初值，特别是在靠近远场边界时
        if region == 2 && x > 0.75 * r_inf_nd  % 显著扩大应用范围
            % 使用更大的epsilon值，确保梯度有足够的数值表达
            epsilon = 1e-3;  % 增大1000倍
            dist_factor = (r_inf_nd - x) / (0.25 * r_inf_nd);  % 0到1的因子，更平滑的过渡
            
            % 为所有梯度变量提供明确的非零初值
            % 温度梯度 - 向外减小
            y_nd(6) = epsilon * dist_factor * (1.0 + 0.1*sin(x*10));  % 添加微小扰动避免精确对称
            
            % 各组分梯度 - 提供物理合理的非零值并增大量级
            y_nd(7) = 0;  % Mg梯度(远场无Mg)
            y_nd(8) = -epsilon * dist_factor * 3.0;  % CO2梯度(向外增加)，强化负梯度
            y_nd(9) = epsilon * dist_factor * 2.0;   % CO梯度(向外减少)，强化正梯度
            y_nd(10) = epsilon * dist_factor * 2.0;  % MgO梯度(向外减少)，强化正梯度
            
            % 确保远场精确到达边界条件 - 在远场边界处特别处理
            if x > 0.98 * r_inf_nd
                % 温度梯度接近于特定的非零小值
                y_nd(6) = epsilon / r_inf_nd^2;
                
                % 组分梯度接近于目标的非零值
                y_nd(8) = -epsilon * 2.0 / r_inf_nd^2;  % CO2
                y_nd(9) = epsilon * 1.5 / r_inf_nd^2;   % CO
                y_nd(10) = epsilon * 1.5 / r_inf_nd^2;  % MgO
            end
        end
    end
end

function rate_info_out = package_rate_info_func(m_dot_mg_total, m_dot_co_surf, r_f, T_f, params)
    mw = params.materials;
    rate_info_out.dmdt_mg = -m_dot_mg_total;
    m_dot_mgo_surf = m_dot_co_surf * (mw.MgO.molar_mass / mw.CO.molar_mass);
    rate_info_out.dmdt_c = m_dot_co_surf * (mw.C.molar_mass / mw.CO.molar_mass);
    m_dot_mg_surf_equiv = m_dot_co_surf * (mw.Mg.molar_mass/mw.CO.molar_mass);
    m_dot_mg_flame = m_dot_mg_total - m_dot_mg_surf_equiv;
    m_dot_mgo_flame = m_dot_mg_flame * (mw.MgO.molar_mass / mw.Mg.molar_mass);
    rate_info_out.dmdt_mgo = m_dot_mgo_surf + m_dot_mgo_flame;
    rate_info_out.r_f = r_f;
    rate_info_out.T_f = T_f;
end 

function plot_initial_guess(solinit)
    % 绘制初始猜测的温度和组分分布图
    global GLOBAL_BVP_DEPS;
    r_f_nd = GLOBAL_BVP_DEPS.r_f_nd_guess;
    r_inf_nd = GLOBAL_BVP_DEPS.r_inf_nd;
    
    % 生成用于绘图的径向网格
    r_plot_nd1 = linspace(1, r_f_nd, 50);
    r_plot_nd2 = linspace(r_f_nd, r_inf_nd, 50);
    r_plot_nd = [r_plot_nd1, r_plot_nd2(2:end)];
    
    % 在每个点上计算初始猜测值
    y_plot = zeros(10, length(r_plot_nd));
    for i = 1:length(r_plot_nd)
        r_i = r_plot_nd(i);
        % 确定当前点所在区域
        if r_i <= r_f_nd
            region_i = 1;
        else
            region_i = 2;
        end
        % 使用初始猜测函数计算初始值
        y_plot(:,i) = initial_guess_func_nd(r_i, region_i);
    end
    
    % 创建新图窗口并绘图
    figure;
    
    % 绘制无量纲温度
    subplot(2,2,1);
    plot(r_plot_nd, y_plot(1,:), 'b-', 'LineWidth', 2);
    hold on;
    plot([r_f_nd, r_f_nd], [0, max(y_plot(1,:))*1.2], 'r--', 'LineWidth', 1);
    title('初始猜测: 无量纲温度分布');
    xlabel('径向距离 r/r_p');
    ylabel('无量纲温度 \theta');
    grid on;
    
    % 绘制反应物质量分数
    subplot(2,2,2);
    plot(r_plot_nd, y_plot(2,:), 'b-', 'LineWidth', 2); % Mg
    hold on;
    plot(r_plot_nd, y_plot(3,:), 'g-', 'LineWidth', 2); % CO2
    plot([r_f_nd, r_f_nd], [0, 1.2], 'r--', 'LineWidth', 1);
    title('初始猜测: 反应物质量分数');
    xlabel('径向距离 r/r_p');
    ylabel('质量分数');
    legend('Mg', 'CO2');
    grid on;
    
    % 绘制产物质量分数
    subplot(2,2,3);
    plot(r_plot_nd, y_plot(4,:), 'b-', 'LineWidth', 2); % CO
    hold on;
    plot(r_plot_nd, y_plot(5,:), 'g-', 'LineWidth', 2); % MgO
    plot([r_f_nd, r_f_nd], [0, 1.2], 'r--', 'LineWidth', 1);
    title('初始猜测: 产物质量分数');
    xlabel('径向距离 r/r_p');
    ylabel('质量分数');
    legend('CO', 'MgO');
    grid on;
    
    % 绘制组分之和，检查质量守恒
    subplot(2,2,4);
    sum_mass = sum(y_plot(2:5,:),1);
    plot(r_plot_nd, sum_mass, 'b-', 'LineWidth', 2);
    hold on;
    plot([r_f_nd, r_f_nd], [0, 1.2], 'r--', 'LineWidth', 1);
    plot(r_plot_nd, ones(size(r_plot_nd)), 'k--', 'LineWidth', 1);
    title('初始猜测: 组分质量分数总和');
    xlabel('径向距离 r/r_p');
    ylabel('总质量分数');
    ylim([0, 1.2]);
    grid on;
    
    % 调整图形布局
    set(gcf, 'Position', [100, 100, 900, 700]);
    
    % 兼容新旧版本MATLAB的标题添加
    try
        sgtitle('边界值问题初始猜测分析'); % 新版MATLAB用sgtitle
    catch
        try
            suptitle('边界值问题初始猜测分析'); % 旧版MATLAB用suptitle
        catch
            % 如果两者都不可用，就不添加总标题
        end
    end
    
    % 保存图形到文件
    saveas(gcf, 'initial_guess_analysis.png');
    fprintf('初始猜测分析图已保存到 initial_guess_analysis.png\n');
end

function analyze_jacobian(bcfun, solinit, p_nd)
    try
        fprintf('正在计算边界条件函数的雅可比矩阵...\n');
        
        % 从solinit中提取所需信息
        x = solinit.x;
        y = solinit.y;
        
        % 找到重复点的位置，确定这是一个多区域问题
        duplicate_points = find(diff(x) == 0);
        
        if ~isempty(duplicate_points)
            fprintf('检测到多区域BVP结构，重复点在位置 %d\n', duplicate_points);
            
            % 构造多区域边界条件输入
            n = length(duplicate_points) + 1;  % 区域数
            yl = zeros(size(y,1), n);  % 每个区域的左边界
            yr = zeros(size(y,1), n);  % 每个区域的右边界
            
            % 第一个区域
            yl(:,1) = y(:,1);
            yr(:,1) = y(:,duplicate_points(1));
            
            % 中间区域（如果有多于2个区域）
            for i = 2:n-1
                yl(:,i) = y(:,duplicate_points(i-1)+1);
                yr(:,i) = y(:,duplicate_points(i));
            end
            
            % 最后一个区域
            yl(:,n) = y(:,duplicate_points(end)+1);
            yr(:,n) = y(:,end);
        else
            fprintf('检测到单区域BVP结构...\n');
            
            % 使用常规方法提取左右边界
            yl = y(:, 1);
            yr = y(:, end);
            
            % 构造多区域形式以适应bc_func_nd
            fprintf('将单区域形式转换为多区域形式以适应bc_func_nd...\n');
            yl_multi = zeros(size(y,1), 2);
            yr_multi = zeros(size(y,1), 2);
            
            yl_multi(:,1) = yl;  % 表面条件
            yr_multi(:,2) = yr;  % 远场条件
            
            % 估计火焰面的值（中间点）
            mid_idx = ceil(length(x)/2);
            y_mid = y(:,mid_idx);
            
            yr_multi(:,1) = y_mid;  % 火焰面左侧
            yl_multi(:,2) = y_mid;  % 火焰面右侧
            
            yl = yl_multi;
            yr = yr_multi;
        end
        
        % 定义偏导数步长
        h = 1e-6;
        
        % 计算基准残差
        fprintf('调用bc_func_nd计算残差...\n');
        base_res = bcfun(yl, yr, p_nd);
        
        % 输出基本信息
        fprintf('边界条件数量: %d\n', length(base_res));
        if ~iscell(solinit.y)
            fprintf('警告: 正在单区域模式下分析多区域边界条件函数。结果可能不准确。\n');
        end
        
        % 余下的雅可比矩阵分析代码不变...
        n_eqs = length(base_res);
        n_vars_y = numel(yl) + numel(yr);
        n_vars_p = length(p_nd);
        
        Jy = zeros(n_eqs, n_vars_y);
        Jp = zeros(n_eqs, n_vars_p);
        
        % 对于多区域输入，我们需要特殊处理扰动
        if size(yl,2) > 1
            % 处理多区域左边界
            idx = 1;
            for r = 1:size(yl,2)
                for i = 1:size(yl,1)
                    yl_perturbed = yl;
                    yl_perturbed(i,r) = yl(i,r) + h;
                    res_perturbed = bcfun(yl_perturbed, yr, p_nd);
                    Jy(:, idx) = (res_perturbed - base_res) / h;
                    idx = idx + 1;
                end
            end
            
            % 处理多区域右边界
            for r = 1:size(yr,2)
                for i = 1:size(yr,1)
                    yr_perturbed = yr;
                    yr_perturbed(i,r) = yr(i,r) + h;
                    res_perturbed = bcfun(yl, yr_perturbed, p_nd);
                    Jy(:, idx) = (res_perturbed - base_res) / h;
                    idx = idx + 1;
                end
            end
        else
            % 原始单区域处理方式
            for i = 1:length(yl)
                yl_perturbed = yl;
                yl_perturbed(i) = yl(i) + h;
                res_perturbed = bcfun(yl_perturbed, yr, p_nd);
                Jy(:, i) = (res_perturbed - base_res) / h;
            end
            
            for i = 1:length(yr)
                yr_perturbed = yr;
                yr_perturbed(i) = yr(i) + h;
                res_perturbed = bcfun(yl, yr_perturbed, p_nd);
                Jy(:, length(yl) + i) = (res_perturbed - base_res) / h;
            end
        end
        
        % 对参数的雅可比矩阵
        for i = 1:length(p_nd)
            p_perturbed = p_nd;
            p_perturbed(i) = p_nd(i) + h;
            res_perturbed = bcfun(yl, yr, p_perturbed);
            Jp(:, i) = (res_perturbed - base_res) / h;
        end
        
        % 分析雅可比矩阵
        cond_Jy = cond(Jy);
        fprintf('雅可比矩阵条件数: %.3e\n', cond_Jy);
        
        [~, S, ~] = svd(Jy);
        s = diag(S);
        min_sv = min(s);
        max_sv = max(s);
        fprintf('最小奇异值: %.3e, 最大奇异值: %.3e, 比值: %.3e\n', min_sv, max_sv, max_sv/min_sv);
        
        % 检查线性相关性
        tol = 1e-10;
        rank_Jy = rank(Jy, tol);
        fprintf('雅可比矩阵的秩: %d (应为 %d)\n', rank_Jy, n_eqs);
        
        if rank_Jy < n_eqs
            fprintf('检测到 %d 个线性相关的行!\n', n_eqs - rank_Jy);
            
            % 使用SVD方法更稳健地识别线性相关行
            [U, S, ~] = svd(Jy);
            s = diag(S);
            
            % 找到接近零的奇异值对应的左奇异向量
            small_sv_idx = find(abs(s) < tol * max(abs(s)));
            
            % 分析每个线性相关的组合
            fprintf('可能线性相关的边界条件组合:\n');
            for i = 1:length(small_sv_idx)
                % 获取对应的左奇异向量，它表示线性组合系数
                v = U(:, small_sv_idx(i));
                
                % 找出贡献最大的几个条件
                [~, idx] = sort(abs(v), 'descend');
                top_contributors = idx(1:min(5, length(idx))); % 选择前5个贡献最大的条件
                
                fprintf('  组合 #%d 主要贡献条件: ', i);
                for j = 1:length(top_contributors)
                    cond_idx = top_contributors(j);
                    fprintf('#%d(系数=%.3f) ', cond_idx, v(cond_idx));
                end
                fprintf('\n');
                
                % 分析边界条件之间的关系
                if length(top_contributors) >= 2
                    main_cond = top_contributors(1);
                    second_cond = top_contributors(2);
                    ratio = v(second_cond) / v(main_cond);
                    fprintf('    条件 #%d 约等于 %.3f * 条件 #%d\n', second_cond, -ratio, main_cond);
                end
            end
            
            % 仍然保留原来的相关性图分析
            figure;
            corr_matrix = corr(Jy');
            imagesc(corr_matrix);
            colorbar;
            title('边界条件雅可比矩阵行相关性分析');
            xlabel('条件索引');
            ylabel('条件索引');
            set(gca, 'XTick', 1:n_eqs, 'YTick', 1:n_eqs);
            saveas(gcf, 'jacobian_correlation.png');
            fprintf('雅可比矩阵相关性热图已保存到 jacobian_correlation.png\n');
        end
        
        % 在SVD分析后添加以下代码
        % 识别对条件数贡献最大的元素
        [max_row, max_col] = find(abs(Jy) == max(abs(Jy(:))));
        fprintf('最大元素(%.3e)位于：行 %d (条件 #%d), 列 %d\n', max(abs(Jy(:))), max_row, max_row, max_col);

        [min_elem_val, min_elem_idx] = min(abs(Jy(abs(Jy)>0)));  % 非零最小元素
        [min_row, min_col] = ind2sub(size(Jy), min_elem_idx);
        fprintf('最小非零元素(%.3e)位于：行 %d (条件 #%d), 列 %d\n', min_elem_val, min_row, min_row, min_col);
        
        % 增强雅可比分析 - 添加SVD可视化
        figure;
        semilogy(s, 'o-');
        title('雅可比矩阵奇异值分布');
        xlabel('索引'); ylabel('奇异值');
        grid on;
        saveas(gcf, 'jacobian_singular_values.png');
        fprintf('雅可比矩阵奇异值分布已保存到 jacobian_singular_values.png\n');
        
        % 分析变量敏感度并输出详细物理含义
        [U, S, V] = svd(Jy);
        fprintf('\n变量敏感度分析 (右奇异向量最大分量):\n');
        
        % 定义变量名称和所属边界
        var_names = {'温度', 'Y_Mg', 'Y_CO2', 'Y_CO', 'Y_MgO', 'd温度/dr', 'dY_Mg/dr', 'dY_CO2/dr', 'dY_CO/dr', 'dY_MgO/dr'};
        regions = {'颗粒表面(左)', '火焰面左侧(右)', '火焰面右侧(左)', '远场边界(右)'};
        params = {'r_f_nd (火焰位置)', 'Pe_m (佩克莱数)'};
        
        % 创建变量编号到物理含义的映射
        var_to_phys = cell(n_vars_y + n_vars_p, 1);
        var_to_region = cell(n_vars_y + n_vars_p, 1);
        var_to_value = zeros(n_vars_y + n_vars_p, 1);
        
        % 填充变量映射
        var_idx = 1;
        % 左边界变量
        for r = 1:size(yl,2)
            for i = 1:size(yl,1)
                var_to_phys{var_idx} = var_names{i};
                var_to_region{var_idx} = regions{2*r-1};
                var_to_value(var_idx) = yl(i,r);
                var_idx = var_idx + 1;
            end
        end
        
        % 右边界变量
        for r = 1:size(yr,2)
            for i = 1:size(yr,1)
                var_to_phys{var_idx} = var_names{i};
                var_to_region{var_idx} = regions{2*r};
                var_to_value(var_idx) = yr(i,r);
                var_idx = var_idx + 1;
            end
        end
        
        % 参数变量
        for i = 1:length(p_nd)
            var_to_phys{var_idx} = params{i};
            var_to_region{var_idx} = '参数';
            var_to_value(var_idx) = p_nd(i);
            var_idx = var_idx + 1;
        end
        
        % 输出前3个最小奇异值对应的右奇异向量
        for i = 1:min(3, size(V,2))
            [~, var_idx] = sort(abs(V(:,end-i+1)), 'descend');
            fprintf('  奇异值 #%d (%.3e) 对应的最敏感变量:\n', size(S,2)-i+1, s(end-i+1));
            for j = 1:min(5, length(var_idx))
                v_idx = var_idx(j);
                fprintf('    #%d: %s at %s (值=%.3e, 系数=%.3f)\n', v_idx, var_to_phys{v_idx}, var_to_region{v_idx}, var_to_value(v_idx), V(v_idx,end-i+1));
            end
        end
        
        % 特别关注变量#40
        if n_vars_y + n_vars_p >= 40
            fprintf('\n特别分析变量 #40:\n');
            fprintf('  物理含义: %s\n', var_to_phys{40});
            fprintf('  所在区域: %s\n', var_to_region{40});
            fprintf('  当前值: %.3e\n', var_to_value(40));
            fprintf('  在最小奇异值向量中的系数: %.3f\n', V(40,end));
            
            % 查找与变量#40关联最强的边界条件
            [~, cond_idx] = sort(abs(Jy(:,40)), 'descend');
            fprintf('  与该变量关联最强的边界条件:\n');
            for j = 1:min(3, length(cond_idx))
                fprintf('    条件 #%d (系数=%.3e)\n', cond_idx(j), Jy(cond_idx(j),40));
            end
        end
        
        % 分析贡献最大的条件对 - 添加变量存在性检查
        fprintf('\n对病态条件贡献最大的边界条件对:\n');
        % 首先检查是否存在small_sv_idx变量
        if ~exist('small_sv_idx', 'var') || isempty(small_sv_idx)
            fprintf('  矩阵是满秩的，没有检测到线性相关行。\n');
            fprintf('  可能的问题来源于数值计算或量级不平衡，而非条件线性相关。\n');
            
            % 分析最小奇异值对应的左奇异向量
            [~, min_sv_idx] = min(s);
            v_min = U(:, min_sv_idx);
            [~, sorted_idx] = sort(abs(v_min), 'descend');
            fprintf('  最小奇异值 %.3e 对应的主要条件: ', s(min_sv_idx));
            for j = 1:min(5, length(sorted_idx))
                fprintf('#%d(系数=%.3f) ', sorted_idx(j), v_min(sorted_idx(j)));
            end
            fprintf('\n');
            
            % 分析具体条件的相互关系和物理含义
            if length(sorted_idx) >= 2
                cond1 = sorted_idx(1);
                cond2 = sorted_idx(2);
                fprintf('  >>> 建议检查条件 #%d 和 #%d 之间的物理关系或数值稳定性\n', cond1, cond2);
            end
        else
            % 原有代码保持不变
            for i = 1:min(2, length(small_sv_idx))
                v = U(:, small_sv_idx(i));
                [~, sorted_idx] = sort(abs(v), 'descend');
                fprintf('  奇异值 #%d (%.3e) 对应的主要条件: ', small_sv_idx(i), s(small_sv_idx(i)));
                for j = 1:min(3, length(sorted_idx))
                    fprintf('#%d(系数=%.3f) ', sorted_idx(j), v(sorted_idx(j)));
                end
                fprintf('\n');
                
                % 分析具体条件的相互关系和物理含义
                if length(sorted_idx) >= 2
                    cond1 = sorted_idx(1);
                    cond2 = sorted_idx(2);
                    fprintf('  >>> 建议检查条件 #%d 和 #%d 之间的物理关系或冗余性\n', cond1, cond2);
                end
            end
        end

        % 计算每个条件的灵敏度
        cond_sensitivity = zeros(n_eqs, 1);
        for i = 1:n_eqs
            Jy_temp = Jy;
            Jy_temp(i,:) = [];  % 移除一行
            cond_sensitivity(i) = cond(Jy_temp);
        end
        [~, most_sensitive] = min(cond_sensitivity);
        fprintf('移除条件 #%d 后条件数最小(%.3e)\n', most_sensitive, cond_sensitivity(most_sensitive));
        
    catch ME
        fprintf('雅可比矩阵分析过程中出错: %s\n', ME.message);
        % 打印更详细的错误信息
        fprintf('错误位置: %s\n', getReport(ME, 'extended'));
    end
end 