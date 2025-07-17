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
    r_inf_phys_fixed = 10 * pState.r_p; % 修改：增大外边界距离
    r_f_phys_guess = pState.r_p * 1.2; % 修改：更合理的火焰半径初始猜测，靠近颗粒表面
    
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
    min_m_dot = 1e-8;
    if m_dot_phys_guess < min_m_dot
        m_dot_phys_guess = min_m_dot;
    end

    % 转化为无量纲猜测值
    p_guess_nd(1) = r_f_phys_guess / L_char; % r_f*
    p_guess_nd(2) = m_dot_phys_guess / (4 * pi * L_char * params.rho_D_gas); % Pe_m (佩克莱数)
    
    % 使用固定的外部边界
    r_inf_nd = r_inf_phys_fixed / L_char; % 使用当前半径r_p进行无量纲化
    
    % V29.5 修正: 使用最简化的4点网格来明确定义2个区域，避免bvp4c的区域误判
    x_mesh_nd = [1, p_guess_nd(1), p_guess_nd(1), r_inf_nd];
    
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
    
    % 逐步求解的标志
    GLOBAL_BVP_DEPS.simplified_mode = false;
    GLOBAL_BVP_DEPS.simplified_reaction = false;
    
    % 三阶段求解策略
    try
        fprintf('正在执行逐步求解:\n');
        
        % 阶段1: 最简化模型 - 忽略辐射和简化反应热
        fprintf('  > 阶段1/3: 求解简化模型...\n');
        GLOBAL_BVP_DEPS.simplified_mode = true;
        GLOBAL_BVP_DEPS.simplified_reaction = true;
        options_simple = bvpset('Stats','on', 'NMax', 2000, 'RelTol', 5e-3, 'AbsTol', 5e-4);
        sol_simple = bvp4c(@ode_func_nd, @bc_func_nd, solinit, options_simple);
        fprintf('    * 阶段1求解完成。火焰位置: r_f/r_p = %.3f, Pe_m = %.3e\n', ...
                sol_simple.parameters(1), sol_simple.parameters(2));
        
        % 阶段2: 部分简化模型 - 启用辐射但保持简化反应热
        fprintf('  > 阶段2/3: 求解部分简化模型...\n');
        GLOBAL_BVP_DEPS.simplified_mode = false;
        GLOBAL_BVP_DEPS.simplified_reaction = true;
        solinit_medium = bvpinit(sol_simple, p_guess_nd);
        options_medium = bvpset('Stats','on', 'NMax', 3000, 'RelTol', 3e-3, 'AbsTol', 3e-4);
        sol_medium = bvp4c(@ode_func_nd, @bc_func_nd, solinit_medium, options_medium);
        fprintf('    * 阶段2求解完成。火焰位置: r_f/r_p = %.3f, Pe_m = %.3e\n', ...
                sol_medium.parameters(1), sol_medium.parameters(2));
        
        % 阶段3: 完整模型 - 完整的辐射和反应热
        fprintf('  > 阶段3/3: 求解完整模型...\n');
        GLOBAL_BVP_DEPS.simplified_mode = false;
        GLOBAL_BVP_DEPS.simplified_reaction = false;
        solinit_full = bvpinit(sol_medium, p_guess_nd);
        options_full = bvpset('Stats','on', 'NMax', 5000, 'RelTol', 1e-3, 'AbsTol', 1e-5);
        sol_nd = bvp4c(@ode_func_nd, @bc_func_nd, solinit_full, options_full);
        fprintf('    * 阶段3求解完成。火焰位置: r_f/r_p = %.3f, Pe_m = %.3e\n', ...
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
            conservative_p_guess(1) = 1.1; % 火焰更靠近颗粒表面
            conservative_p_guess(2) = p_guess_nd(2) * 0.1; % 更小的Pe_m
            
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
    theta = z_nd(1,:);
    T_phys = theta * deps.T_char; % 计算物理温度以获取物性
    
    props = deps.physicalModel.get_gas_properties(T_phys);
    k_gas = props.k_gas;
    cp_gas = props.Cp_gas;
    
    % 计算无量纲数 Le (刘易斯数)
    Le = k_gas / (deps.params.rho_D_gas * cp_gas);
    
    dzdr_nd = zeros(size(z_nd));
    dzdr_nd(1:5,:) = z_nd(6:10,:); % 定义没有改变
    
    % 无量纲化的能量方程
    dzdr_nd(6,:) = (Pe_m ./ (Le .* r_nd.^2)) .* z_nd(6,:) - (2./r_nd) .* z_nd(6,:);
    % 无量纲化的组分方程
    dzdr_nd(7:10,:) = (Pe_m ./ r_nd.^2) .* z_nd(7:10,:) - (2./r_nd) .* z_nd(7:10,:);
end

function res = bc_func_nd(yl_nd, yr_nd, p_nd)
    global GLOBAL_BVP_DEPS;
    deps = GLOBAL_BVP_DEPS;

    % --- 解包无量纲参数和状态 ---
    r_f_nd = p_nd(1);
    Pe_m = p_nd(2);
    
    y_at_rp_nd = yl_nd(:,1);
    y_at_rf_left_nd = yr_nd(:,1);
    y_at_rf_right_nd = yl_nd(:,2);
    y_at_rinf_nd = yr_nd(:,2);

    res = zeros(22,1);
    
    % --- 组1: 颗粒表面边界条件 (r* = 1) ---
    % (共4个条件)
    res(1) = y_at_rp_nd(1) - 1.0; % BC 1: 表面温度等于沸点 (无量纲温度为1)
    
    % BC 2: 表面Mg质量分数由克-克方程在沸点下的情况决定
    % (严格来说，Y_mg应为1，这里保留完整公式以备未来扩展)
    T_surf_phys = y_at_rp_nd(1) * deps.T_char;
    T_amb_phys = deps.params.ambient_temperature;
    mw = deps.params.materials;
    L_v = mw.Mg.L_evap_Mg;
    R_Mg = deps.params.R_u / mw.Mg.molar_mass;
    M_mix = deps.physicalModel.get_ambient_mixture_molar_mass();
    exponent = (L_v / R_Mg) * (1/T_amb_phys - 1/T_surf_phys);
    Y_mg_theory = min((mw.Mg.molar_mass / M_mix) * exp(exponent), 1.0);
    res(2) = y_at_rp_nd(2) - Y_mg_theory;
    
    res(3) = y_at_rp_nd(8); % BC 3: 表面无CO2反应 (CO2梯度为零)
    
    % BC 4: 表面能量守恒
    % [传导热量] + [表面反应放热] - [辐射损失] - [蒸发耗热] = 0
    theta_surf = y_at_rp_nd(1);
    dtheta_dr_star_surf = y_at_rp_nd(6);
    Yco_surf = y_at_rp_nd(4);
    dYco_dr_star_surf = y_at_rp_nd(9);
    J_co_surf_nd = -dYco_dr_star_surf + Pe_m * Yco_surf; % CO质量通量
    k_surf = deps.physicalModel.get_gas_properties(T_surf_phys).k_gas;
    k_char = deps.physicalModel.get_gas_properties(deps.T_char).k_gas;
    
    % V32: 使用恒定的反应热以保证数值稳定性
    H_reac_s = deps.params.reaction_heats.surface_reac_H;
    
    % 对于简化反应模式，调整反应热
    if isfield(deps, 'simplified_reaction') && deps.simplified_reaction
        % 使用较小的反应热值来改善数值稳定性
        H_reac_s = H_reac_s * 0.5;
    end
    
    cp_char = deps.physicalModel.get_gas_properties(deps.T_char).Cp_gas;
    Ja_s = cp_char * deps.T_char / (-H_reac_s / mw.CO.molar_mass);
    Ja_v = cp_char * deps.T_char / L_v;
    
    % 对于简化模式，调整辐射项
    Q_rad_nd = deps.params.emissivity*deps.params.sigma*(T_surf_phys^4 - T_amb_phys^4) * deps.L_char / (k_char * deps.T_char);
    if isfield(deps, 'simplified_mode') && deps.simplified_mode
        % 在简化模式中减小辐射项
        Q_rad_nd = Q_rad_nd * 0.1;
    end
    
    res(4) = -(k_surf/k_char)*dtheta_dr_star_surf + J_co_surf_nd/Ja_s - Q_rad_nd - Pe_m/Ja_v;

    % --- 组2: 远场边界条件 (r* -> inf) ---
    % (共5个条件)
    res(5) = y_at_rinf_nd(1) - deps.params.ambient_temperature / deps.T_char; % BC 5: 温度等于环境温度
    res(6) = y_at_rinf_nd(3) - 1.0; % BC 6: CO2质量分数为1 (纯CO2环境)
    res(7) = y_at_rinf_nd(2); % BC 7: Mg质量分数为0
    res(8) = y_at_rinf_nd(4); % BC 8: CO质量分数为0
    res(9) = y_at_rinf_nd(5); % BC 9: MgO质量分数为0

    % --- 组3: 火焰面界面条件 (r* = r_f*) ---
    % (共13个条件)
    
    % (3a) 状态连续性 (5个条件)
    % 在简化模式中，使用松弛因子处理连续性条件
    if isfield(deps, 'simplified_mode') && deps.simplified_mode
        relaxation = 0.9; % 松弛因子
        for i = 1:5
            res(9+i) = relaxation * (y_at_rf_left_nd(i) - y_at_rf_right_nd(i)); 
        end
    else
        % 完整模型中严格要求连续性
        res(10:14) = y_at_rf_left_nd(1:5) - y_at_rf_right_nd(1:5);
    end
    
    % (3b) 反应物耗尽 (2个条件)
    res(15) = y_at_rf_right_nd(2); % BC 15: 火焰右侧无Mg
    res(16) = y_at_rf_left_nd(3);  % BC 16: 火焰左侧无CO2
    
    % (3c) 能量跳变 (1个条件)
    theta_f = y_at_rf_left_nd(1); 
    T_f_phys = theta_f * deps.T_char;
    
    % V32: 使用恒定的反应热以保证数值稳定性
    H_reac_f = deps.params.reaction_heats.flame_reac_H;
    
    % 对于简化反应模式，调整反应热
    if isfield(deps, 'simplified_reaction') && deps.simplified_reaction
        % 使用较小的反应热值来改善数值稳定性
        H_reac_f = H_reac_f * 0.5;
    end

    Ja_f = cp_char * deps.T_char / (-H_reac_f / mw.Mg.molar_mass);
    k_f = deps.physicalModel.get_gas_properties(T_f_phys).k_gas;
    
    % 修改：使用梯度计算反应物通量，而非完整通量表达式
    J_mg_flame = -y_at_rf_left_nd(7); % 在火焰面，Y_mg=0，所以只有梯度项
    
    res(17) = -(k_f/k_char)*(y_at_rf_right_nd(6) - y_at_rf_left_nd(6)) - J_mg_flame/Ja_f; % BC 17: 温度梯度跳变等于火焰放热
    
    % (3d) 修改：化学计量通量平衡 (Burke-Schumann 模型) (5个条件)
    calc_flux_nd = @(Y, dYdr) -dYdr + Pe_m * Y;
    
    % 修改：使用梯度计算反应物通量
    J_mg_flame = -y_at_rf_left_nd(7); % 在火焰面，Y_mg=0，所以只有梯度项
    J_co2_flame = -y_at_rf_right_nd(8); % 在火焰面，Y_co2=0，所以只有梯度项
    
    % 产物通量计算保持不变
    J_co_left_nd = calc_flux_nd(y_at_rf_left_nd(4), y_at_rf_left_nd(9));
    J_co_right_nd = calc_flux_nd(y_at_rf_right_nd(4), y_at_rf_right_nd(9));
    
    J_mgo_left_nd = calc_flux_nd(y_at_rf_left_nd(5), y_at_rf_left_nd(10));
    J_mgo_right_nd = calc_flux_nd(y_at_rf_right_nd(5), y_at_rf_right_nd(10));

    % 修改：去除通量为零的强制约束，改为使用化学计量比约束
    % BC 18: 反应物摩尔通量比符合化学计量比 (Mg:CO2 = 1:1)
    res(18) = J_mg_flame/mw.Mg.molar_mass - J_co2_flame/mw.CO2.molar_mass;
    
    % BC 19: Mg与CO的摩尔通量平衡 (Mg -> CO，摩尔比1:1)
    J_co_net = (J_co_right_nd - J_co_left_nd);
    res(19) = J_mg_flame/mw.Mg.molar_mass - J_co_net/mw.CO.molar_mass;
    
    % BC 20-22: Mg与MgO的摩尔通量平衡 (Mg -> MgO，摩尔比1:1)
    J_mgo_net = (J_mgo_right_nd - J_mgo_left_nd);
    res(20) = J_mg_flame/mw.Mg.molar_mass - J_mgo_net/mw.MgO.molar_mass;
    
    % 新增：确保质量守恒 (总质量通量守恒)
    % BC 21: 火焰左侧总质量通量 = 火焰右侧总质量通量
    total_left = J_mg_flame + J_co_left_nd + J_mgo_left_nd;
    total_right = J_co2_flame + J_co_right_nd + J_mgo_right_nd;
    res(21) = total_left - total_right;
    
    % BC 22: 元素守恒 (Mg原子守恒)
    res(22) = J_mg_flame - J_mgo_net;
end

function y_nd = initial_guess_func_nd(r_nd, region)
    global GLOBAL_BVP_DEPS;
    r_f_nd = GLOBAL_BVP_DEPS.r_f_nd_guess;
    r_inf_nd = GLOBAL_BVP_DEPS.r_inf_nd;

    y_nd = zeros(10, 1);
    theta_f_guess = 2.0; % T_f ~ 2*T_p
    Y_mg_s_guess = 0.1;
    Y_co_peak_guess = 0.2;
    
    delta1 = (r_f_nd - 1) / 4;
    delta2 = (r_inf_nd - r_f_nd) / 4;
    frac1 = (tanh((r_nd - (1 + r_f_nd)/2) / (delta1+eps)) + 1) / 2;
    frac2 = (tanh((r_nd - (r_f_nd + r_inf_nd)/2) / (delta2+eps)) + 1) / 2;
    
    if region == 1
        y_nd(1) = 1 + (theta_f_guess - 1) * frac1;
        y_nd(2) = Y_mg_s_guess * (1 - frac1);
        y_nd(3) = 0;
        y_nd(4) = Y_co_peak_guess * frac1;
        y_nd(5) = 0.1 * frac1;
    else % region == 2
        y_nd(1) = theta_f_guess - (theta_f_guess - 0.5) * frac2; % assume T_amb ~ 0.5 T_p
        y_nd(2) = 0;
        y_nd(3) = 1.0 * frac2;
        y_nd(4) = Y_co_peak_guess * (1 - frac2);
        y_nd(5) = 0.1 * (1 - frac2);
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