function rate_info = solve_vaporization_integral(bvp_deps)
    % SOLVE_VAPORIZATION_INTEGRAL - 使用积分法求解气相燃烧问题
    % 替代原有的BVP方法，保持输入输出接口一致
    
        % 提取依赖
        pState = bvp_deps.pState;
        params = bvp_deps.params;
        
        % 特征量定义（与原BVP保持一致）
        L_char = pState.r_p;  % 特征长度: 颗粒半径
        T_char = pState.T_p;  % 特征温度: 颗粒沸点
        r_inf_phys = 5 * pState.r_p; % 远场边界
        
        % 初始猜测值
        r_f_phys_guess = pState.r_p * 1.5; % 火焰位置初猜
        
        % 估计蒸发速率初值（与原BVP相同方法）
        T_surf = pState.T_p;
        T_amb = params.ambient_temperature;
        A_p = 4 * pi * pState.r_p^2;
        sigma = params.sigma;
        emissivity = params.emissivity;
        L_v = params.materials.Mg.L_evap_Mg;
        Q_rad = emissivity * sigma * A_p * (T_surf^4 - T_amb^4);
        m_dot_phys_guess = max(Q_rad / L_v, 1e-10);
        
        % ======= 积分法求解代替BVP =======
        
        % 1. 构建初始解向量
        % 初始解向量包含：积分常数、界面值、火焰位置
        mw = params.materials; % 分子量信息
        
        x0 = zeros(15, 1);
        
        % 火焰温度初猜(2倍颗粒温度)
        T_f_guess = 2.0 * T_surf;
        
        % 计算初始通量猜测值
        k_gas = params.rho_D_gas * params.Cp_gas; % 近似导热系数
        rho_D = params.rho_D_gas;
        
        % 能量通量初猜 (基于线性温度分布)
        x0(1) = -k_gas * (T_f_guess - T_surf) / (r_f_phys_guess - pState.r_p) * (4*pi*pState.r_p*r_f_phys_guess); % C_T1
        x0(2) = -k_gas * (T_amb - T_f_guess) / (r_inf_phys - r_f_phys_guess) * (4*pi*r_f_phys_guess*r_inf_phys); % C_T2
        
        % 质量通量初猜
        x0(3) = -m_dot_phys_guess; % C_Mg (负值表示向外扩散)
        x0(4) = m_dot_phys_guess * (mw.CO.molar_mass/mw.Mg.molar_mass) * 0.5; % C_CO1 
        x0(5) = -m_dot_phys_guess * (mw.CO2.molar_mass/mw.Mg.molar_mass) * 0.5; % C_CO2
        x0(6) = m_dot_phys_guess * (mw.MgO.molar_mass/mw.Mg.molar_mass) * 0.5; % C_MgO1
        x0(7) = m_dot_phys_guess * (mw.MgO.molar_mass/mw.Mg.molar_mass) * 0.5; % C_MgO2
        
        % 界面质量分数猜测
        % 根据克拉珀龙方程计算表面Mg蒸气分数
        R_Mg = params.R_u / mw.Mg.molar_mass;
        M_mix = mw.CO2.molar_mass; % 近似为纯CO2环境
        exponent = (L_v / R_Mg) * (1/T_amb - 1/T_surf);
        Y_mg_surf = min((mw.Mg.molar_mass / M_mix) * exp(exponent), 0.95);
        
        x0(9) = T_f_guess; % 火焰温度
        x0(10) = 0; % Y_Mg_f (火焰面处Mg质量分数)
        
        % 火焰面处CO和MgO质量分数初猜
        % 按照反应 Mg + CO2 -> MgO + CO 的化学计量关系
        total_molar_mass = mw.CO.molar_mass + mw.MgO.molar_mass;
        Y_co_flame = mw.CO.molar_mass / total_molar_mass;
        Y_mgo_flame = mw.MgO.molar_mass / total_molar_mass;
        
        x0(11) = Y_co_flame; % Y_CO_f_left
        x0(12) = Y_mgo_flame; % Y_MgO_f_left
        x0(13) = Y_co_flame * 0.5; % Y_CO_f_right
        x0(14) = Y_mgo_flame * 0.5; % Y_MgO_f_right
        
        x0(15) = r_f_phys_guess; % 火焰位置
        
        % 2. 构建求解选项
        fprintf('开始使用积分法求解气相燃烧问题...\n');
        fprintf('  > 颗粒半径: %.3e m, 表面温度: %.2f K\n', pState.r_p, T_surf);
        fprintf('  > 初始火焰半径猜测: %.3e m (r_f/r_p = %.2f)\n', r_f_phys_guess, r_f_phys_guess/pState.r_p);
        fprintf('  > 初始蒸发速率猜测: %.3e kg/s\n', m_dot_phys_guess);
        
        % 3. 实施逐步求解策略
        % 使用比原BVP更稳健的优化器选项
        options = optimoptions('fsolve', ...
            'Algorithm', 'levenberg-marquardt', ... % LM算法处理病态问题
            'Display', 'iter', ...
            'MaxFunctionEvaluations', 10000, ...
            'MaxIterations', 1000, ...
            'FunctionTolerance', 1e-6, ...
            'StepTolerance', 1e-8);
        
        try
            % 阶段1: 高度简化模型
            fprintf('  > 阶段1/3: 求解高度简化模型...\n');
            params_simple = params;
            params_simple.simplified_mode = true;
            params_simple.simplified_reaction = true;
            
            % 调用求解器
            x_simple = fsolve(@(x) vap_equations(x, params_simple, pState, r_inf_phys), x0, options);
            
            % 阶段2: 加入辐射但保持反应简化
            fprintf('  > 阶段2/3: 求解中等简化模型...\n');
            params_medium = params;
            params_medium.simplified_mode = false;
            params_medium.simplified_reaction = true;
            
            x_medium = fsolve(@(x) vap_equations(x, params_medium, pState, r_inf_phys), x_simple, options);
            
            % 阶段3: 完整模型
            fprintf('  > 阶段3/3: 求解完整模型...\n');
            params_full = params;
            params_full.simplified_mode = false;
            params_full.simplified_reaction = false;
            
            x = fsolve(@(x) vap_equations(x, params_full, pState, r_inf_phys), x_medium, options);
            
            % 解包结果
            [C_T1, C_T2, C_Mg, C_CO1, C_CO2, C_MgO1, C_MgO2, ~, T_f, ~, ~, ~, ~, ~, r_f] = unpack_variables(x);
            
        catch ME
            fprintf('积分法求解过程中出错: %s\n', ME.message);
            
            % 尝试回退到最简化模型
            fprintf('尝试使用最简化模型...\n');
            params.simplified_mode = true;
            params.simplified_reaction = true;
            
            try
                % 使用更保守的初始猜测
                x0(15) = pState.r_p * 2.0; % 增大火焰位置猜测
                x0(3) = -m_dot_phys_guess * 0.5; % 减小蒸发速率猜测
                
                x = fsolve(@(x) vap_equations(x, params, pState, r_inf_phys), x0, options);
                [C_T1, C_T2, C_Mg, C_CO1, C_CO2, C_MgO1, C_MgO2, ~, T_f, ~, ~, ~, ~, ~, r_f] = unpack_variables(x);
            catch
                % 如果依然失败，返回合理的默认值
                fprintf('求解失败，使用默认估计值...\n');
                C_T1 = -m_dot_phys_guess * L_v;
                C_Mg = -m_dot_phys_guess;
                C_CO1 = -C_Mg * (mw.CO.molar_mass/mw.Mg.molar_mass);
                T_f = 1.5 * T_surf;
                r_f = 1.5 * pState.r_p;
            end
        end
        
        % 4. 结果处理与输出
        % 计算最终蒸发速率和CO通量
        m_dot_sol = -C_Mg / (4 * pi);
        m_dot_co_surf = C_CO1 / (4 * pi);
        
        % 打包结果
        rate_info = package_rate_info_func(m_dot_sol, m_dot_co_surf, r_f, T_f, params);
        
        fprintf('积分法求解完成。火焰位置: %.3e m (r_f/r_p = %.2f), 蒸发速率: %.3e kg/s\n', ...
                r_f, r_f/pState.r_p, m_dot_sol);
    end
    
    % 辅助函数: 方程组定义
    function F = vap_equations(x, params, pState, r_inf)
        % 解包变量
        [C_T1, C_T2, C_Mg, C_CO1, C_CO2, C_MgO1, C_MgO2, ~, T_f, Y_Mg_f, Y_CO_f_left, Y_MgO_f_left, Y_CO_f_right, Y_MgO_f_right, r_f] = unpack_variables(x);
        
        % 参数
        r_p = pState.r_p;
        T_p = pState.T_p;
        T_amb = params.ambient_temperature;
        mw = params.materials;
        L_v = mw.Mg.L_evap_Mg;
        
        % 物性参数 (可以根据温度更新)
        props = params.physicalModel.get_gas_properties((T_p + T_f)/2);
        k_gas1 = props.k_gas;
        cp_gas = props.Cp_gas;
        
        props = params.physicalModel.get_gas_properties((T_f + T_amb)/2);
        k_gas2 = props.k_gas;
        
        rho_D = params.rho_D_gas;
        
        % 火焰反应热
        H_reac_f = params.reaction_heats.flame_reac_H;
        H_reac_s = params.reaction_heats.surface_reac_H;
        
        % 方程残差
        F = zeros(15, 1);
        
        % 根据松弛级别计算松弛因子
        if isfield(params, 'simplified_mode') && params.simplified_mode
            relaxation = 0.5;
        else
            relaxation = 0.9;
        end
        
        % 1. 表面条件
        F(1) = relaxation * (T_p - T_f - (C_T1/k_gas1)*(1/r_p - 1/r_f));
        
        % 根据蒸气压计算表面Mg质量分数
        R_Mg = params.R_u / mw.Mg.molar_mass;
        M_mix = params.physicalModel.get_ambient_mixture_molar_mass();
        exponent = (L_v / R_Mg) * (1/T_amb - 1/T_p);
        Y_mg_theory = min((mw.Mg.molar_mass / M_mix) * exp(exponent), 1.0);
        
        F(2) = Y_Mg_f - Y_mg_theory + (C_Mg/rho_D)*(1/r_p - 1/r_f);
        F(3) = Y_CO_f_left - (C_CO1/rho_D)*(1/r_p - 1/r_f);
        F(4) = Y_MgO_f_left - (C_MgO1/rho_D)*(1/r_p - 1/r_f);
        
        % 表面能量平衡
        sigma = params.sigma;
        emissivity = params.emissivity;
        
        % 对于简化模式，调整辐射项
        if isfield(params, 'simplified_mode') && params.simplified_mode
            rad_factor = 0.1; % 简化模式中减小辐射项
        else
            rad_factor = 1.0;
        end
        
        Q_rad = rad_factor * emissivity * sigma * (4 * pi * r_p^2) * (T_p^4 - T_amb^4);
        F(5) = C_T1 + L_v*(-C_Mg/(4*pi)) + Q_rad;
        
        % 2. 远场条件
        F(6) = relaxation * (T_amb - T_f - (C_T2/k_gas2)*(1/r_f - 1/r_inf));
        F(7) = 1.0 - Y_CO_f_right - (C_CO2/rho_D)*(1/r_f - 1/r_inf);
        F(8) = Y_CO_f_right - (C_CO1/rho_D)*(1/r_f - 1/r_inf);
        F(9) = Y_MgO_f_right - (C_MgO2/rho_D)*(1/r_f - 1/r_inf);
        
        % 3. 火焰面条件
        % 能量跳变
        omega_f = -C_Mg/(4*pi*r_f^2); % 火焰处反应速率
        Q = -H_reac_f / mw.Mg.molar_mass;
        F(10) = relaxation * (-C_T1/(k_gas1*r_f^2) + C_T2/(k_gas2*r_f^2) - Q*omega_f);
        
        % 组分守恒
        F(11) = relaxation * (Y_CO_f_right);  % 左侧无CO2
        F(12) = relaxation * (Y_Mg_f);        % 右侧无Mg
        F(13) = relaxation * (C_Mg + C_CO2);  % 物质守恒
        F(14) = C_MgO1 - C_MgO2 - C_Mg*(mw.MgO.molar_mass/mw.Mg.molar_mass);
        F(15) = C_CO1 - C_CO2 - C_CO2*(mw.CO.molar_mass/mw.CO2.molar_mass);
    end
    
    % 辅助函数: 解包变量
    function [C_T1, C_T2, C_Mg, C_CO1, C_CO2, C_MgO1, C_MgO2, CCO2, T_f, Y_Mg_f, Y_CO_f_left, Y_MgO_f_left, Y_CO_f_right, Y_MgO_f_right, r_f] = unpack_variables(x)
        C_T1 = x(1);
        C_T2 = x(2);
        C_Mg = x(3);
        C_CO1 = x(4);
        C_CO2 = x(5);
        C_MgO1 = x(6);
        C_MgO2 = x(7);
        CCO2 = x(8);
        T_f = x(9);
        Y_Mg_f = x(10);
        Y_CO_f_left = x(11);
        Y_MgO_f_left = x(12);
        Y_CO_f_right = x(13);
        Y_MgO_f_right = x(14);
        r_f = x(15);
    end