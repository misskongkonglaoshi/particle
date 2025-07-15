classdef VaporizationStage < handle
    % VaporizationStage (气化阶段求解器) - V5 (最终版: 统一BVP求解)
    %
    % 核心功能:
    % 对于一个处于沸点的、给定状态的颗粒，通过一次BVP调用，求解一个在物理上
    % 完全自洽的稳态燃烧模型，同时得到反应速率和火焰结构。
    %
    % 求解逻辑:
    % 使用bvp4c求解一个包含未知参数的边界值问题。
    % 未知数 = [10个状态变量(T, Yi, d/dr), 火焰半径(r_f), 总蒸发速率(m_dot)]
    % 方程数 = [10个ODE, 12个边界/界面条件]
    % 通过这种方式，所有的物理定律（能量守恒、质量守恒、化学计量）都在
    % 一个统一的框架内被同时满足。

    properties
        params          % 参数对象
        physicalModel   % 物理模型对象
        thermo          % 热力学数据读取器
    end

    methods
        function obj = VaporizationStage(params, physicalModel)
            % 构造函数
            obj.params = params;
            obj.physicalModel = physicalModel;
            obj.thermo = physicalModel.thermo_reader;
        end

        function rate_info = solve_reaction_rates(obj, pState)
            % --- 主求解函数 ---
            % 直接调用统一的BVP求解器
            try
                rate_info = obj.solve_unified_bvp(pState);
            catch ME
                fprintf('统一BVP求解失败: %s\n', ME.message);
                
                % V12 调试: 在BVP失败时输出关键尺度参数
                r_p = pState.r_p;
                r_inf = r_p + obj.params.flam_thickness;
                fprintf('  > 失败时关键参数: r_p = %.3e m, r_inf = %.3e m (尺度差异: %.1f 倍)\n', ...
                        r_p, r_inf, r_inf/r_p);

                % 在某些极端条件下(如颗粒过小)，BVP可能无解。
                % 在这种情况下，返回零速率，让ODE求解器终止。
                rate_info.dmdt_mg = 0;
                rate_info.dmdt_mgo = 0;
                rate_info.dmdt_c = 0;
                rate_info.r_f = pState.r_p;
                rate_info.T_f = pState.T_p;
            end
        end

        function rate_info = solve_unified_bvp(obj, pState)
            % --- BVP求解器 (求解所有未知量) ---
            
            r_p = pState.r_p;
            r_inf = r_p + obj.params.flam_thickness; % V11 恢复: 使用固定厚度

            % --- 对未知参数提供合理的初始猜测 ---
            % p(1) = r_f, p(2) = m_dot_total
            r_f_guess = r_p * 1.5;
            m_dot_guess = 4 * pi * r_p * 1e-4;

            % --- 准备BVP的初始结构 (回到多区域方案A) ---
            mesh1 = linspace(r_p, r_f_guess, 15);
            mesh2 = linspace(r_f_guess, r_inf, 15);

            % 多区域模式下 initial_guess 需要接收 region 参数
            y_init_fun = @(r, region) obj.initial_guess(r, region, pState, r_f_guess);
            solinit    = bvpinit({mesh1, mesh2}, y_init_fun, [r_f_guess, m_dot_guess]);

            % --- 创建句柄并调用求解器 (多区域) ---
            ode_handle = @(r, z, region, p) obj.bvp_ode(r, z, region, p(2));
            bc_handle = @(yl, yr, p) obj.bvp_bc(yl, yr, p, pState);
            
            options = bvpset('NMax', 6000, 'RelTol', 1e-4, 'AbsTol', 1e-7);
            sol = bvp4c(ode_handle, bc_handle, solinit, options);
            
            % --- 从解中提取结果并打包 ---
            r_f_sol = sol.parameters(1);
            m_dot_sol = sol.parameters(2);
            
            % 计算表面的CO消耗速率以完成打包
            y_surf = sol.y(:,1);
            dYco_dr_surf = y_surf(9);
            Yco_surf = y_surf(4);
            rho_D_avg = obj.params.rho_D_gas;
            A_p = 4 * pi * pState.r_p^2;
            m_dot_co_surf = -(A_p * rho_D_avg) * dYco_dr_surf + m_dot_sol * Yco_surf;

            rate_info = obj.package_rate_info(m_dot_sol, m_dot_co_surf, sol, r_f_sol);
        end
        
        function dzdr = bvp_ode(obj, r, z, region, m_dot) % 修正: 明确使用region参数
            % BVP的常微分方程 (稳态输运方程)
            % V16 修正: 明确使用 region 参数，这对多区域 bvp4c 至关重要
            
            % 提取当前状态变量
            T = z(1,:);
            props = obj.physicalModel.get_gas_properties(T);
            k_gas = props.k_gas;
            cp_gas = props.Cp_gas;
            rho_D = obj.params.rho_D_gas;
            
            % 初始化导数矩阵
            dzdr = zeros(size(z));
            
            % 根据区域应用相应的方程
            if region == 1  % 区域1: r_p -> r_f
                % 导数关系
                dzdr(1:5,:) = z(6:10,:);
                
                % 能量方程
                dzdr(6,:) = (m_dot * cp_gas ./ (4 * pi * r.^2 .* k_gas)) .* z(6,:) - (2./r) .* z(6,:);
                
                % 组分方程
                common_term = m_dot ./ (4 * pi * r.^2 * rho_D);
                dzdr(7:10,:) = common_term .* z(7:10,:) - (2./r) .* z(7:10,:);
            else  % 区域2: r_f -> r_inf
                % 导数关系
                dzdr(1:5,:) = z(6:10,:);
                
                % 能量方程
                dzdr(6,:) = (m_dot * cp_gas ./ (4 * pi * r.^2 .* k_gas)) .* z(6,:) - (2./r) .* z(6,:);
                
                % 组分方程
                common_term = m_dot ./ (4 * pi * r.^2 * rho_D);
                dzdr(7:10,:) = common_term .* z(7:10,:) - (2./r) .* z(7:10,:);
            end
        end

        function res = bvp_bc(obj, yleft, yright, p, pState)
            % BVP的边界/界面条件 (12个条件 for 10 ODEs + 2 parameters)
            r_f = p(1);
            m_dot = p(2);
            r_p = pState.r_p;

            % --- 提取边界值 (多区域标准) ---
            % yleft / yright 维度: nState × nRegion
            % region 1 : r_p → r_f ; region 2 : r_f → r_inf
            y_at_rp       = yleft(:,1);   % r = r_p
            y_at_rf_left  = yright(:,1);  % r → r_f  (左侧)
            y_at_rf_right = yleft(:,2);   % r_f ← r  (右侧)
            y_at_rinf     = yright(:,2);  % r = r_inf

            r_inf = pState.r_p + obj.params.flam_thickness;
            rho_D_avg = obj.params.rho_D_gas;
            res = zeros(12,1);

            %% 1-4 颗粒表面 r = r_p
            res(1) = y_at_rp(1) - pState.T_p;
            B_m = exp(m_dot / (4*pi*r_p*rho_D_avg)) - 1;
            Y_mg_s = B_m/(1+B_m);
            res(2) = y_at_rp(2) - Y_mg_s;
            res(3) = y_at_rp(8);
            T_surf = y_at_rp(1);
            dTdr_surf  = y_at_rp(6);
            dYco_dr_s  = y_at_rp(9);
            Yco_surf   = y_at_rp(4);
            k_gas_surf = obj.physicalModel.get_gas_properties(T_surf).k_gas;
            A_p        = 4*pi*r_p^2;
            Q_cond     = -k_gas_surf * dTdr_surf * A_p;
            m_dot_co_s = -(A_p*rho_D_avg)*dYco_dr_s + m_dot*Yco_surf;
            H_reac_s   = obj.thermo.get_reaction_enthalpy('Mg(g)+CO(g)=MgO(s)+C(s)',T_surf);
            Q_reac_s   = m_dot_co_s * (-H_reac_s / obj.params.materials.CO.molar_mass);
            Q_rad      = obj.params.emissivity*obj.params.sigma*(T_surf^4-obj.params.ambient_temperature^4)*A_p;
            Q_evap     = m_dot*obj.params.materials.Mg.latent_heat;
            res(4)     = (Q_cond + Q_reac_s) - (Q_rad + Q_evap);

            %% 5-7 远场 r = r_inf
            res(5) = y_at_rinf(1) - obj.params.ambient_temperature;
            res(6) = y_at_rinf(3) - 1.0;
            res(7) = y_at_rinf(4);

            %% 8-12 火焰界面 r = r_f
            res(8)  = y_at_rf_left(1) - y_at_rf_right(1);
            res(9)  = y_at_rf_left(2);
            res(10) = y_at_rf_right(3);
            mw_mg  = obj.params.materials.Mg.molar_mass;
            mw_co2 = obj.params.materials.CO2.molar_mass;
            J_mg  = -(4*pi*r_f^2*rho_D_avg)*y_at_rf_left(7) + m_dot*y_at_rf_left(2);
            J_co2 = -(4*pi*r_f^2*rho_D_avg)*y_at_rf_right(8)+ m_dot*y_at_rf_right(3);
            res(11)= J_mg/mw_mg + J_co2/mw_co2;
            T_f    = y_at_rf_left(1);
            H_reac_f = obj.thermo.get_reaction_enthalpy('Mg(g)+CO2(g)=MgO(g)+CO(g)',T_f);
            m_dot_mg_f = J_mg;
            Q_flame = m_dot_mg_f*(-H_reac_f/mw_mg);
            k_gas_f = obj.physicalModel.get_gas_properties(T_f).k_gas;
            res(12)= (y_at_rf_right(6)-y_at_rf_left(6))*k_gas_f + Q_flame/(4*pi*r_f^2);
        end
                
        function y = initial_guess(obj, r, region, pState_ig, r_f_ig) %#ok<INUSL>
            % BVP的初始猜测值 (V15修正: 确保为网格点返回矩阵)
            y = zeros(10, length(r)); % y的大小应为 (状态数 x 网格点数)
            T_f_guess = 2500;
            T_p = pState_ig.T_p;
            T_amb = obj.params.ambient_temperature;
            r_p = pState_ig.r_p;
            r_inf = r_p + obj.params.flam_thickness; 
            Y_mg_s_guess = 0.1; % Mg表面质量分数的猜测值

            if region == 1 % 区域1: r_p -> r_f
                frac = (r - r_p) / (r_f_ig - r_p + eps); % 加eps避免除零
                y(1,:) = T_p + (T_f_guess - T_p) .* frac;
                y(2,:) = Y_mg_s_guess .* (1 - frac);
                y(3,:) = 0; % MATLAB会自动广播
                y(4,:) = 0.1 .* frac .* (1-frac); % CO
                y(5,:) = 0.1 .* frac; % MgO
            else % 区域2: r_f -> r_inf
                frac = (r - r_f_ig) / (r_inf - r_f_ig + eps); % 加eps避免除零
                y(1,:) = T_f_guess - (T_f_guess - T_amb) .* frac;
                y(2,:) = 0;
                y(3,:) = 1.0 .* frac;
                y(4,:) = 0.1 .* (1-frac);
                y(5,:) = 0.1; % MATLAB会自动广播
            end
        end
        
        function rate_info_out = package_rate_info(obj, m_dot_mg_total, m_dot_co_surf, sol, r_f)
            % --- 辅助函数: 将最终结果打包成结构体 ---
            mw_mg = obj.params.materials.Mg.molar_mass;
            mw_mgo = obj.params.materials.MgO.molar_mass;
            mw_c = obj.params.materials.C.molar_mass;
            mw_co = obj.params.materials.CO.molar_mass;

            rate_info_out.dmdt_mg = -m_dot_mg_total;
            
            % 计算表面生成的MgO和C的速率
            m_dot_mgo_surf = m_dot_co_surf * (mw_mgo / mw_co);
            rate_info_out.dmdt_c = m_dot_co_surf * (mw_c / mw_co);
            
            % 计算火焰面生成的MgO的速率
            m_dot_mg_surf_equiv = m_dot_co_surf * (mw_mg/mw_co);
            m_dot_mg_flame = m_dot_mg_total - m_dot_mg_surf_equiv;
            m_dot_mgo_flame = m_dot_mg_flame * (mw_mgo / mw_mg);

            rate_info_out.dmdt_mgo = m_dot_mgo_surf + m_dot_mgo_flame;
            
            % 存储火焰面信息
            rate_info_out.r_f = r_f;
            sol_at_rf = deval(sol, r_f);
            rate_info_out.T_f = sol_at_rf(1);
        end
    end
end