classdef VaporizationStage < handle
    properties
        params          % Simulation parameters
        physicalModel   % 物理模型计算器
        thermo          % 热力学数据读取器
    end
    
    methods
        function obj = VaporizationStage(params, physicalModel)
            % 构造函数
            if nargin > 0
            obj.params = params;
            obj.physicalModel = physicalModel;
                % 使用从上层传递下来的共享ThermoReader实例
                obj.thermo = physicalModel.thermo_reader; 
            end
        end
        
        function rate_info = solve(obj, particleState)
            % SOLVE: 为给定的颗粒状态求解准稳态BVP问题,并返回关键速率
            % 不再负责时间推进，仅作为ode45的"速率计算函数"的核心
            
            % 1. 求解准稳态气相场
            sol_bvp = obj.solve_radial_bvp(particleState);

            % 2. 从解中提取并返回关键结果
            % params = [m_dot_mg_reac, r_f, beta]
            rate_info.m_dot_mg_reac = sol_bvp.parameters(1);
            rate_info.r_f = sol_bvp.parameters(2);
            rate_info.beta = sol_bvp.parameters(3);
        end
        
        function sol = solve_radial_bvp(obj, pState)
            % --- 步骤 1: 定义物理常数和边界 ---
            r_p = pState.r_p;
            T_p = pState.T_p;
            r_inf = r_p + obj.params.flam_thickness;
            
            % 从参数文件加载摩尔质量 (kg/mol)
            M_Mg = obj.params.materials.Mg.molar_mass;
            M_CO2 = obj.params.materials.CO2.molar_mass;
            M_MgO = obj.params.materials.MgO.molar_mass;
            M_CO = obj.params.materials.CO.molar_mass;
            
            % 反应化学计量系数 (质量基)
            w_Mg = 1; % 参考
            w_CO2 = M_CO2 / M_Mg;
            w_MgO = M_MgO / M_Mg;
            w_CO = M_CO / M_Mg;
            
            % 从热力学数据文件加载生成焓 (J/mol)
            Hf_Mg = obj.thermo.get_property('Mg', 'Hf298');
            Hf_CO2 = obj.thermo.get_property('CO2', 'Hf298');
            Hf_MgO = obj.thermo.get_property('MgO', 'Hf298');
            Hf_CO = obj.thermo.get_property('CO', 'Hf298');
            
            % 计算火焰面反应热 (J/kg of Mg)
            Q_reac_f = (Hf_MgO + Hf_CO - Hf_Mg - Hf_CO2) / M_Mg;

            % --- 步骤 2: 设置BVP求解器 ---
            % 定义新的未知参数: [Mg蒸发率, 火焰半径, 产物内流比例]
            m_dot_mg_guess = 4e-7;    % kg/s
            r_f_guess = r_p + 0.5 * (r_inf - r_p); 
            beta_guess = 0.5;         % 0-1
            params_guess = [m_dot_mg_guess; r_f_guess; beta_guess];

            % 定义求解网格 (分段BVP需要元胞数组)
            mesh = {linspace(r_p, r_f_guess, 40), linspace(r_f_guess, r_inf, 40)};
            solinit = bvpinit(mesh, @initial_guess, params_guess);
            
            bvp_options = bvpset('Stats','on', 'RelTol', 1e-4, 'AbsTol', 1e-6, 'NMax', 5000);

            % 调用 bvp4c 求解
            sol = bvp4c(@bvp_ode, @bvp_bc, solinit, bvp_options);

            % --- 嵌套函数定义 ---

            function dzdr = bvp_ode(r, z, region, params)
                % 一阶ODE系统，在整个求解域内形式统一
                % z(1-5): T, Y_Mg, Y_MgO, Y_CO, Y_CO2
                % z(6-10): dTdr, dY_Mg_dr, dY_MgO_dr, dY_CO_dr, dY_CO2_dr
                
                % 解包基本未知参数
                m_dot_mg_reac = params(1);
                beta          = params(3);

                % 计算衍生的总产物生成率 (kg/s)
                m_prod_mgo_total = m_dot_mg_reac * w_MgO;
                m_prod_co_total  = m_dot_mg_reac * w_CO;

                % 计算在空间中恒定的总对流质量流率 (kg/s)
                % m_dot_total = (Mg蒸发) - (向内流动的产物)
                m_dot_total = m_dot_mg_reac - beta * (m_prod_mgo_total + m_prod_co_total);
                
                % 获取温度依赖的物性
                T = z(1);
                props = obj.physicalModel.get_gas_properties(T);
                k_gas = props.k_gas;
                cp_gas = props.cp_gas;
                rho_D = 0.0005; % 简化假设 rho*D, TODO: 使用更精确的模型

                % 求解所有10个微分方程
                dzdr = zeros(10, 1);
                dzdr(1:5) = z(6:10); % z' = z'
                
                % 能量方程: d(dT/dr)/dr = ...
                dzdr(6) = (m_dot_total * cp_gas / (4 * pi * k_gas)) * (z(6)/r^2) - (2/r) * z(6);
                
                % 组分方程: d(dY_i/dr)/dr = ... (i = Mg, MgO, CO, CO2)
                common_term = (m_dot_total / (4 * pi * rho_D));
                dzdr(7)  = common_term * (z(7)/r^2)  - (2/r) * z(7);
                dzdr(8)  = common_term * (z(8)/r^2)  - (2/r) * z(8);
                dzdr(9)  = common_term * (z(9)/r^2)  - (2/r) * z(9);
                dzdr(10) = common_term * (z(10)/r^2) - (2/r) * z(10);
            end

            function res = bvp_bc(yleft, yright, params)
                % 为三点BVP定义的23个边界条件
                
                % --- A. 解包所有变量 ---
                % 1. 解包未知参数
                m_dot_mg_reac = params(1);
                r_f           = params(2);
                beta          = params(3);
                
                % 2. 解包边界点状态向量
                T_p_sol = yleft(1,1); Y_p_sol = yleft(2:5,1); d_p_sol = yleft(6:10,1);
                T_f_L = yright(1,1); Y_f_L = yright(2:5,1); d_f_L = yright(6:10,1);
                T_f_R = yleft(1,2);  Y_f_R = yleft(2:5,2);  d_f_R = yleft(6:10,2);
                T_inf_sol = yright(1,2); Y_inf_sol = yright(2:5,2);

                % 3. 获取边界点的物性
                props_p = obj.physicalModel.get_gas_properties(T_p_sol);
                k_p = props_p.k_gas; rho_D_p = props_p.rho_gas * obj.params.D_gas;
                props_f = obj.physicalModel.get_gas_properties(T_f_L);
                k_f = props_f.k_gas; rho_D_f = props_f.rho_gas * obj.params.D_gas;

                % 4. 计算衍生的质量流率
                m_prod_mgo_total = m_dot_mg_reac * w_MgO;
                m_prod_co_total  = m_dot_mg_reac * w_CO;
                m_dot_total = m_dot_mg_reac - beta * (m_prod_mgo_total + m_prod_co_total);

                % --- B. 定义通量计算辅助函数 ---
                flux = @(r, Y, dY, rho_D_local) (Y * m_dot_total / (4*pi*r^2) - rho_D_local .* dY);
                
                % --- C. 定义23个边界条件的残差 ---
                res = zeros(23, 1);

                % --- BCs at r = r_p (颗粒表面, 6个条件) ---
                res(1) = T_p_sol - T_p; % 1. 温度
                
                % 确保 M_species 是列向量以匹配 Y_p_sol 的维度
                M_species = [M_Mg; M_MgO; M_CO; M_CO2];
                MW_mix_p = 1 / sum(Y_p_sol ./ M_species);
                p_sat_Mg = obj.thermo.get_p_sat(T_p_sol);
                Y_Mg_sat = (p_sat_Mg / obj.params.ambient_pressure) * (M_Mg / MW_mix_p);
                res(2) = Y_p_sol(1) - Y_Mg_sat; % 2. Mg饱和蒸气压
                
                res(4) = Y_p_sol(2); % 4. MgO完美沉积 (Y_MgO = 0)
                res(5) = Y_p_sol(3); % 5. CO完美反应 (Y_CO = 0)
                
                J_CO2_p = flux(r_p, Y_p_sol(4), d_p_sol(5), rho_D_p);
                res(6) = J_CO2_p; % 6. CO2在表面为惰性 (通量为0)

                % 3. 表面能量守恒
                H_dep_MgO = 1.4e7; H_reac_CO_surf = 6.2e6; % TODO: move to params
                J_MgO_p = flux(r_p, Y_p_sol(2), d_p_sol(3), rho_D_p);
                J_CO_p  = flux(r_p, Y_p_sol(3), d_p_sol(4), rho_D_p);
                J_Mg_p = m_dot_mg_reac / (4*pi*r_p^2); % 定义Mg蒸发通量

                q_cond_out = -k_p * d_p_sol(1);
                q_rad_out = obj.params.emissivity * obj.params.sigma * (T_p_sol^4 - obj.params.ambient_temperature^4);
                q_evap_out = J_Mg_p * obj.params.materials.Mg.L_evap_Mg;
                q_dep_MgO_in = -J_MgO_p * H_dep_MgO;
                q_reac_CO_in = -J_CO_p * H_reac_CO_surf;
                res(3) = q_cond_out + q_rad_out + q_evap_out - q_dep_MgO_in - q_reac_CO_in;

                % --- BCs at r = r_inf (无穷远, 5个条件) ---
                res(7) = T_inf_sol - obj.params.ambient_temperature; % T = T_inf
                res(8) = Y_inf_sol(1);    % Y_Mg = 0
                res(9) = Y_inf_sol(2);    % Y_MgO = 0
                res(10) = Y_inf_sol(3);   % Y_CO = 0
                res(11) = Y_inf_sol(4) - 1.0; % Y_CO2 = 1

                % --- BCs at r = r_f (火焰面, 12个条件) ---
                % 连续性条件
                res(12) = T_f_L - T_f_R;         % 温度连续
                res(15) = Y_f_L(2) - Y_f_R(2);   % Y_MgO连续
                res(16) = Y_f_L(3) - Y_f_R(3);   % Y_CO连续
                
                % 火焰面定义 (反应物耗尽)
                res(13) = Y_f_L(1); % 左侧Mg耗尽 (Y_Mg(r_f-) = 0)
                res(14) = Y_f_R(4); % 右侧CO2耗尽 (Y_CO2(r_f+) = 0)
                
                % 计算火焰面处的通量
                J_E_L = -k_f * d_f_L(1); J_E_R = -k_f * d_f_R(1);
                J_Mg_L = flux(r_f, Y_f_L(1), d_f_L(2), rho_D_f);
                J_MgO_L = flux(r_f, Y_f_L(2), d_f_L(3), rho_D_f); J_MgO_R = flux(r_f, Y_f_R(2), d_f_R(3), rho_D_f);
                J_CO_L = flux(r_f, Y_f_L(3), d_f_L(4), rho_D_f); J_CO_R = flux(r_f, Y_f_R(3), d_f_R(4), rho_D_f);
                J_CO2_L = flux(r_f, Y_f_L(4), d_f_L(5), rho_D_f);
                J_CO2_R = flux(r_f, Y_f_R(4), d_f_R(5), rho_D_f);
                
                % 通量跳跃/守恒条件
                source_E_f = (m_dot_mg_reac * Q_reac_f) / (4*pi*r_f^2);
                res(17) = (J_E_L - J_E_R) - source_E_f; % 能量通量跳跃
                
                source_MgO_f = m_prod_mgo_total / (4*pi*r_f^2);
                res(19) = (J_MgO_R - J_MgO_L) - source_MgO_f; % MgO通量跳跃
                
                source_CO_f = m_prod_co_total / (4*pi*r_f^2);
                res(20) = (J_CO_R - J_CO_L) - source_CO_f; % CO通量跳跃
                
                res(21) = J_CO2_L; % 火焰面左侧CO2通量为0
                
                % 定义产物流向的条件
                res(22) = J_MgO_L - (-beta * source_MgO_f); % 流向内部的MgO通量
                res(23) = J_CO_L - (-beta * source_CO_f);  % 流向内部的CO通量
                
                % 化学计量比条件 (替换了旧的独立通量条件)
                % 此处需要两个独立的残差来约束两个独立的通量，但我们只有一个res(18)了。
                % 最核心的约束是化学计量比。
                m_reac_co2_total = m_dot_mg_reac * w_CO2;
                
                % 定义一个组合残差，强制两个通量与各自的源项的化学计量比相等
                err_mg = J_Mg_L / (-m_dot_mg_reac / (4*pi*r_f^2)); % 应为1
                err_co2 = J_CO2_R / (-m_reac_co2_total / (4*pi*r_f^2)); % 应为1
                res(18) = err_mg - err_co2; % 强制比率相等
                end
                
            function y = initial_guess(r, k, p)
                % 为BVP求解器提供状态向量的初始猜测
                % y(1-5): T, Y_Mg, Y_MgO, Y_CO, Y_CO2
                % y(6-10): 各自的导数
                
                r_f_g = p(2); % r_f is now the second parameter
                
                y = zeros(10,1);
                
                % 温度猜测: 从T_p线性下降到T_amb
                T_f_guess = 2000;
                
                % 根据区域索引 k 提供猜测
                if k == 1 % 区域 1: r_p -> r_f
                    y(1) = T_p - (T_p - T_f_guess) * (r - r_p) / (r_f_g - r_p);
                    y(2) = 0.1 * (1 - (r - r_p) / (r_f_g - r_p)); % Y_Mg
                    y(3) = 0.05 * ((r - r_p) / (r_f_g - r_p));    % Y_MgO
                    y(4) = 0.05 * ((r - r_p) / (r_f_g - r_p));    % Y_CO
                    y(5) = 0; % Y_CO2
                else % 区域 2: r_f -> r_inf
                    y(1) = T_f_guess - (T_f_guess - obj.params.ambient_temperature) * (r - r_f_g) / (r_inf - r_f_g);
                    y(2) = 0; % Y_Mg
                    y(3) = 0.1 * (1 - (r - r_f_g) / (r_inf - r_f_g)); % Y_MgO
                    y(4) = 0.1 * (1 - (r - r_f_g) / (r_inf - r_f_g)); % Y_CO
                    y(5) = 1.0 * (r - r_f_g) / (r_inf - r_f_g); % Y_CO2
                end
                
                % 导数猜测可以设为0
            end
        end
    end
end