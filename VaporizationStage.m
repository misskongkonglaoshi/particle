classdef VaporizationStage < handle
    % VaporizationStage 气相燃烧阶段求解器 (拟稳态BVP方法)

    properties
        params          % Simulation parameters
        physicalModel   % 物理模型计算器
        
        % Data logging
        results         % Struct to store time-history of the solution
    end
    
    methods
        function obj = VaporizationStage(params, physicalModel)
            % 构造函数
            obj.params = params;
            obj.physicalModel = physicalModel;
            
            % Initialize results struct
            obj.results.time = [];
            obj.results.particleStateHistory = [];
        end
        
        function [t_history, state_history] = solve(obj, particleState, t_start)
            % SOLVE: 启动求解器
            % 输入:
            %   particleState: 当前的颗粒状态对象
            %   t_start: 阶段开始时间

            % 1. 设置ODE求解器参数
            t_end = t_start + obj.params.t_combustion;
            tspan = [t_start, t_end];

            % 初始状态向量 y0 = [颗粒温度; 镁核质量; 氧化镁壳质量]
            y0 = [particleState.T_p; particleState.m_mg; particleState.m_mgo];
            
            % 设置ode45选项, 控制精度并确保质量不为负
            options = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, 'NonNegative', [2,3]);

            % 2. 调用 ode45 求解
            fprintf('开始求解气相燃烧阶段...\n');
            
            % 创建一个临时的ParticleState对象用于ODE内部计算，避免污染原始对象
            transientState = particleState.copy();
            
            [t, y] = ode45(@(t, y) particle_ode_system(t, y, transientState), tspan, y0, options);
            fprintf('气相燃烧阶段求解完成。\n');

            % 3. 处理并存储结果
            t_history = t;
            num_steps = length(t);
            state_history(1, num_steps) = ParticleState(obj.params); % Pre-allocate
            
            for i = 1:num_steps
                % 从历史记录中填充每个时间点的颗粒状态
                tempState = ParticleState(obj.params);
                tempState.T_p = y(i, 1);
                tempState.m_mg = y(i, 2);
                tempState.m_mgo = y(i, 3);
                tempState.update_geometry_from_mass(); % 根据质量更新半径
                state_history(i) = tempState;
            end
            
            % 4. 将最终状态写回传入的particleState对象
            if ~isempty(y)
                particleState.T_p = y(end, 1);
                particleState.m_mg = y(end, 2);
                particleState.m_mgo = y(end, 3);
                particleState.update_geometry_from_mass();
            end

            % --- 嵌套函数定义 ---
            function dydt = particle_ode_system(t, y, pState)
                % --- 步骤 A: 解包当前状态 ---
                % 从 y 向量更新瞬时颗粒状态
                pState.T_p   = y(1);
                pState.m_mg  = y(2);
                pState.m_mgo = y(3);
                pState.update_geometry_from_mass(); % 确保半径等参数是最新的
                
                % 如果颗粒燃尽，则停止计算
                if pState.m_mg < 1e-15
                    dydt = [0; 0; 0];
                    return;
                end

                % --- 步骤 B: 求解准稳态径向BVP ---
                sol_bvp = solve_radial_bvp(pState);

                % --- 步骤 C: 从BVP解计算表面通量 ---
                    % 使用 deval 计算颗粒表面的解的导数 (梯度)
                    z_prime_surface = deval(sol_bvp, pState.r_p, 1);
                    
                    % --- 用户需实现: 从梯度计算物理通量 ---
                    gas_props = obj.physicalModel.get_gas_properties(pState.T_p);
                    k_gas = gas_props.k_gas;

                    T_gas_gradient = z_prime_surface(1); % 假设 z(1) 是气体温度
                    heat_flux_conduction = -k_gas * T_gas_gradient; 
                    
                    % 质量通量 (蒸发/反应)
                    mass_flux_mg_evaporation = 0; % 占位符: 应由BVP解的组分梯度计算
                    mass_flux_mgo_formation = 0; % 占位符: 应由BVP解的组分梯度计算
                    
                    % --- 步骤 D: 基于通量计算状态导数 (能量和质量平衡) ---
                    surface_area = 4 * pi * pState.r_p^2;
                    
                    % 能量平衡
                    heat_loss_radiation = obj.params.emissivity * obj.params.sigma * (pState.T_p^4 - obj.params.ambient_temperature^4);
                    Q_net = (heat_flux_conduction - heat_loss_radiation) * surface_area;
                    Cp_total = obj.physicalModel.calculate_total_heat_capacity(pState);
                    
                    if Cp_total > 1e-9
                        dTp_dt = Q_net / Cp_total;
                    else
                        dTp_dt = 0;
                    end

                    % 质量平衡
                    dm_mg_dt = -mass_flux_mg_evaporation * surface_area; 
                    dm_mgo_dt = mass_flux_mgo_formation * surface_area;

                dydt = [dTp_dt; dm_mg_dt; dm_mgo_dt];
            end

            function sol = solve_radial_bvp(pState)
                r_p = pState.r_p;
                T_p = pState.T_p;

                % --- 步骤 1: 定义边界并为未知参数提供初始猜测 ---
                r_inf = r_p + obj.params.flam_thickness;
                
                % 为未知参数提供初始猜测值
                r_f_guess = r_p + 0.4 * (r_inf - r_p); 
                m_dot_guess = 1e-5; % kg/s, 总质量流率的初始猜测
                
                % 确保火焰面猜测值在有效范围内
                if r_f_guess <= r_p + 1e-6 || r_f_guess >= r_inf - 1e-6
                    error('火焰面初始猜测 r_f = %.4e 无效', r_f_guess);
                end

                % --- 步骤 2: 为多点BVP设置求解器 (8变量系统) ---
                % 变量顺序: z = [T, Y_Mg, Y_MgO, Y_CO, Y_CO2, Y_O2, m_dot, r_f]^T
                solinit = bvpinit([r_p, r_f_guess, r_inf], ...
                                  @(r,k) initial_guess(r, k, r_p, r_inf, T_p, m_dot_guess));
                
                bvp_options = bvpset('Stats','off', 'RelTol', 1e-4, 'NMax', 5000);

                sol = bvp4c(@(r,z,region) bvp_ode(r, z, region), ...
                            @(yleft,yright) bvp_bc(yleft, yright, T_p), ...
                            solinit, bvp_options);

                % --- 嵌套函数定义 (已适配8变量系统) ---

                function dzdr = bvp_ode(r, z, region)
                    % --- 径向输运方程 ---
                    % 解包变量
                    T     = z(1);%火焰温度
                    Y_Mg  = z(2);
                    Y_MgO = z(3);
                    Y_CO  = z(4);
                    Y_CO2 = z(5);
                    Y_O2  = z(6);
                    m_dot = z(7); % 总质量流率 (常数)
                    r_f   = z(8); % 火焰半径 (常数)

                    % 获取物性 (占位符)
                    gas_props = obj.physicalModel.get_gas_properties(T);
                    k_gas = gas_props.k_gas;
                    cp_gas = gas_props.cp_gas;
                    rho_gas = gas_props.rho_gas;
                    D_species = 1e-5; % 扩散系数 (占位符)
                    
                    % 初始化导数和源项
                    dTdr = 0;
                    dYdr = zeros(5, 1);
                    omega_T = 0; % 能量源项 (J/m^3/s)
                    omega_Y = zeros(5, 1); % 组分源项 (kg/m^3/s)

                    % 分区域的化学反应源项 (占位符)
                    switch region
                        case 1 % 区域 1: r_p -> r_f (Mg + CO2 -> MgO + CO)
                            % omega_Y(Mg) = ...
                            % omega_Y(CO2) = ...
                            % omega_Y(MgO) = ...
                            % omega_Y(CO) = ...
                        case 2 % 区域 2: r_f -> r_inf (CO + 0.5 O2 -> CO2)
                            % omega_Y(CO) = ...
                            % omega_Y(O2) = ...
                            % omega_Y(CO2) = ...
                    end

                    % 输运方程 (球坐标, 拟稳态)
                    % 能量方程: d/dr(r^2*k*dT/dr) - (m_dot*cp/4pi)*dT/dr + r^2*omega_T = 0
                    % 此为二阶方程，bvp4c需要一阶形式。需引入 q = r^2*k*dT/dr
                    % 这里为简化，直接写出一阶形式的最终结构
                    % (实际实现时需要转换或使用不同的数值方法)
                    
                    % --- 简化版占位符方程 ---
                    convection_coeff = m_dot * cp_gas / (4 * pi * r^2);
                    diffusion_coeff_T = k_gas;
                    dTdr = (convection_coeff * T_something - r^2*omega_T) / diffusion_coeff_T; % 这是一个非常简化的例子
                    
                    % 组分方程:
                    % dYdr(i) = ...
                    
                    % 返回导数向量
                    dzdr_phys = [dTdr; dYdr];
                    dzdr = [dzdr_phys; 0; 0]; % m_dot 和 r_f 的导数为0
                end

                function res = bvp_bc(yleft, yright, T_p_bvp)
                    % --- 多点边界条件 (8变量系统 -> 16个条件) ---
                    % 变量顺序: z = [T, Y_Mg, Y_MgO, Y_CO, Y_CO2, Y_O2, m_dot, r_f]^T
                    res = zeros(16, 1);

                    % 从解向量中获取未知参数 (在所有点上都相同)
                    m_dot = yleft(7,1);
                    r_f   = yleft(8,1);
                    
                    % --- 物理常数 (占位符) ---
                    L_evap_Mg = 6.0e6; % J/kg, Mg的蒸发潜热

                    % --- 获取表面和火焰面的导数值 (用于通量/定义条件) ---
                    % 表面 (r=r_p), 状态 yleft(:,1)
                    z_p = yleft(:,1);
                    r_p = pState.r_p; % 从外部pState获取精确的r_p
                    dzdr_p = bvp_ode(r_p, z_p, 1);
                    dTdr_p = dzdr_p(1);

                    % 火焰面 (r=r_f), 状态 yright(:,1)
                    z_f = yright(:,1);
                    dzdr_f = bvp_ode(r_f, z_f, 1);
                    dTdr_f = dzdr_f(1);

                    % --- 边界条件 ---
                    % -- 边界 1: 颗粒表面 r = r_p (6个条件) --
                    res(1) = yleft(1,1) - T_p_bvp; % 1. 温度 = 颗粒温度
                    res(2) = yleft(2,1) - calculate_saturation_pressure(T_p_bvp); % 2. Mg质量分数 = 饱和值
                    res(3) = yleft(6,1); % 3. 占位符: 表面O2质量分数为0 (应为通量=0)
                    % 以下3个是更真实的通量平衡占位符，但需要物性
                    res(4) = m_dot * yleft(4,1); % 4. 占位符: CO表面通量为0 (m*Y - rho*D*dY/dr = 0)
                    res(5) = m_dot * yleft(5,1); % 5. 占位符: CO2表面通量为0
                    res(6) = (4*pi*r_p^2 * 1.0 * dTdr_p) - m_dot * L_evap_Mg; % 6. 能量平衡 -> 求解 m_dot (k=1.0是占位符)

                    % -- 边界 2: 火焰面 r = r_f (7个条件) --
                    res(7) = dTdr_f; % 7. 火焰定义 -> 求解 r_f (dT/dr = 0)
                    res(8)  = yright(1,1) - yleft(1,2); % 8. 温度连续
                    res(9)  = yright(2,1) - yleft(2,2); % 9. Y_Mg 连续 (物理上应为Y_Mg(rf-)=0和通量平衡)
                    res(10) = yright(3,1) - yleft(3,2); % 10. Y_MgO 连续
                    res(11) = yright(4,1) - yleft(4,2); % 11. Y_CO 连续
                    res(12) = yright(5,1) - yleft(5,2); % 12. Y_CO2 连续
                    res(13) = yright(6,1) - yleft(6,2); % 13. Y_O2 连续

                    % -- 边界 3: 远场 r = r_inf (3个条件) --
                    res(14) = yright(1,2) - obj.params.ambient_temperature; % 14. 温度 = 环境温度
                    res(15) = yright(6,2) - obj.params.ambient_gas_composition.O2; % 15. O2 = 环境值
                    res(16) = yright(5,2) - obj.params.ambient_gas_composition.CO2; % 16. CO2 = 环境值 (若环境中有)
                end
                
                function z_guess = initial_guess(r, k, r_p_bvp, r_inf_bvp, T_p_bvp, m_dot_guess_val)
                     % --- 分区域的初始猜测 (8变量系统) ---
                     
                     % 内部猜测
                     r_f_guess = r_p_bvp + 0.4 * (r_inf_bvp - r_p_bvp);
                     T_flame_guess = 3000;
                     Y_Mg_surf_guess = calculate_saturation_pressure(T_p_bvp);

                     % 初始化所有物理变量
                     phys_vars_guess = zeros(6, 1);
                     
                     switch k
                         case 1 % 区域 1: r_p -> r_f
                             p = (r - r_p_bvp) / (r_f_guess - r_p_bvp);
                             T_guess = T_p_bvp + (T_flame_guess - T_p_bvp) * 4 * (p - p^2);
                             Y_Mg_guess = Y_Mg_surf_guess * (1 - p);
                             Y_MgO_guess = 0.1 * p; % 线性增长
                             Y_CO_guess = 0.2 * p; % 线性增长
                             Y_CO2_guess = 0.1 * (1-p); % 线性消耗
                             Y_O2_guess = 0; % 假设区域1无氧气
                             phys_vars_guess = [T_guess; Y_Mg_guess; Y_MgO_guess; Y_CO_guess; Y_CO2_guess; Y_O2_guess];

                         case 2 % 区域 2: r_f -> r_inf
                             p = (r - r_f_guess) / (r_inf_bvp - r_f_guess);
                             T_guess = T_flame_guess + p * (obj.params.ambient_temperature - T_flame_guess);
                             Y_Mg_guess = 0;
                             Y_MgO_guess = 0.1 * (1-p); % 从火焰处开始衰减
                             Y_CO_guess = 0.2 * (1-p); % 从火焰处开始消耗
                             Y_CO2_guess = (0.2) + p * (obj.params.ambient_gas_composition.CO2 - 0.2); % 线性恢复
                             Y_O2_guess = p * obj.params.ambient_gas_composition.O2; % 从0恢复到环境值
                             phys_vars_guess = [T_guess; Y_Mg_guess; Y_MgO_guess; Y_CO_guess; Y_CO2_guess; Y_O2_guess];
                     end
                     
                     % 组合成包含未知参数的完整猜测向量
                     z_guess = [phys_vars_guess; m_dot_guess_val; r_f_guess];
                end
            end
        end
    end
end

function Y_sat = calculate_saturation_pressure(T)
    % Placeholder for calculating saturation pressure of Mg vapor.
    % This should be a proper physical correlation.
    Y_sat = 1e-3 * exp(-3000/T + 5.0); % Dummy correlation
end