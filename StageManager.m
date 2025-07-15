classdef StageManager < handle
    % StageManager 阶段管理器
    % 管理不同阶段的切换和求解 (动态调度器)
    
    properties
        model         % 颗粒模型 (ParticleModel)
        params        % 参数对象 (Parameters)
        physicalModel % 物理计算模型 (PhysicalModel)
        thermo_reader % 热力学数据读取器 (ThermoReader)
        
        % 结果存储
        results       % 包含所有历史记录的结构体
    end
    
    methods
        function obj = StageManager(params)
            % 构造函数：初始化阶段管理器
            if nargin > 0
                obj.params = params;
                
                % 初始化依赖的模型
                obj.thermo_reader = ThermoReader('thermo.dat', obj.params);
                obj.physicalModel = PhysicalModel(obj.params, obj.thermo_reader);
                obj.model = ParticleModel(obj.params);
                
                % 初始化结果存储
                obj.results.time = [];
                obj.results.particleStateHistory = [];
                obj.results.stageHistory = {};
            end
        end
        
        function results = solve(obj)
            % 求解整个过程
            
            % 清空历史记录
            obj.results.time = [];
            obj.results.particleStateHistory = [];
            obj.results.stageHistory = {};
            
            current_time = 0;
            
            % 循环求解各阶段
            while current_time < obj.params.total_time
                % 1. 确定当前阶段
                current_stage_name = obj.model.determineStage();
                
                % 2. 根据阶段动态创建求解器并求解
                fprintf('当前时间: %.4fs, 进入 %s 阶段...\n', current_time, current_stage_name);
                
                t_stage = [];
                state_stage_history = ParticleState.empty;
                
                switch current_stage_name
                    case {'preheating', 'liquid_heating'}
                        if strcmp(current_stage_name, 'preheating')
                            solver = PreheatingStage(obj.params, obj.physicalModel);
                            target_temp = obj.params.materials.Mg.melting_point;
                            fprintf('--- 进入固相加热阶段 ---\n');
                            fprintf('  (目标: 熔点 %.1f K)\n', target_temp);
                        else % liquid_heating
                            solver = HeatingStage(obj.params, obj.physicalModel);
                            target_temp = obj.params.materials.Mg.ignition_temp;
                            fprintf('--- 进入液相加热阶段 ---\n');
                            fprintf('  (目标: 点火温度 %.1f K)\n', target_temp);
                        end
                        
                        % 前置条件检查: 如果当前温度已达到或超过目标，则跳过求解
                        if obj.model.particleState.T_p >= target_temp
                            fprintf('  (跳过: 当前温度 %.1f K >= 目标温度 %.1f K)\n', obj.model.particleState.T_p, target_temp);
                            t_stage = [current_time, current_time];
                            T_stage = [obj.model.particleState.T_p, obj.model.particleState.T_p];
                        else
                            initial_state_template = obj.model.particleState.copy(); % 保存阶段初始状态
                            [t_stage, T_stage] = solver.solve(obj.model.particleState, current_time, target_temp);
                        end

                        % 直接根据HeatingStage的输出构建正确的状态历史
                        num_steps = length(T_stage);
                        if ~isempty(t_stage)
                            state_stage_history(1, num_steps) = ParticleState(); % 预分配
                            for i = 1:num_steps
                                tempState = initial_state_template.copy(); % 从初始状态模板复制
                                tempState.T_p = T_stage(i);                % 仅更新当前步的温度
                                state_stage_history(i) = tempState;
                            end
                        end

                    case 'melting'
                        solver = MeltingStage(obj.params, obj.physicalModel);

                        % --- 开始修改 (2/2) ---
                        initial_state_template = obj.model.particleState.copy(); % 保存阶段初始状态
                        [t_stage, X_stage] = solver.solve(obj.model.particleState, current_time);
                        
                        % 直接根据MeltingStage的输出构建正确的状态历史
                        num_steps = length(X_stage);
                        if ~isempty(t_stage)
                            state_stage_history(1, num_steps) = ParticleState(); % 预分配
                            for i = 1:num_steps
                                tempState = initial_state_template.copy();      % 从初始状态模板复制
                                tempState.melted_fraction = X_stage(i);         % 更新当前步的熔化分数
                                tempState.T_p = obj.params.materials.Mg.melting_point; % 确保温度为熔点
                                state_stage_history(i) = tempState;
                            end
                        end
                        % --- 结束修改 ---
                        
                    case 'vaporization'
                        % 使用ode45进行自适应时间步长求解
                        vaporization_solver = VaporizationStage(obj.params, obj.physicalModel);
                        
                        % --- 开始修改 (1/3): 扩展初始状态向量y0以包含温度 ---
                        y0 = [obj.model.particleState.m_mg, obj.model.particleState.m_mgo, obj.model.particleState.T_p];
                        % --- 结束修改 ---

                        t_span = [current_time, obj.params.total_time];
                        
                        ode_options = odeset('RelTol', 1e-5, 'Events', @(t,y) obj.vaporization_events(t,y));
                        
                        [t_stage, y_stage] = ode45(@(t,y) obj.vaporization_ode(t, y, vaporization_solver), t_span, y0, ode_options);
                        
                        % 将质量和温度历史转换为状态历史
                        state_stage_history = obj.convert_mass_to_state_history(y_stage);
                end

                if isempty(t_stage)
                    fprintf('阶段 %s 未产生有效的时间步, 仿真可能已停滞。\n', current_stage_name);
                    break;
                end
                
                current_time = t_stage(end);
                fprintf('阶段 %s 完成, 结束时间: %.4fs\n', current_stage_name, current_time);
                
                % 3. 记录结果
                obj.recordResults(t_stage, state_stage_history, current_stage_name);

                % 4. 关键: 用当前阶段的最终状态更新主模型
                obj.model.particleState = state_stage_history(end);
                
                % 5. 检查是否需要结束
                if current_time >= obj.params.total_time
                    break;
                end
            end
            
            results = obj.prepareResults();
        end

        function dydt = vaporization_ode(obj, t, y, solver)
            % --- 开始修改 (2/3): 重写ODE函数以包含温度耦合 ---
            % y(1) -> m_mg, y(2) -> m_mgo, y(3) -> T_p
            
            % 1. 从当前状态y构建一个临时的ParticleState对象
            tempState = obj.model.particleState.copy(); % 复制以保留其他状态
            tempState.m_mg = y(1);
            tempState.m_mgo = y(2);
            tempState.T_p = y(3);
            tempState.update_geometry_from_mass(); % 根据新质量更新半径

            % 2. 调用BVP求解器获取瞬时反应速率
            try
                rate_info = solver.solve(tempState);
                dmdt_mg = -rate_info.m_dot_mg_reac; % 消耗率，所以为负
            catch ME
                fprintf('t=%.4f s 时BVP求解失败: %s\n', t, ME.message);
                dmdt_mg = 0; % 求解失败时，速率为0
            end

            % 3. 根据化学计量关系计算m_mgo的变化率
            mw_mg = obj.params.materials.Mg.molar_mass;
            mw_mgo = obj.params.materials.MgO.molar_mass;
            if dmdt_mg < 0
                dmol_mg_dt = -dmdt_mg / mw_mg; % Mg消耗的摩尔速率 (正值)
                dmdt_mgo = dmol_mg_dt * mw_mgo; % MgO生成的质量速率
            else
                dmdt_mgo = 0;
            end
            
            % 4. 计算温度变化率 dT/dt = (Q_reac - Q_loss) / C_p_total
            % 4.1 计算反应热 Q_reac = d(mol_mg)/dt * H_reac
            % 反应: Mg + 0.5 CO2 -> MgO + 0.5 C (假设)
            % 需要计算反应焓变 (J/mol_Mg)
            H_mg_g = obj.thermo_reader.calculate_H('Mg', tempState.T_p); % 气态Mg
            H_co2 = obj.thermo_reader.calculate_H('CO2', tempState.T_p);
            H_mgo_g = obj.thermo_reader.calculate_H('MgO', tempState.T_p); % 气态MgO
            H_c_s = 0; % 固态C的生成焓近似为0
            H_reac_per_mol_mg = H_mgo_g + 0.5 * H_c_s - H_mg_g - 0.5 * H_co2;
            
            if dmdt_mg < 0
                Q_reac = dmol_mg_dt * (-H_reac_per_mol_mg); % 放热为正
            else
                Q_reac = 0;
            end

            % 4.2 计算热损失 Q_loss
            q_conv = obj.params.h_conv * (tempState.T_p - obj.params.ambient_temperature);
            q_rad = obj.params.emissivity * obj.params.sigma * (tempState.T_p^4 - obj.params.ambient_temperature^4);
            A_p = 4 * pi * tempState.r_p^2;
            Q_loss = (q_conv + q_rad) * A_p;

            % 4.3 计算总热容 C_p_total
            Cp_total = obj.physicalModel.calculate_total_heat_capacity(tempState);
            
            % 4.4 计算dT/dt
            if Cp_total > 1e-9
                dTdt = (Q_reac - Q_loss) / Cp_total;
            else
                dTdt = 0;
            end

            dydt = [dmdt_mg; dmdt_mgo; dTdt];
            
            fprintf('  t=%.4f s, m_mg=%.2e, T_p=%.1f K, dm/dt=%.2e, dT/dt=%.2e\n', t, y(1), y(3), dmdt_mg, dTdt);
            % --- 结束修改 ---
        end
        
        function [value, isterminal, direction] = vaporization_events(obj, t, y)
            % ODE事件函数: 当Mg质量耗尽时停止
            value = y(1);      % 当y(1) (m_mg) 小于一个阈值时触发
            isterminal = 1;    % 触发后停止积分
            direction = -1;    % 从正向负穿过时触发
        end

        function recordResults(obj, t, state_history, stage_name)
            % 记录结果
            if isempty(t)
                return;
            end
            
            % 排除重复的第一个时间点
            if ~isempty(obj.results.time) && t(1) == obj.results.time(end)
                t = t(2:end);
                state_history = state_history(2:end);
            end

            if isempty(t)
                return;
            end
            
            obj.results.time = [obj.results.time; t(:)];
            obj.results.particleStateHistory = [obj.results.particleStateHistory, state_history];
            obj.results.stageHistory = [obj.results.stageHistory; repmat({stage_name}, length(t), 1)];
        end

        function state_history = convert_mass_to_state_history(obj, y_stage)
            % --- 开始修改 (3/3): 更新辅助函数以处理包含温度的状态向量 ---
            % y_stage 是一个 N x 3 的矩阵, [m_mg, m_mgo, T_p]
            num_steps = size(y_stage, 1);
            if num_steps == 0
                state_history = ParticleState.empty;
                return;
            end
            state_history(1, num_steps) = ParticleState(); % 使用默认构造函数进行预分配
            
            % 获取此阶段开始时的其他状态作为模板
            initial_state_template = obj.model.particleState.copy();
            
            for i = 1:num_steps
                tempState = initial_state_template.copy();
                tempState.m_mg = y_stage(i, 1);
                tempState.m_mgo = y_stage(i, 2);
                tempState.T_p = y_stage(i, 3);
                tempState.update_geometry_from_mass(); % 重要: 根据新质量更新几何
                
                state_history(i) = tempState;
            end
            % --- 结束修改 ---
        end



        %%%%%  本质即从运行计算的结果中读取想要的参数结果 转化为易于出图的结构形式。
        function results = prepareResults(obj)
            % 辅助函数: 将历史记录转换为易于绘图和分析的格式
            num_entries = length(obj.results.time);
            if num_entries == 0
                results = struct();
                return;
            end

            % 初始化数组
            results.time = obj.results.time;
            results.stage = obj.results.stageHistory;
            results.temperature = zeros(num_entries, 1);
            results.melted_fraction = zeros(num_entries, 1);
            results.mass_mg = zeros(num_entries, 1);
            results.mass_mgo = zeros(num_entries, 1);
            results.mass_c = zeros(num_entries, 1);
            results.radius = zeros(num_entries, 1);
            results.oxide_thickness = zeros(num_entries, 1); % 提取氧化层厚度

            % 填充数据
            for i = 1:num_entries
                state = obj.results.particleStateHistory(i);
                results.temperature(i) = state.T_p;
                results.melted_fraction(i) = state.melted_fraction;
                results.mass_mg(i) = state.m_mg;
                results.mass_mgo(i) = state.m_mgo;
                results.mass_c(i) = state.m_c;
                results.radius(i) = state.r_p;
                results.oxide_thickness(i) = state.oxide_thickness; % 提取氧化层厚度
            end
        end
    end
end 