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
                        solver = HeatingStage(obj.params, obj.physicalModel);
                        if strcmp(current_stage_name, 'preheating')
                            target_temp = obj.params.materials.Mg.melting_point;
                            fprintf('  (目标: 熔点 %.1f K)\n', target_temp);
                        else % liquid_heating
                            target_temp = obj.params.materials.Mg.ignition_temp;
                            fprintf('  (目标: 点火温度 %.1f K)\n', target_temp);
                        end
                        [t_stage, T_stage] = solver.solve(obj.model.particleState, current_time, target_temp);
                        
                        % --- 开始修改 ---
                        % 直接根据HeatingStage的输出构建正确的状态历史
                        num_steps = length(T_stage);
                        state_stage_history(1, num_steps) = ParticleState(); % 预分配
                        for i = 1:num_steps
                            tempState = obj.model.particleState.copy(); % 复制通用状态
                            tempState.T_p = T_stage(i);                 % 更新当前步的温度
                            state_stage_history(i) = tempState;
                        end
                        % --- 结束修改 ---

                    case 'melting'
                        solver = MeltingStage(obj.params, obj.physicalModel);
                        [t_stage, X_stage] = solver.solve(obj.model.particleState, current_time);
                        
                        % --- 开始修改 ---
                        % 直接根据MeltingStage的输出构建正确的状态历史
                        num_steps = length(X_stage);
                        state_stage_history(1, num_steps) = ParticleState(); % 预分配
                        for i = 1:num_steps
                            tempState = obj.model.particleState.copy(); % 复制通用状态
                            tempState.melted_fraction = X_stage(i);     % 更新当前步的熔化分数
                            state_stage_history(i) = tempState;
                        end
                        % --- 结束修改 ---
                        
                    case 'vaporization'
                        % 使用ode45进行自适应时间步长求解
                        vaporization_solver = VaporizationStage(obj.params, obj.physicalModel);
                        % 状态向量 y = [m_mg, m_mgo]
                        y0 = [obj.model.particleState.m_mg, obj.model.particleState.m_mgo];
                        t_span = [current_time, obj.params.total_time];
                        
                        ode_options = odeset('RelTol', 1e-5, 'Events', @(t,y) obj.vaporization_events(t,y));
                        
                        [t_stage, y_stage] = ode45(@(t,y) obj.vaporization_ode(t, y, vaporization_solver), t_span, y0, ode_options);
                        
                        % 将质量历史转换为状态历史
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
            
            % 准备返回结果
            results = obj.prepareResults();
        end

        function dydt = vaporization_ode(obj, t, y, solver)
            % ODE函数句柄，用于ode45
            % y(1) -> m_mg (颗粒核心Mg的质量)
            % y(2) -> m_mgo (颗粒上MgO的质量)
            
            % 1. 从当前状态y构建一个临时的ParticleState对象
            tempState = obj.model.particleState.copy(); % 复制以保留其他状态
            tempState.m_mg = y(1);
            tempState.m_mgo = y(2);
            tempState.update_geometry_from_mass(); % 根据新质量更新半径

            % 2. 调用BVP求解器获取瞬时反应速率
            try
                rate_info = solver.solve(tempState);
                dmdt_mg = -rate_info.m_dot_mg_reac; % 消耗率，所以为负
                % fprintf('反应mg消耗量：%.4f \n',dmdt_mg);
            catch ME
                fprintf('t=%.4f s 时BVP求解失败: %s\n', t, ME.message);
                %fprintf('m_dot_mg_reac=%.4f s 时BVP求解失败: %s\n',m_dot_mg_reac, ME.message);
                dmdt_mg = 0; % 求解失败时，速率为0，让ode45减小步长重试
            end

            % 3. 根据化学计量关系计算m_mgo的变化率
            % 反应: 2Mg + CO2 -> 2MgO + C, 摩尔质量比为 1:1
            mw_mg = obj.params.materials.Mg.molar_mass;
            mw_mgo = obj.params.materials.MgO.molar_mass;
            dmol_mg_dt = -dmdt_mg / mw_mg; % Mg消耗的摩尔速率 (正值)
            dmdt_mgo = dmol_mg_dt * mw_mgo; % MgO生成的质量速率

            dydt = [dmdt_mg; dmdt_mgo];
            
            fprintf('  t=%.4f s, m_mg=%.2e kg, m_mgo=%.2e kg, dm_mg/dt=%.2e kg/s\n', t, y(1), y(2), dmdt_mg);
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
            % 辅助函数: 将质量历史转换为状态历史
            % y_stage 是一个 N x 2 的矩阵, [m_mg, m_mgo]
            num_steps = size(y_stage, 1);
            state_history(1, num_steps) = ParticleState(); % 使用默认构造函数进行预分配
            
            % 获取此阶段开始时的其他状态作为模板
            initial_state_template = obj.model.particleState.copy();
            
            for i = 1:num_steps
                tempState = initial_state_template.copy();
                tempState.m_mg = y_stage(i, 1);
                tempState.m_mgo = y_stage(i, 2);
                tempState.update_geometry_from_mass(); % 重要: 根据新质量更新几何
                
                state_history(i) = tempState;
            end
        end

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

            % 填充数据
            for i = 1:num_entries
                state = obj.results.particleStateHistory(i);
                results.temperature(i) = state.T_p;
                results.melted_fraction(i) = state.melted_fraction;
                results.mass_mg(i) = state.m_mg;
                results.mass_mgo(i) = state.m_mgo;
                results.mass_c(i) = state.m_c;
                results.radius(i) = state.r_p;
            end
        end
    end
end 