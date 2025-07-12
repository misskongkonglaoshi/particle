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
            obj.params = params;
            
            % 初始化依赖的模型
            obj.thermo_reader = ThermoReader('thermo.dat');
            obj.physicalModel = PhysicalModel(obj.params, obj.thermo_reader);
            obj.model = ParticleModel(obj.params);
            
            % 初始化结果存储
            obj.results.time = [];
            obj.results.particleStateHistory = [];
            obj.results.stageHistory = {};
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
                stage_name = obj.model.determineStage();
                
                % 2. 根据阶段动态创建求解器并求解
                fprintf('当前时间: %.4fs, 进入 %s 阶段...\n', current_time, stage_name);
                
                t_stage = [];
                state_stage_history = [];
                
                switch stage_name
                    case 'preheating'
                        solver = PreheatingStage(obj.params, obj.physicalModel);
                        [t_stage, T_stage, current_time] = solver.solve(obj.model.particleState, current_time);
                        % 将温度历史转换为状态历史
                        state_stage_history = obj.convert_T_to_state_history(t_stage, T_stage);

                    case 'melting'
                        solver = MeltingStage(obj.params, obj.physicalModel);
                        [t_stage, X_stage, current_time] = solver.solve(obj.model.particleState, current_time);
                        % 将熔化分数历史转换为状态历史
                        state_stage_history = obj.convert_X_to_state_history(t_stage, X_stage);
                        
                    case 'heating'
                        solver = HeatingStage(obj.params, obj.physicalModel);
                        [t_stage, T_stage, current_time] = solver.solve(obj.model.particleState, current_time);
                        % 将温度历史转换为状态历史
                        state_stage_history = obj.convert_T_to_state_history(t_stage, T_stage);
                        
                    case 'vaporization'
                        solver = VaporizationStage(obj.params, obj.physicalModel);
                        [t_stage, state_stage_history] = solver.solve(obj.model.particleState, current_time);
                        current_time = obj.params.total_time; % 气相燃烧后直接结束
                end
                
                fprintf('阶段 %s 完成, 结束时间: %.4fs\n', stage_name, current_time);
                
                % 3. 记录结果
                obj.recordResults(t_stage, state_stage_history, stage_name);
                
                % 如果已经到了总时间或气相燃烧阶段，结束循环
                if strcmp(stage_name, 'vaporization') || current_time >= obj.params.total_time
                    break;
                end
            end
            
            % 准备返回结果
            results = obj.results;
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

        function state_history = convert_T_to_state_history(obj, t_stage, T_stage)
            % 辅助函数：将温度历史转换为状态历史
            num_steps = length(t_stage);
            state_history(1, num_steps) = ParticleState(obj.params);
            for i = 1:num_steps
                tempState = obj.model.particleState.copy();
                tempState.T_p = T_stage(i);
                state_history(i) = tempState;
            end
        end

        function state_history = convert_X_to_state_history(obj, t_stage, X_stage)
            % 辅助函数：将熔化分数历史转换为状态历史
            num_steps = length(t_stage);
            state_history(1, num_steps) = ParticleState(obj.params);
            for i = 1:num_steps
                tempState = obj.model.particleState.copy();
                tempState.melted_fraction = X_stage(i);
                state_history(i) = tempState;
            end
        end
    end
end 