classdef UnifiedHeatingStage < handle
    % UnifiedHeatingStage 统一的加热/熔化阶段求解器 (基于焓方法)
    
    properties
        params
        physicalModel
        initialStateTemplate % 用于在事件函数中作为模板
    end
    
    methods
        function obj = UnifiedHeatingStage(params, physicalModel)
            obj.params = params;
            obj.physicalModel = physicalModel;
        end
        
        function [t_history, state_history] = solve(obj, particleState, t_start, target_temp)
            % 求解统一的加热/熔化阶段, 直到达到目标温度
            
            % 将初始状态保存为模板，供事件函数使用
            obj.initialStateTemplate = particleState;
            
            % 1. 根据初始状态计算初始总焓
            H0 = obj.physicalModel.get_enthalpy_from_state(particleState);
            
            % 2. 定义求解总焓H的ODE
            ode_func = @(t, H) obj.ode_system(H, particleState);

            % 3. 定义一个用于平滑输出的时间向量，而不是简单的起止点
            output_dt = 1e-4; % 输出步长 (s), 可根据需要调整以获得更平滑的曲线
            tspan = t_start:output_dt:obj.params.total_time;

            % 4. 设置事件函数: 当温度达到目标时停止 (通过焓来判断)
            options = odeset('Events', @(t, H) obj.target_temp_event(t, H, target_temp));
            
            % 5. 调用ODE求解器求解 H(t)
            [t_history, H_history] = ode45(ode_func, tspan, H0, options);
            
            % 6. 将焓的历史记录转换为完整的状态历史记录
            num_steps = length(t_history);
            if num_steps > 0
                state_history(1, num_steps) = ParticleState(); % 预分配
                for i = 1:num_steps
                    state_history(i) = obj.physicalModel.get_state_from_enthalpy(H_history(i), particleState);
                end
            else
                state_history = ParticleState.empty;
            end
        end
        
        function dHdt = ode_system(obj, H, particleState)
            % 统一加热阶段的ODE系统: dH/dt = Q_total
            
            % a) 从当前的总焓H反算出瞬时状态, 主要是为了得到温度
            tempState = obj.physicalModel.get_state_from_enthalpy(H, particleState);
            
            % b) 根据这个瞬时状态(主要是温度)计算净热流
            dHdt = obj.physicalModel.calculate_heat_flux(tempState);
        end
        
        function [value, isterminal, direction] = target_temp_event(obj, ~, H, target_temp)
            % 事件函数: 当根据焓算出的温度达到目标温度时触发
            
            % a) 从当前的总焓H反算出瞬时状态 (使用保存的模板)
            tempState = obj.physicalModel.get_state_from_enthalpy(H, obj.initialStateTemplate);
            
            % b) 比较当前温度和目标温度
            value = tempState.T_p - target_temp;
            isterminal = 1;
            direction = 1;
        end
    end
end 