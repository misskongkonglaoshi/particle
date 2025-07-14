classdef HeatingStage < handle
    % HeatingStage 熔融后加热阶段求解器
    
    properties
        params         % 参数对象
        physicalModel  % 物理模型计算器
    end
    
    methods
        function obj = HeatingStage(params, physicalModel)
            % 构造函数
            if nargin > 0
                obj.params = params;
                obj.physicalModel = physicalModel;
            end
        end
        
        function [t_history, T_history, t_end] = solve(obj, particleState, t_start, target_temp)
            % 求解加热阶段，直到达到目标温度
            %
            % 输入:
            %   particleState: 当前的颗粒状态对象
            %   t_start: 阶段开始时间
            %   target_temp: 该阶段要达到的目标温度
            % 输出:
            %   t_history: 该阶段的时间历史
            %   T_history: 该阶段的温度历史
            %   t_end: 阶段结束时间点

            % 确保熔化分数是1
            particleState.melted_fraction = 1.0;

            % 定义求解的ODE
            ode_func = @(t, T) obj.ode_system(T, particleState);

            % 设置求解时间跨度和初始条件
            tspan = [t_start, obj.params.total_time];
            y0 = particleState.T_p;

            % --- 诊断代码 ---
            if isempty(y0) || ~isnumeric(y0) || ~isscalar(y0) || isnan(y0)
                error('HeatingStage:InvalidInitialTemperature', ...
                      '传递给ode45的初始温度y0无效。值为: %s。请检查Parameters.m或上游计算。', ...
                      mat2str(y0));
            end
            % --- 结束诊断 ---

            % 设置事件函数，当温度达到动态指定的目标温度时停止
            options = odeset('Events', @(t, y) obj.reach_target_temp_event(t, y, target_temp));
            
            % 调用ODE求解器
            [t_history, T_history, te, ~, ~] = ode45(ode_func, tspan, y0, options);
            
            % 更新状态并返回结束时间
            if ~isempty(te)
                t_end = te(1);
                particleState.T_p = T_history(end);
            else
                t_end = t_history(end);
                particleState.T_p = T_history(end);
            end
        end
        
        function dTdt = ode_system(obj, T, particleState)
            % 加热阶段的ODE系统
            % fprintf('开始加热阶段ode...\n');
            % 更新瞬时温度以进行计算
            particleState.T_p = T;
            
            % 计算热通量 (与预热阶段相同)
            q_conv = obj.params.h_conv * (obj.params.ambient_temperature - T);
            q_rad = obj.params.emissivity * obj.params.sigma * (obj.params.ambient_temperature^4 - T^4);
            
            A_p = 4 * pi * particleState.r_p^2;
            Q_total = (q_conv + q_rad) * A_p;
            
            % 计算热容
            Cp = obj.physicalModel.calculate_total_heat_capacity(particleState);
            
            if Cp > 1e-9
                dTdt = Q_total / Cp;
            else
                dTdt = 0;
            end
        end
        
        function [value, isterminal, direction] = reach_target_temp_event(obj, ~, y, target_temp)
            % 事件函数: 当温度达到目标温度时触发
            value = y(1) - target_temp;
            isterminal = 1; % 触发后终止求解
            direction = 1;  % 从负到正穿过时触发
        end
    end
end 