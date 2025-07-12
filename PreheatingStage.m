classdef PreheatingStage < handle
    % PreheatingStage 预热阶段求解器
    
    properties
        params         % 参数对象
        physicalModel  % 物理模型计算器
    end
    
    methods
        function obj = PreheatingStage(params, physicalModel)
            % 构造函数
            obj.params = params;
            obj.physicalModel = physicalModel;
        end
        
        function [t_history, T_history, t_end] = solve(obj, particleState, t_start)
            % 求解预热阶段
            %
            % 输入:
            %   particleState: 当前的颗粒状态对象
            %   t_start: 阶段开始时间
            % 输出:
            %   t_history: 该阶段的时间历史
            %   T_history: 该阶段的温度历史
            %   t_end: 阶段结束时间点

            % 定义求解的ODE
            ode_func = @(t, T) obj.ode_system(T, particleState);

            % 设置求解时间跨度和初始条件
            tspan = [t_start, obj.params.total_time];
            y0 = particleState.T_p;

            % 设置事件函数，当温度达到熔点时停止
            options = odeset('Events', @(t, y) obj.reach_melting_point_event(t, y));
            
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
            % 预热阶段的ODE系统
            
            % 更新瞬时温度以进行计算
            particleState.T_p = T;
            
            % 计算热通量
            q_conv = obj.params.h_conv * (obj.params.ambient_temperature - T);
            q_rad = obj.params.emissivity * obj.params.sigma * (obj.params.ambient_temperature^4 - T^4);
            
            % 计算总热量和温度变化率
            A_p = 4 * pi * particleState.r_p^2;
            Q_total = (q_conv + q_rad) * A_p;
            
            Cp = obj.physicalModel.calculate_total_heat_capacity(particleState);
            
            if Cp > 1e-9
                dTdt = Q_total / Cp;
            else
                dTdt = 0;
            end
        end
        
        function [value, isterminal, direction] = reach_melting_point_event(obj, ~, y)
            % 事件函数: 当温度达到熔点时触发
            value = y(1) - obj.params.materials.Mg.melting_point; % 当value=0时事件触发
            isterminal = 1; % 触发后终止求解
            direction = 1;  % 从负到正穿过时触发
        end
    end
end 