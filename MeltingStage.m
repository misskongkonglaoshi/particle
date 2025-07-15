classdef MeltingStage < handle
    % MeltingStage 熔融阶段求解器
    
    properties
        params         % 参数对象
        physicalModel  % 物理模型计算器
    end
    
    methods
        function obj = MeltingStage(params, physicalModel)
            % 构造函数
            if nargin > 0
                obj.params = params;
                obj.physicalModel = physicalModel;
            end
        end
        
        function [t_history, X_history, t_end] = solve(obj, particleState, t_start)
            % 求解熔融阶段 (熔化分数随时间的变化)
            %
            % 输入:
            %   particleState: 当前的颗粒状态对象
            %   t_start: 阶段开始时间
            % 输出:
            %   t_history: 该阶段的时间历史
            %   X_history: 该阶段的熔化分数历史
            %   t_end: 阶段结束时间点

            % 在此阶段, 温度保持在熔点
            particleState.T_p = obj.params.materials.Mg.melting_point;

            % 定义求解的ODE (求解熔化分数X的变化率)
            ode_func = @(t, X) obj.ode_system(X, particleState);

            % 设置求解时间跨度和初始条件
            tspan = [t_start, obj.params.total_time];
            y0 = particleState.melted_fraction;

            % 设置事件函数，当熔化分数达到1时停止
            options = odeset('Events', @(t, y) obj.fully_melted_event(t, y));
            
            % 调用ODE求解器
            [t_history, X_history, te, ~, ~] = ode45(ode_func, tspan, y0, options);
            
            % 更新状态并返回结束时间
            if ~isempty(te)
                t_end = te(1);
                particleState.melted_fraction = X_history(end);
            else
                t_end = t_history(end);
                particleState.melted_fraction = X_history(end);
            end
            % 确保温度是熔点
            particleState.T_p = obj.params.materials.Mg.melting_point;
        end
        
        function dXdt = ode_system(obj, X, particleState)
            % 熔融阶段的ODE系统
            fprintf('开始求解熔融阶段ode...\n');
            % 计算用于熔化的热通量
            q_conv = obj.params.h_conv * (obj.params.ambient_temperature - particleState.T_p);
            q_rad = obj.params.emissivity * obj.params.sigma * (obj.params.ambient_temperature^4 - particleState.T_p^4);
            
            A_p = 4 * pi * particleState.r_p^2;
            Q_in = (q_conv + q_rad) * A_p;
            
            % 计算熔化所需总能量
            L_m = obj.params.materials.Mg.latent_heat; % 熔化潜热
            m_total_mg = particleState.m_mg; % 只有Mg熔化
            
            % 熔化分数的变化率 dX/dt = Q_in / (m_mg * L_m)
            if m_total_mg > 1e-12 && L_m > 0
                dXdt = Q_in / (m_total_mg * L_m);
            else
                dXdt = 0;
            end
        end
        
        function [value, isterminal, direction] = fully_melted_event(obj, ~, y)
            % 事件函数: 当熔化分数达到1时触发
            value = y(1) - 1.0; % 当value=0时事件触发
            isterminal = 1;     % 触发后终止求解
            direction = 1;      % 从负到正穿过时触发
        end
    end
end 