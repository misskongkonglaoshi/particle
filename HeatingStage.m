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
            ode_func = @(t, T) obj.ode_system(t, T, particleState);

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
            options = odeset('Events', @(t, y) obj.target_temp_event(t, y, target_temp));
            
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
        
        function dTdt = ode_system(obj, t, T, particleState)
            % 更新粒子状态以便进行热容计算
            particleState.T_p = T;
            
            % 计算总热通量 (对流 + 辐射)
            Q_flux = obj.physicalModel.calculate_heat_flux(particleState);

            % 定义反应放热的接口 (目前为0)
            Q_reaction = 0;
            
            % 总热量
            Q_total = Q_flux + Q_reaction;

            % 计算总热容 (J/K)
            Cp_total = obj.physicalModel.calculate_total_heat_capacity(particleState);
            
            % 根据正确的物理公式 dT/dt = Q_total / Cp_total 计算
            if Cp_total < 1e-9
                dTdt = 0;
            else
                dTdt = Q_total / Cp_total;
            end
            
            % 确保返回的是列向量
            dTdt = dTdt(:);
        end
        
        function [value, isterminal, direction] = target_temp_event(~, t, y, target_temp)
            % 事件函数: 当温度达到目标温度时触发
            value = y(1) - target_temp;
            isterminal = 1; % 触发后终止求解
            direction = 1;  % 从负到正穿过时触发
        end
    end
end 