classdef PreheatingStage < handle
    % PreheatingStage 预热阶段求解器 (重构后)
    
    properties
        params         % 参数对象
        physicalModel  % 物理模型计算器
    end
    
    methods
        function obj = PreheatingStage(params, physicalModel)
            % 构造函数
            if nargin > 0
                obj.params = params;
                obj.physicalModel = physicalModel;
            end
        end
        
        function [t_history, T_history, t_end] = solve(obj, particleState, t_start, target_temp)
            % 求解预热阶段, 直到达到目标温度
            
            % 定义求解的ODE
            ode_func = @(t, T) obj.ode_system(t, T, particleState);

            % 设置求解时间跨度和初始条件
            tspan = [t_start, obj.params.total_time];
            y0 = particleState.T_p;

            % 设置事件函数，当温度达到动态指定的目标温度时停止
            options = odeset('Events', @(t, y) obj.target_temp_event(t, y, target_temp));
            
            % 调用ODE求解器
            [t_history, T_history, te, ~, ~] = ode45(ode_func, tspan, y0, options);
            
            % 更新状态并返回结束时间
            if ~isempty(te)
                t_end = te(1);
            else
                t_end = t_history(end);
            end
            % 确保最终状态被设置
            if ~isempty(T_history)
                particleState.T_p = T_history(end);
            end
        end
        
        function dTdt = ode_system(obj, t, T, particleState)
            % 预热阶段的ODE系统
            
            % 更新瞬时温度以进行计算
            particleState.T_p = T;
            
            % 使用物理模型计算总热通量
            Q_total = obj.physicalModel.calculate_heat_flux(particleState);
            
            % 计算总热容 (J/K)
            Cp_total = obj.physicalModel.calculate_total_heat_capacity(particleState);

            % 根据物理公式 dT/dt = Q_total / Cp_total 计算温度变化率
            if Cp_total < 1e-9 % 避免除以零
                dTdt = 0;
            else
                dTdt = Q_total / Cp_total;
            end

            % 确保返回的是列向量
            dTdt = dTdt(:);
        end
        
        function [value, isterminal, direction] = target_temp_event(obj, t, y, target_temp)
            % 事件函数: 当温度达到目标温度时触发
            value = y(1) - target_temp; % 当value=0时事件触发
            isterminal = 1; % 触发后终止求解
            direction = 1;  % 从负到正穿过时触发
        end
    end
end 