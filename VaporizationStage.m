classdef VaporizationStage < handle
    % VaporizationStage (气化阶段求解器) - V5 (最终版: 统一BVP求解)
    %
    % 核心功能:
    % 对于一个处于沸点的、给定状态的颗粒，通过一次BVP调用，求解一个在物理上
    % 完全自洽的稳态燃烧模型，同时得到反应速率和火焰结构。
    %
    % 求解逻辑:
    % 使用bvp4c求解一个包含未知参数的边界值问题。
    % 未知数 = [10个状态变量(T, Yi, d/dr), 火焰半径(r_f), 总蒸发速率(m_dot)]
    % 方程数 = [10个ODE, 12个边界/界面条件]
    % 通过这种方式，所有的物理定律（能量守恒、质量守恒、化学计量）都在
    % 一个统一的框架内被同时满足。

    properties
        params          % 参数对象
        physicalModel   % 物理模型对象
        thermo          % 热力学数据读取器
    end

    methods
        function obj = VaporizationStage(params, physicalModel)
            % 构造函数
            obj.params = params;
            obj.physicalModel = physicalModel;
            obj.thermo = physicalModel.thermo_reader;
        end

        function rate_info = solve_reaction_rates(obj, pState)
            % --- 主求解函数 ---
            % 将所有依赖项打包到一个结构体中
            bvp_deps.pState = pState;
            bvp_deps.params = obj.params;
            bvp_deps.physicalModel = obj.physicalModel;
            bvp_deps.thermo = obj.thermo;
            
            try
                % 调用完全解耦的外部BVP求解器
                rate_info = solve_vaporization_bvp(bvp_deps);
            catch ME
                fprintf('--- BVP 求解器失败 ---\n');
                fprintf('错误信息: %s\n', ME.message);
                
                % 在BVP失败时输出关键尺度参数
                r_p = pState.r_p;
                r_inf = 20 * r_p; % 与solve_vaporization_bvp.m中保持一致
                fprintf('  > 关键参数: r_p = %.3e m, r_inf = %.3e m (尺度差异: %.1f 倍)\n', ...
                        r_p, r_inf, r_inf/r_p);

                fprintf('详细错误堆栈:\n');
                for i = 1:length(ME.stack)
                    fprintf('  > 文件: %s, 行: %d, 函数: %s\n', ...
                        ME.stack(i).file, ME.stack(i).line, ME.stack(i).name);
        end
                
                % 重新抛出错误以终止整个计算
                rethrow(ME);
            end
        end
    end
end