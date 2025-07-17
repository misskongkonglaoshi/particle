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
        vaporization_rate_info_cache % V12: 用于在气化阶段ODE求解中缓存BVP结果
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
                obj.results.flame_radius = [];      % 新增: 初始化火焰半径历史
                obj.results.flame_temperature = []; % 新增: 初始化火焰温度历史
            end
        end
        
        function results = solve(obj)
            % 求解整个过程
            
            % 清空历史记录
            obj.results.time = [];
            obj.results.particleStateHistory = [];
            obj.results.stageHistory = {};
            obj.results.flame_radius = [];      % 新增: 清空火焰半径历史
            obj.results.flame_temperature = []; % 新增: 清空火焰温度历史
            
            current_time = 0;
            
            % 循环求解各阶段
            while current_time < obj.params.total_time
                % 1. 确定当前阶段
                current_stage_name = obj.model.determineStage();
                
                % 2. 根据阶段动态创建求解器并求解
                fprintf('当前时间: %.4fs, 进入 %s 阶段...\n', current_time, current_stage_name);
                
                t_stage = [];
                state_stage_history = ParticleState.empty;
                flame_info_history = []; % 新增: 为火焰信息历史记录提供默认空值

                switch current_stage_name
                    case 'heating_and_melting'
                        % 使用统一的、基于焓的求解器
                        solver = UnifiedHeatingStage(obj.params, obj.physicalModel);
                            target_temp = obj.params.materials.Mg.ignition_temp;
                        fprintf('--- 进入统一加热/熔化阶段 ---\n');
                            fprintf('  (目标: 点火温度 %.1f K)\n', target_temp);
                        
                        % 前置条件检查
                        if obj.model.particleState.T_p >= target_temp
                            fprintf('  (跳过: 当前温度 %.1f K >= 目标温度 %.1f K)\n', obj.model.particleState.T_p, target_temp);
                            t_stage = [current_time, current_time];
                            state_stage_history = [obj.model.particleState, obj.model.particleState];
                        else
                            [t_stage, state_stage_history] = solver.solve(obj.model.particleState, current_time, target_temp);
                        end
                        
                    case 'vaporization'
                        % --- V12 修改: 增加OutputFcn来缓存BVP结果 ---

                        % 1. 将颗粒温度锁定在沸点
                        T_boil = obj.params.materials.Mg.ignition_temp;
                        obj.model.particleState.T_p = T_boil;
                        fprintf('--- 进入气相燃烧阶段, 颗粒温度锁定在 %.1f K ---\n', T_boil);

                        % 2. 初始化求解器和结果缓存
                        vaporization_solver = VaporizationStage(obj.params, obj.physicalModel);
                        obj.vaporization_rate_info_cache = {}; % V12: 初始化缓存
                        
                        % 3. 设置ODE的初始状态向量 [m_mg, m_mgo, m_c, r_c, r_p]
                        pState = obj.model.particleState;
                        y0 = [pState.m_mg, pState.m_mgo, pState.m_c, pState.r_c, pState.r_p];

                        % 4. 设置求解时间区间、事件和OutputFcn
                        t_span = [current_time, obj.params.total_time];
                        output_fcn_handle = @(t,y,flag) obj.vaporization_output_fcn(t, y, flag, vaporization_solver, T_boil);
                        ode_options = odeset('RelTol', 1e-5, 'Events', @(t,y) obj.vaporization_events(t,y), 'OutputFcn', output_fcn_handle);
                        
                        % 5. 求解状态向量变化的ODE
                        [t_stage, y_stage] = ode45(@(t,y) obj.vaporization_ode(t, y, vaporization_solver, T_boil), t_span, y0, ode_options);
                        
                        % 6. 将状态向量历史转换为状态对象历史, 并记录火焰信息
                        [state_stage_history, flame_info_history] = obj.convert_state_vector_to_state_history(t_stage, y_stage, T_boil, obj.vaporization_rate_info_cache);
                        
                        % 7. 清理缓存
                        obj.vaporization_rate_info_cache = {};

                        % --- 结束修改 ---
                end

                if isempty(t_stage)
                    fprintf('阶段 %s 未产生有效的时间步, 仿真可能已停滞。\n', current_stage_name);
                    break;
                end
                
                current_time = t_stage(end);
                fprintf('阶段 %s 完成, 结束时间: %.4fs\n', current_stage_name, current_time);
                
                % 3. 记录结果 (传入火焰信息)
                obj.recordResults(t_stage, state_stage_history, current_stage_name, flame_info_history);

                % 4. 关键: 用当前阶段的最终状态更新主模型
                obj.model.particleState = state_stage_history(end);
                
                % 5. 检查是否需要结束
                if current_time >= obj.params.total_time
                    break;
                end
            end
            
            results = obj.prepareResults();
        end

        function dydt = vaporization_ode(obj, t, y, solver, T_p_const)
            % --- V11 重构: 求解包含质量和半径的5变量状态向量 ---
            % y(1)=m_mg, y(2)=m_mgo, y(3)=m_c, y(4)=r_c, y(5)=r_p
            
            % 为防止数值问题，增加保护性判断
            if y(1) < 1e-15 || y(4) < 1e-9 || y(5) < 1e-9
                dydt = zeros(5,1);
                return;
            end

            % 1. 从当前状态向量y构建一个临时的ParticleState对象
            % 注意：这里不再调用update_geometry_from_mass
            tempState = obj.model.particleState.copy();
            tempState.m_mg = y(1);
            tempState.m_mgo = y(2);
            tempState.m_c = y(3);
            tempState.r_c = y(4);
            tempState.r_p = y(5);
            tempState.oxide_thickness = tempState.r_p - tempState.r_c;
            tempState.T_p = T_p_const;
            
            % 2. 调用BVP求解器获取质量变化率
            try
                rate_info = solver.solve_reaction_rates(tempState);
                dmdt_mg  = rate_info.dmdt_mg;
                dmdt_mgo = rate_info.dmdt_mgo;
                dmdt_c   = rate_info.dmdt_c;
            catch ME
                fprintf('t=%.4f s 时反应速率求解失败: %s\n', t, ME.message);
                % 重新抛出错误，这将中止ODE求解器
                rethrow(ME);
            end

            % 3. 新增核心逻辑: 计算半径变化率
            rho_mg = obj.params.materials.Mg.density;
            rho_mgo = obj.params.materials.MgO.density;
            rho_c = obj.params.materials.C.density;

            % 核心半径变化率 (基于Mg消耗)
            dVdt_mg = dmdt_mg / rho_mg;
            drc_dt = dVdt_mg / (4 * pi * y(4)^2 + eps);

            % 外部颗粒半径变化率 (基于产物沉积)
            dVdt_mgo = dmdt_mgo / rho_mgo;
            dVdt_c = dmdt_c / rho_c;
            dVdt_product_deposition = dVdt_mgo + dVdt_c;
            drp_dt = dVdt_product_deposition / (4 * pi * y(5)^2 + eps);

            % 4. 组合成完整的导数向量
            dydt = [dmdt_mg; dmdt_mgo; dmdt_c; drc_dt; drp_dt];
            
            fprintf('  t=%.4f s, m_mg=%.2e, r_p=%.2e, drp/dt=%.2e\n', t, y(1), y(5), drp_dt);
        end
        
        function [value, isterminal, direction] = vaporization_events(obj, t, y)
            % ODE事件函数: 当Mg质量耗尽时停止
            value = y(1);      % 当y(1) (m_mg) 小于一个阈值时触发
            isterminal = 1;    % 触发后停止积分
            direction = -1;    % 从正向负穿过时触发
        end

        function recordResults(obj, t, state_history, stage_name, flame_info_history)
            % --- 开始修改: 增加对火焰信息的记录 ---
            if nargin < 5
                flame_info_history = []; % 向后兼容
            end

            if isempty(t)
                return;
            end
            
            % 排除重复的第一个时间点
            if ~isempty(obj.results.time) && t(1) == obj.results.time(end)
                t = t(2:end);
                state_history = state_history(2:end);
                if ~isempty(flame_info_history)
                    flame_info_history = flame_info_history(2:end);
                end
            end

            if isempty(t)
                return;
            end
            
            obj.results.time = [obj.results.time; t(:)];
            obj.results.particleStateHistory = [obj.results.particleStateHistory, state_history];
            obj.results.stageHistory = [obj.results.stageHistory; repmat({stage_name}, length(t), 1)];
            
            % 记录火焰信息
            if ~isempty(flame_info_history)
                % V12: 确保flame_info_history是一个结构体数组
                if isstruct(flame_info_history)
                obj.results.flame_radius = [obj.results.flame_radius; [flame_info_history.r_f]'];
                obj.results.flame_temperature = [obj.results.flame_temperature; [flame_info_history.T_f]'];
                end
            else
                % 如果没有火焰信息，用NaN填充
                obj.results.flame_radius = [obj.results.flame_radius; nan(length(t), 1)];
                obj.results.flame_temperature = [obj.results.flame_temperature; nan(length(t), 1)];
            end
            % --- 结束修改 ---
        end

        function [state_history, flame_info_history] = convert_state_vector_to_state_history(obj, t_stage, y_stage, T_p_const, rate_info_cache)
            % --- V12 修改: 不再重复求解BVP，而是使用缓存的结果 ---
            % y_stage 是一个 N x 5 的矩阵 [m_mg, m_mgo, m_c, r_c, r_p]
            num_steps = size(y_stage, 1);
            if num_steps == 0
                state_history = ParticleState.empty;
                flame_info_history = [];
                return;
            end

            % V12: 检查缓存大小与步数是否匹配
            if num_steps ~= length(rate_info_cache)
                warning('StageManager:cacheMismatch', ...
                        'ode45返回的步数 (%d) 与OutputFcn缓存的速率信息数 (%d) 不匹配。火焰信息可能不准确。', ...
                        num_steps, length(rate_info_cache));
            end

            state_history(1, num_steps) = ParticleState(); % 预分配
            flame_info_history(1, num_steps) = struct('r_f', NaN, 'T_f', NaN);
            
            for i = 1:num_steps
                % 1. 从状态向量填充ParticleState对象
                currentState = y_stage(i, :);
                tempState = obj.build_temp_state_from_vector(currentState, T_p_const);
                state_history(i) = tempState;
                
                % 2. V12: 直接从缓存中读取火焰信息
                if i <= length(rate_info_cache)
                    rate_info = rate_info_cache{i};
                    flame_info_history(i).r_f = rate_info.r_f;
                    flame_info_history(i).T_f = rate_info.T_f;
                end
            end
        end


        % --- V12 新增: ODE求解相关的辅助函数 ---
        function status = vaporization_output_fcn(obj, t, y, flag, solver, T_p_const)
            % 作为ode45的OutputFcn, 在每个成功的时间步计算并缓存BVP结果
            status = 0; % 默认继续积分
            try
                switch flag
                    case 'init'
                        % 在积分开始时调用
                        % y是初始状态向量y0
                        tempState = obj.build_temp_state_from_vector(y, T_p_const);
                        rate_info = solver.solve_reaction_rates(tempState);
                        obj.vaporization_rate_info_cache{1} = rate_info;
                    case ''
                        % 在每个积分步后调用
                        % y可能包含多个时间步的结果，y的每一列是一个时间点的状态
                        for i = 1:size(y, 2)
                           tempState = obj.build_temp_state_from_vector(y(:,i), T_p_const);
                           rate_info = solver.solve_reaction_rates(tempState);
                           obj.vaporization_rate_info_cache{end+1} = rate_info;
                        end
                end
            catch ME
                fprintf('OutputFcn在t=%.4f s计算BVP时失败: %s. 终止积分。\n', t(1), ME.message);
                status = 1; % 返回1以终止积分
            end
        end

        function tempState = build_temp_state_from_vector(obj, y_vec, T_p_const)
             % 辅助函数: 从状态向量构建一个临时的ParticleState对象
             tempState = ParticleState(obj.params);
             tempState.m_mg  = y_vec(1);
             tempState.m_mgo = y_vec(2);
             tempState.m_c   = y_vec(3);
             tempState.r_c   = y_vec(4);
             tempState.r_p   = y_vec(5);
             tempState.oxide_thickness = tempState.r_p - tempState.r_c;
             tempState.T_p = T_p_const;
             tempState.melted_fraction = 1.0;
        end
        % --- 结束 V12 新增 ---


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
            results.stage = obj.results.stageHistory; % 保留原始的计算阶段
            results.physical_stage = cell(num_entries, 1); % 新增物理阶段描述
            results.temperature = zeros(num_entries, 1);
            results.melted_fraction = zeros(num_entries, 1);
            results.mass_mg = zeros(num_entries, 1);
            results.mass_mgo = zeros(num_entries, 1);
            results.mass_c = zeros(num_entries, 1); % 新增: 碳质量
            results.radius = zeros(num_entries, 1);
            results.oxide_thickness = zeros(num_entries, 1);
            results.flame_radius = []; % 新增: 火焰半径
            results.flame_temperature = []; % 新增: 火焰温度

            % 填充数据
            for i = 1:num_entries
                state = obj.results.particleStateHistory(i);
                results.temperature(i) = state.T_p;
                results.melted_fraction(i) = state.melted_fraction;
                results.mass_mg(i) = state.m_mg;
                results.mass_mgo(i) = state.m_mgo;
                results.mass_c(i) = state.m_c; % 新增
                results.radius(i) = state.r_p;
                results.oxide_thickness(i) = state.oxide_thickness;
                
                % 根据温度和熔化分数追认物理阶段
                if strcmp(results.stage{i}, 'heating_and_melting')
                    T = state.T_p;
                    X = state.melted_fraction;
                    T_melt = obj.params.materials.Mg.melting_point;

                    if T < T_melt - 1e-6
                        results.physical_stage{i} = 'Preheating';
                    elseif abs(T - T_melt) < 1e-6 && X < 1.0
                        results.physical_stage{i} = 'Melting';
                    else
                        results.physical_stage{i} = 'Liquid Heating';
                    end
                else
                    % 对于其他阶段(如vaporization), 直接使用计算阶段名
                    results.physical_stage{i} = results.stage{i};
                end
            end
            
            % 直接从记录中附加火焰数据
            if isfield(obj.results, 'flame_radius')
                results.flame_radius = obj.results.flame_radius;
                results.flame_temperature = obj.results.flame_temperature;
            end
        end
    end
end 