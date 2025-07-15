classdef PhysicalModel < handle
    % PhysicalModel 物理模型计算中心
    % 集中管理所有与物理、化学、热力学相关的计算。
    
    properties
        params          % 参数对象
        thermo_reader   % 热力学数据读取器
    end
    
    methods
        function obj = PhysicalModel(params, thermo_reader)
            % 构造函数
            %
            % 输入:
            %   params: 参数对象
            %   thermo_reader: ThermoReader 实例
            if nargin > 0
            obj.params = params;
            obj.thermo_reader = thermo_reader;
            end
        end
        
        function cp_total = calculate_total_heat_capacity(obj, particleState)
            % 计算颗粒的总热容 (J/K)
            %
            % 输入:
            %   particleState: 当前的颗粒状态对象
            % 输出:
            %   cp_total: 总热容 (J/K)
            
            % 计算各组分的摩尔数
            n_mg = particleState.m_mg / obj.params.materials.Mg.molar_mass;
            n_mgo = particleState.m_mgo / obj.params.materials.MgO.molar_mass;
            n_total = n_mg + n_mgo;
            
            if n_total > 1e-20
                % 计算摩尔分数
                x_mg = n_mg / n_total;
                x_mgo = n_mgo / n_total;
                
                % 计算当前温度下的摩尔比热容
                cp_mg_mol = obj.calculate_specific_heat('Mg', particleState.T_p);
                cp_mgo_mol = obj.calculate_specific_heat('MgO', particleState.T_p);
                
                % 按摩尔分数加权计算总的摩尔比热容
                cp_mol = x_mg * cp_mg_mol + x_mgo * cp_mgo_mol;
                
                % 转换为总热容
                cp_total = cp_mol * n_total;
            else
                cp_total = 0;
            end
            
            % 添加碳的贡献 (使用常数比热)
            cp_total = cp_total + particleState.m_c * obj.params.materials.C.heat_capacity;
        end
        
        function cp = calculate_specific_heat(obj, species, T)
            % 使用thermo_reader计算特定温度下的摩尔比热容
            %
            % 输入:
            %   species: 物质名称 ('Mg' 或 'MgO')
            %   T: 温度 (K)
            % 输出:
            %   cp: 摩尔比热容 (J/(mol·K))
            
            % 如果没有提供thermo_reader或缺少物种数据，则使用常数比热
            if isempty(obj.thermo_reader) || ~obj.thermo_reader.has_species(species)
                cp = obj.params.materials.(species).heat_capacity * obj.params.materials.(species).molar_mass;
                return;
            end
            
            % 直接调用ThermoReader中的方法进行计算
            cp = obj.thermo_reader.calculate_Cp(species, T);
            
            % 防御性修复: 确保返回的是一个标量
            % 即使 calculate_Cp 错误地返回一个向量, 这里也会将其求和为标量
            if ~isscalar(cp)
                cp = sum(cp);
            end
        end
        
        function Q_total = calculate_heat_flux(obj, particleState)
            % 计算颗粒与环境之间的总热通量 (W)
            %
            % 输入:
            %   particleState: 当前的颗粒状态对象
            % 输出:
            %   Q_total: 总热通量 (W), >0表示吸热

            T_p = particleState.T_p;
            A_p = 4 * pi * particleState.r_p^2;

            % 对流换热
            q_conv = obj.params.h_conv * (obj.params.ambient_temperature - T_p);
            
            % 辐射换热
            q_rad = obj.params.emissivity * obj.params.sigma * (obj.params.ambient_temperature^4 - T_p^4);
            
            % 总热通量 (W)
            Q_total = (q_conv + q_rad) * A_p;
        end
        
        function props = get_gas_properties(obj, T_gas)
            % 获取在特定温度下的气体属性
            % TODO: 实现更复杂的温度依赖模型
            % 目前, 返回在Parameters中定义的常数值
            props = obj.params.gas_properties;
        end

        function H = get_enthalpy_from_state(obj, particleState)
            % 根据颗粒状态（温度和熔化分数）计算其总焓 (J)
            % 参考点: 固态物质在0K时焓为0
            
            T = particleState.T_p;
            X = particleState.melted_fraction;
            T_melt = obj.params.materials.Mg.melting_point;
            
            % 我们需要一个在0K到T_melt之间的平均总热容
            % 为简单起见, 我们使用在初始温度下的热容作为整个固相的常数热容
            tempState_solid = particleState.copy();
            tempState_solid.T_p = obj.params.initial_temperature;
            Cp_total_solid = obj.calculate_total_heat_capacity(tempState_solid);

            % 熔点时固相的焓
            H_solid_at_melt = Cp_total_solid * T_melt;
            
            % 熔化潜热
            L_m_total = particleState.m_mg * obj.params.materials.Mg.latent_heat;

            if X < 1e-9 % 固相
                H = Cp_total_solid * T;
            elseif X >= 1e-9 && X < 1 % 熔化中
                H = H_solid_at_melt + X * L_m_total;
            else % 液相
                tempState_liquid = particleState.copy();
                tempState_liquid.T_p = T_melt; % 使用熔点温度估算液相热容
                Cp_total_liquid = obj.calculate_total_heat_capacity(tempState_liquid);
                
                H_liquid_at_melt = H_solid_at_melt + L_m_total;
                H = H_liquid_at_melt + Cp_total_liquid * (T - T_melt);
            end
        end

        function pStateOut = get_state_from_enthalpy(obj, H, pStateIn)
            % 根据总焓 H 反算颗粒的状态（温度和熔化分数）
            pStateOut = pStateIn.copy();

            % 使用与get_enthalpy_from_state一致的逻辑和近似
            tempState_solid = pStateIn.copy();
            tempState_solid.T_p = obj.params.initial_temperature;
            Cp_total_solid = obj.calculate_total_heat_capacity(tempState_solid);
            
            T_melt = obj.params.materials.Mg.melting_point;
            
            H_solid_at_melt = Cp_total_solid * T_melt;
            L_m_total = pStateIn.m_mg * obj.params.materials.Mg.latent_heat;
            
            if L_m_total < 1e-12 % 防止除零
                H_liquid_at_melt = H_solid_at_melt;
            else
                H_liquid_at_melt = H_solid_at_melt + L_m_total;
            end

            if H < H_solid_at_melt
                pStateOut.T_p = H / Cp_total_solid;
                pStateOut.melted_fraction = 0;
            elseif H >= H_solid_at_melt && H <= H_liquid_at_melt
                pStateOut.T_p = T_melt;
                if L_m_total > 1e-12
                    pStateOut.melted_fraction = (H - H_solid_at_melt) / L_m_total;
                else
                    pStateOut.melted_fraction = 1.0;
                end
            else
                tempState_liquid = pStateIn.copy();
                tempState_liquid.T_p = T_melt;
                Cp_total_liquid = obj.calculate_total_heat_capacity(tempState_liquid);
                
                pStateOut.T_p = T_melt + (H - H_liquid_at_melt) / Cp_total_liquid;
                pStateOut.melted_fraction = 1.0;
            end
        end

    end
end 