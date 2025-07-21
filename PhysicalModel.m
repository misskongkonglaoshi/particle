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
            n_c = particleState.m_c / obj.params.materials.C.molar_mass;
            n_total = n_mg + n_mgo + n_c;
            
            if n_total > 1e-20
                % 计算摩尔分数
                x_mg = n_mg / n_total;
                x_mgo = n_mgo / n_total;
                x_c = n_c / n_total;
                % 计算当前温度下的摩尔比热容
                cp_mg_mol = obj.calc_cp_solid_mg(particleState.T_p);
                cp_mgo_mol = obj.calc_cp_solid_mgo(particleState.T_p);
                cp_c_mol = obj.calc_cp_solid_c(particleState.T_p);  
                % 按摩尔分数加权计算总的摩尔比热容
                cp_mol = x_mg * cp_mg_mol + x_mgo * cp_mgo_mol + x_c * cp_c_mol;
                
                % 转换为总热容
                cp_total = cp_mol * n_total;
            else
                cp_total = 0;
            end
        end
        
        function cp = calc_cp_solid_mg(~, T)
            % 计算固态镁的摩尔比热容 (J/(mol·K))
            % 使用线性多项式 Cp = a + b*T
            t = T /1000;
            if T < 923
                a = 26.541;  % J/(mol·K)
                b = -1.533;
                c = 8.062;
                d = 0.572;
                e = -0.174 % J/(mol·K²)
            elseif T < 1363
                a = 34.309;  % J/(mol·K)
                b = -7.471e-10;
                c = 6.146e-10;
                d = -1.59e-10;
                e = -1.152e-11 % J/(mol·K²)
            end
            cp = a + b * t + c * t^2 + d * t^3 + e /t^2;
        end
        
        function cp = calc_cp_solid_mgo(~, T)
            % 计算氧化镁的摩尔比热容 (J/(mol·K))
            t = T /1000;
            a = 47.26;  % J/(mol·K)
            b = 5.682;
            c = -0.873;
            d = 0.104;
            e = -1.054 ; % J/(mol·K²)
            cp = a + b * t + c * t^2 + d * t^3 + e /t^2;
        end
        
        function cp = calc_cp_solid_c(~, T)
            % 计算固态碳的摩尔比热容 (J/(mol·K))
            cp = 8.6186;
        end
        
        function cp = calculate_specific_heat(obj, species, T)
            % 使用thermo_reader计算特定温度下的摩尔比热容
            %
            % 输入:
            %   species: 物质名称
            %   T: 温度 (K)
            % 输出:
            %   cp: 摩尔比热容 (J/(mol·K))
            
            % 根据不同物种调用对应的比热计算函数
            switch species
                case 'Mg'
                    cp = obj.calc_cp_solid_mg(T);
                case 'MgO'
                    cp = obj.calc_cp_solid_mgo(T);
                case 'C'
                    % 注意：这里返回的是摩尔比热容，需要转换
                    molar_mass_c = 0.012; % kg/mol
                    cp = obj.calc_cp_solid_c(T) * molar_mass_c;
                otherwise
                    % 对于气相物种，尝试使用thermo_reader
                    if ~isempty(obj.thermo_reader) && obj.thermo_reader.has_species(species)
                        % 使用NASA 7系数多项式计算
                        cp = obj.thermo_reader.calculate_Cp(species, T);
                    else
                        % 回退到参数中的常数值
                        cp = obj.params.materials.(species).heat_capacity * obj.params.materials.(species).molar_mass;
                    end
            end
            
            % 防御性修复: 确保返回的是一个标量
            if ~isscalar(cp)
                cp = sum(cp);
            end
        end
        
        function cp = calc_cp_gas(obj, species, T)
            % 计算气体的摩尔比热容 (J/(mol·K))
            % 使用NASA 7系数多项式
            
            if ~isempty(obj.thermo_reader) && obj.thermo_reader.has_species(species)
                cp = obj.thermo_reader.calculate_Cp(species, T);
            else
                % 如果没有数据，使用默认值
                switch species
                    case 'CO2'
                        % 简化的CO2比热多项式
                        cp = 22.26 + 5.981e-2*T - 3.501e-5*T^2 + 7.469e-9*T^3;
                    case 'CO'
                        % 简化的CO比热多项式
                        cp = 25.56 + 6.096e-3*T + 4.054e-6*T^2 - 2.671e-9*T^3;
                    case 'O2'
                        % 简化的O2比热多项式
                        cp = 25.48 + 1.520e-2*T - 7.155e-6*T^2 + 1.312e-9*T^3;
                    case 'N2'
                        % 简化的N2比热多项式
                        cp = 28.58 - 3.330e-3*T + 1.035e-5*T^2 - 3.729e-9*T^3;
                    otherwise
                        % 默认值
                        cp = 29.1; % 近似值 (J/(mol·K))
                end
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
            props = obj.params.gas_properties;
            
            % 计算温度依赖的热导率
            props.k_gas = obj.calc_k_gas_mixture(T_gas);
            
            % 计算温度依赖的比热容
            props.Cp_gas = obj.calc_cp_gas_mixture(T_gas);
        end
        
        function k_mix = calc_k_gas_mixture(obj, T)
            % 计算混合气体的热导率 (W/(m·K))
            gas_comp = obj.params.ambient_gas_composition;
            gas_fields = fieldnames(gas_comp);
            
            k_mix = 0;
            total_fraction = 0;
            
            for i = 1:length(gas_fields)
                gas_name = gas_fields{i};
                mole_fraction = gas_comp.(gas_name);
                
                if mole_fraction > 0
                    k_species = obj.calc_k_gas_species(gas_name, T);
                    k_mix = k_mix + mole_fraction * k_species;
                    total_fraction = total_fraction + mole_fraction;
                end
            end
            
            if total_fraction > 1e-9
                k_mix = k_mix / total_fraction;
            else
                k_mix = 0.026; % 默认值 (W/(m·K))
            end
        end
        
        function k = calc_k_gas_species(~, species, T)
            % 计算单个气体组分的热导率 (W/(m·K))
            switch species
                case 'CO2'
                    % CO2热导率多项式系数
                    a = -0.01183;
                    b = 1.0174e-4;
                    c = -2.2242e-8;
                    k = a + b * T + c * T^2;
                case 'O2'
                    % O2热导率多项式系数
                    a = 0.0074;
                    b = 7.0e-5;
                    k = a + b * T;
                    fprintf('出现O2,更新氧气热导率计算');
                case 'N2'
                    % N2热导率多项式系数
                    a = 0.0066;
                    b = 6.5e-5;
                    k = a + b * T;
                    fprintf('出现O2,更新氧气热导率计算');
                case 'CO'
                    % CO热导率多项式系数
                    a = 0.0059;
                    b = 6.3e-5;
                    k = a + b * T;
                    fprintf('出现O2,更新氧气热导率计算');
                case 'Mg'
                    % 气态Mg热导率（简化）
                    fprintf('出现O2,更新氧气热导率计算');
                    k = 0.01;
                otherwise
                    k = 0.026; % 默认值
            end
        end
        
        function cp_mix = calc_cp_gas_mixture(obj, T)
            % 计算混合气体的质量比热容 (J/(kg·K))
            gas_comp = obj.params.ambient_gas_composition;
            gas_fields = fieldnames(gas_comp);
            
            cp_mol_mix = 0;
            total_fraction = 0;
            
            for i = 1:length(gas_fields)
                gas_name = gas_fields{i};
                mole_fraction = gas_comp.(gas_name);
                
                if mole_fraction > 0
                    cp_mol_species = obj.calc_cp_gas(gas_name, T);
                    cp_mol_mix = cp_mol_mix + mole_fraction * cp_mol_species;
                    total_fraction = total_fraction + mole_fraction;
                end
            end
            
            if total_fraction > 1e-9
                cp_mol_mix = cp_mol_mix / total_fraction;
                
                % 转换为质量比热容
                M_mix = obj.get_ambient_mixture_molar_mass();
                cp_mix = cp_mol_mix / M_mix;
            else
                cp_mix = 1000; % 默认值 (J/(kg·K))
            end
        end

        function M_mix = get_ambient_mixture_molar_mass(obj)
            % 根据环境气体组分计算加权平均的摩尔质量 (kg/mol)
            M_mix = 0;
            total_fraction = 0;
            
            gas_fields = fieldnames(obj.params.ambient_gas_composition);
            for i = 1:length(gas_fields)
                gas_name = gas_fields{i};
                mole_fraction = obj.params.ambient_gas_composition.(gas_name);
                
                if mole_fraction > 0
                    if isfield(obj.params.materials, gas_name)
                        molar_mass_gas = obj.params.materials.(gas_name).molar_mass;
                        M_mix = M_mix + mole_fraction * molar_mass_gas;
                        total_fraction = total_fraction + mole_fraction;
                    else
                        warning('PhysicalModel:gasNotFound', ...
                            '在 materials 结构体中未找到气体 "%s" 的数据。', gas_name);
                    end
                end
            end
            
            if total_fraction > 1e-9 % 避免除零
                M_mix = M_mix / total_fraction;
            else
                warning('PhysicalModel:noGasComposition', ...
                    '环境气体组分总和为零，无法计算平均摩尔质量。');
                M_mix = 28.97e-3; % 回退到空气的近似值
            end
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