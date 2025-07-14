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
        end
        
        function props = get_gas_properties(obj, T_gas)
            % 获取在特定温度下的气体属性
            % TODO: 实现更复杂的温度依赖模型
            % 目前, 返回在Parameters中定义的常数值
            props = obj.params.gas_properties;
        end

    end
end 