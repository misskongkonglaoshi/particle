classdef Parameters < handle
    % Parameters 仿真参数管理类
    
    properties
        % 物理常数
        R_u = 8.314;         % 通用气体常数 (J/(mol·K))
        sigma = 5.67e-8;     % Stefan-Boltzmann常数 (W/(m²·K⁴))
        
        % 材料物理参数
        materials = struct(...
            'Mg', struct(...
                'density', 1740, ...        % 密度 (kg/m^3)
                'heat_capacity', 1020, ...  % 比热容 (J/(kg·K))
                'latent_heat', 8.7e5, ...  % 熔化潜热 (J/kg)
                'L_evap_Mg', 5.25e6, ...    % 蒸发潜热 (J/kg)
                'molar_mass', 24.305e-3, ...% 摩尔质量 (kg/mol)
                'melting_point', 923, ...   % 熔点 (K)
                'Hf298', 147100, ...        % 气相标准生成焓 (J/mol)
                'ignition_temp', 1366 ...   % 点火温度 (K)
            ), ...
            'MgO', struct(...
                'density', 3580, ...        % 密度 (kg/m^3)
                'heat_capacity', 937, ...   % 比热容 (J/(kg·K))
                'molar_mass', 40.304e-3, ...% 摩尔质量 (kg/mol)
                'thermal_conductivity', 30, ...  % 热导率 (W/(m·K))
                'Hf298', -19400, ...        % 气相标准生成焓 (J/mol)
                'diffusivity', 1e-12 ...    % 扩散系数 (m^2/s)
            ), ...
            'C', struct(...
                'density', 2267, ...
                'heat_capacity', 709, ...
                'molar_mass', 12.011e-3 ...
            ), ...
            'CO', struct(...
                'molar_mass', 28.010e-3, ...
                'Hf298', -110500 ...        % 标准生成焓 (J/mol)
            ), ...
            'CO2', struct(...
                'molar_mass', 44.010e-3, ...
                'Hf298', -393500 ...        % 标准生成焓 (J/mol)
            ) ...
        );
        
        % 蒸气压计算参数 
        antoine_coeffs_mg = struct(...
            'A', 7.378, ...
            'B', 1605, ...
            'C', 211 ...
        );

        % 颗粒物理参数
        initial_diameter = 50e-6;  % 初始颗粒直径 (m)
        initial_temperature = 300;  % 初始颗粒温度 (K)
        initial_oxide_thickness = 5e-7;  % 初始氧化层厚度 (m)，默认为100nm
        
        % 环境参数
        ambient_temperature = 1500;  % 环境温度 (K)
        ambient_pressure = 101325;   % 环境压力 (Pa)
        emissivity = 0.8;           % 发射率
        h_conv = 50;                % 对流换热系数 (W/(m²·K))
        
        % 气体组成 (摩尔分数)
        ambient_gas_composition = struct(...
            'CO2', 1.0, ...  % 纯CO2环境
            'O2', 0.0, ...
            'N2', 0.0 ...
        );
        
        % 仿真控制参数
        time_step = 1e-6;          % 时间步长 (s)
        total_time = 1;          % 总仿真时间 (s)
        t_combustion = 0.05;       % 气相燃烧求解时长 (s)
        output_interval = 1000;     % 输出间隔（步数）
        
        % 气相计算初始设置火焰面 1.5倍初始直径
        flam_thickness = 2.5e-5;      % 气膜厚度 (m)
        
        % 气体物理属性 (示例值, 应根据实际情况调整)
        gas_properties = struct(...
            'k_gas', 0.1, ...       % 气体混合物平均热导率 (W/(m·K))
            'rho_gas', 0.2, ...     % 气体混合物平均密度 (kg/m^3)
            'Cp_gas', 1500 ...      % 气体混合物平均比热 (J/(kg·K))
        );
        
        % 求解器配置
        solver_options = struct(...
            'RelTol', 1e-6, ...
            'AbsTol', 1e-8, ...
            'MaxStep', 1e-5 ...
        );
    end
    
    methods
        function obj = Parameters()
            % 构造函数：初始化参数对象
        end
        
        function validate(obj)
            % 验证参数的有效性
            assert(obj.initial_diameter > 0, '初始直径必须大于0');
            assert(obj.initial_temperature > 0, '初始温度必须大于0');
            assert(obj.materials.Mg.melting_point > 0, '熔点必须大于0');
            assert(obj.materials.Mg.ignition_temp > obj.materials.Mg.melting_point, ...
                '点火温度必须高于熔点');
            assert(obj.ambient_temperature > 0, '环境温度必须大于0');
            assert(obj.ambient_pressure > 0, '环境压力必须大于0');
            assert(obj.emissivity > 0 && obj.emissivity <= 1, '发射率必须在0到1之间');
            assert(obj.time_step > 0, '时间步长必须大于0');
            assert(obj.total_time > 0, '总仿真时间必须大于0');
            assert(obj.output_interval > 0, '输出间隔必须大于0');
            assert(obj.initial_oxide_thickness >= 0, '初始氧化层厚度必须大于等于0');
        end
    end
end 