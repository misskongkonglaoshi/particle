classdef ThermoReader < handle
    % ThermoReader 读取CHEMKIN-III格式的热力学数据文件并计算性质
    
    properties
        thermo_file
        params
        species_data
        Ru = 8.31446261815324; % J/(mol·K)
    end
    
    methods
        function obj = ThermoReader(thermo_file, params)
            if nargin > 0
                obj.thermo_file = thermo_file;
                if nargin > 1
                    obj.params = params;
                end
                obj.species_data = struct();
                obj.read_thermo_file();
            end
        end
        
        function read_thermo_file(obj)
            fid = fopen(obj.thermo_file, 'r');
            if fid == -1
                error('无法打开文件: %s', obj.thermo_file);
            end
            
            % 跳过文件头直到THERMO关键字
            line = fgetl(fid);
            while ischar(line) && ~strncmpi(strtrim(line), 'THERMO', 6)
                line = fgetl(fid);
            end
            
            % 读取并忽略全局温度行
            fgetl(fid);
            
            % 循环读取物种数据
            while true
                line1 = fgetl(fid);
                if ~ischar(line1) || strncmpi(strtrim(line1), 'END', 3)
                    break;
                end
                
                line2 = fgetl(fid);
                line3 = fgetl(fid);
                line4 = fgetl(fid);

                try
                    % --- 开始修改: 使用sprintf填充字符串以确保长度 ---
                    line2 = sprintf('%-75s', line2);
                    line3 = sprintf('%-75s', line3);
                    line4 = sprintf('%-60s', line4); % 第4行有4个系数，需要60个字符
                    % --- 结束修改 ---

                    % 解析第一行
                    species_name = strtrim(line1(1:18));
                    % --- 开始修改: 扩大温度读取范围 ---
                    temp_ranges_str = sscanf(line1(46:73), '%f');
                    % --- 结束修改 ---
                    T_low = temp_ranges_str(1);
                    T_high = temp_ranges_str(2);
                    T_mid = temp_ranges_str(3);

                    % --- 开始修改: 使用健壮的固定宽度解析方法 ---
                    coeffs_str_cell = {};
                    % 第2行: 5个系数
                    for i=1:5
                        coeffs_str_cell{end+1} = line2((i-1)*15+1 : i*15);
                    end
                    % 第3行: 5个系数
                    for i=1:5
                        coeffs_str_cell{end+1} = line3((i-1)*15+1 : i*15);
                    end
                    % 第4行: 4个系数
                    for i=1:4
                        coeffs_str_cell{end+1} = line4((i-1)*15+1 : i*15);
                    end
                    coeffs_vals = str2double(coeffs_str_cell);
                    % --- 结束修改 ---
                    
                    if length(coeffs_vals) ~= 14
                        warning('物种 %s 的系数数量不为14，跳过。', species_name);
                        continue;
                    end
                    
                    % NASA 7系数多项式
                    high_T_coeffs = coeffs_vals(1:7);
                    low_T_coeffs = coeffs_vals(8:14);
                    
                    obj.species_data.(species_name) = struct(...
                        'T_low', T_low, 'T_high', T_high, 'T_mid', T_mid, ...
                        'high_T_coeffs', high_T_coeffs, ...
                        'low_T_coeffs', low_T_coeffs);
                catch ME
                    warning('解析物种 %s 时出错: %s', species_name, ME.message);
                end
            end
            
            fclose(fid);
        end
        
        function coeffs = get_coeffs(obj, species, T)
            if ~isfield(obj.species_data, species)
                error('缺少物种 %s 的热力学数据', species);
            end
            data = obj.species_data.(species);
            if T > data.T_mid
                coeffs = data.high_T_coeffs;
            else
                coeffs = data.low_T_coeffs;
            end
        end

        function Cp = calculate_Cp(obj, species, T)
            % 计算定压比热容 (J/mol·K)
            a = obj.get_coeffs(species, T);
            T_vec = [1; T; T.^2; T.^3; T.^4];
            Cp = (a(1:5) * T_vec) * obj.Ru;
        end
        
        function H = calculate_H(obj, species, T)
            % 计算焓 (J/mol)
            a = obj.get_coeffs(species, T);
            T_vec = [1; T/2; T.^2/3; T.^3/4; T.^4/5; 1./T];
            H = (a(1:6) * T_vec) * obj.Ru * T;
        end
        
        function S = calculate_S(obj, species, T)
            % 计算熵 (J/mol·K)
            a = obj.get_coeffs(species, T);
            T_vec = [log(T); T; T.^2/2; T.^3/3; T.^4/4];
            S = (a(1:5) * T_vec + a(7)) * obj.Ru;
        end
        
        function p_sat = get_p_sat(obj, T_kelvin)
            if isempty(obj.params)
                 error('需要params对象来计算饱和蒸气压');
            end
            T_celsius = T_kelvin - 273.15;
            A = obj.params.antoine_coeffs_mg.A;
            B = obj.params.antoine_coeffs_mg.B;
            C = obj.params.antoine_coeffs_mg.C;
            log10_p_mmHg = A - B / (T_celsius + C);
            p_mmHg = 10^log10_p_mmHg;
            p_sat = p_mmHg * 133.322;
        end
        
        function has_data = has_species(obj, species)
            has_data = isfield(obj.species_data, species);
        end

        function value = get_property(obj, species, property_name)
            % 通用属性获取函数
            % 优先从thermo.dat计算，然后回退到params.m
            
            % 尝试从CHEMKIN数据计算
            if obj.has_species(species)
                switch lower(property_name)
                    case 'hf298' % 298.15K下的标准生成焓
                        value = obj.calculate_H(species, 298.15);
                        return;
                    case 'cp'
                        % 注意: Cp是温度依赖的, 这里只为示例，实际应调用calculate_Cp
                        % value = obj.calculate_Cp(species, T); 
                        % 此处需要一个温度参数, 因此不直接处理
                end
            end
            
            % 如果无法计算, 尝试从params回退
            if ~isempty(obj.params) && isfield(obj.params.materials, species) ...
               && isfield(obj.params.materials.(species), property_name)
                value = obj.params.materials.(species).(property_name);
                return;
            end
            
            % 如果都找不到，则报错
            error('ThermoReader:PropertyNotFound', ...
                  '无法在thermo.dat或Parameters.m中找到物种 "%s" 的属性 "%s"。', ...
                  species, property_name);
        end
    end
end 