classdef ThermoReader < handle
    % ThermoReader 读取热力学数据文件并计算热力学性质
    
    properties
        thermo_file     % thermo.dat文件路径
        species_data    % 物种数据存储
        Ru = 8.314      % 通用气体常数 (J/mol·K)
    end
    
    methods
        function obj = ThermoReader(thermo_file)
            % 构造函数：初始化热力学数据读取器
            % 输入:
            %   thermo_file: thermo.dat文件的路径
            
            obj.thermo_file = thermo_file;
            obj.species_data = struct();
            obj.read_thermo_file();
        end
        
        function read_thermo_file(obj)
            % 读取thermo.dat文件并解析数据
            
            % 读取文件
            fid = fopen(obj.thermo_file, 'r');
            if fid == -1
                error('无法打开文件: %s', obj.thermo_file);
            end
            
            % 逐行读取数据
            while ~feof(fid)
                % 读取物种信息行
                line = fgetl(fid);
                if isempty(strtrim(line)) || line(1) == '!'
                    continue;
                end
                
                % 解析物种名称和相态
                species_name = strtrim(line(1:18));
                if length(line) >= 45
                    phase = line(45);
                else
                    phase = 'G';  % 默认为气体
                end
                
                try
                    % 读取三行系数
                    line1 = fgetl(fid);
                    cp_coeffs = sscanf(line1, '%f');
                    
                    line2 = fgetl(fid);
                    h_coeffs = sscanf(line2, '%f');
                    
                    line3 = fgetl(fid);
                    s_coeffs = sscanf(line3, '%f');
                    
                    % 存储物种数据
                    obj.species_data.(species_name) = struct(...
                        'phase', phase, ...
                        'cp_coeffs', cp_coeffs, ...
                        'h_coeffs', h_coeffs, ...
                        's_coeffs', s_coeffs);
                catch ME
                    warning('解析%s的热力学数据时出错：%s', species_name, ME.message);
                end
            end
            
            fclose(fid);
        end
        
        function Cp = calculate_Cp(obj, species, T)
            % 计算定压比热容
            % 输入:
            %   species: 物种名称
            %   T: 温度 (K)
            % 输出:
            %   Cp: 定压比热容 (J/mol·K)
            
            if ~isfield(obj.species_data, species)
                error('缺少物种 %s 的热力学系数，无法计算比热容', species);
            end
            
            a = obj.species_data.(species).cp_coeffs;
            T_powers = [1, T, T^2, T^3, T^4];
            Cp = sum(a .* T_powers) * obj.Ru;
        end
        
        function H = calculate_H(obj, species, T)
            % 计算焓值
            % 输入:
            %   species: 物种名称
            %   T: 温度 (K)
            % 输出:
            %   H: 焓值 (J/mol)
            
            if ~isfield(obj.species_data, species)
                error('缺少物种 %s 的热力学系数，无法计算焓值', species);
            end
            
            a = obj.species_data.(species).h_coeffs;
            T_powers = [1, 1/2*T, 1/3*T^2, 1/4*T^3, 1/5*T^4];
            H = sum(a .* T_powers) * obj.Ru * T;
        end
        
        function S = calculate_S(obj, species, T)
            % 计算熵
            % 输入:
            %   species: 物种名称
            %   T: 温度 (K)
            % 输出:
            %   S: 熵 (J/mol·K)
            
            if ~isfield(obj.species_data, species)
                error('缺少物种 %s 的热力学系数，无法计算熵', species);
            end
            
            a = obj.species_data.(species).s_coeffs;
            T_powers = [log(T), T, T^2/2, T^3/3, T^4/4];
            S = sum(a .* T_powers) * obj.Ru;
        end
        
        function has_data = has_species(obj, species)
            % 检查是否存在指定物种的热力学数据
            % 输入:
            %   species: 物种名称
            % 输出:
            %   has_data: 是否存在该物种的数据
            
            has_data = isfield(obj.species_data, species);
        end
    end
end 