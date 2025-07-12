classdef Visualization
    % Visualization 结果可视化类
    
    methods (Static)
        function plotResults(results)
            % 绘制仿真结果
            %
            % 输入:
            %   results: 仿真结果结构体
            
            % 创建一个新的图形窗口
            figure('Name', 'Mg颗粒燃烧仿真结果', 'Position', [100, 100, 1200, 800]);
            
            % 子图1: 温度历史
            subplot(3, 2, 1);
            plot(results.time, results.temperature, 'LineWidth', 1.5);
            title('温度历史');
            xlabel('时间 (s)');
            ylabel('温度 (K)');
            grid on;
            
            % 子图2: 熔融分数历史
            subplot(3, 2, 2);
            plot(results.time, results.melted_fraction, 'LineWidth', 1.5);
            title('熔融分数历史');
            xlabel('时间 (s)');
            ylabel('熔融分数');
            ylim([0 1]);
            grid on;
            
            % 子图3: 阶段历史
            subplot(3, 2, 3);
            
            % 将阶段名称转换为数值以便绘图
            [unique_stages, ~, stage_indices] = unique(results.stage);
            
            % 绘制阶段历史
            plot(results.time, stage_indices, 'LineWidth', 1.5);
            title('阶段历史');
            xlabel('时间 (s)');
            yticks(1:length(unique_stages));
            yticklabels(unique_stages);
            grid on;
            
            % 子图4: 氧化层厚度历史
            subplot(3, 2, 4);
            if isfield(results, 'oxide_thickness')
                % 将氧化层厚度转换为微米
                plot(results.time, results.oxide_thickness * 1e6, 'LineWidth', 1.5);
                title('氧化层厚度历史');
                xlabel('时间 (s)');
                ylabel('氧化层厚度 (µm)');
                grid on;
            else
                text(0.5, 0.5, '氧化层数据不可用', 'HorizontalAlignment', 'center');
                title('氧化层厚度历史');
                axis off;
            end
            
            % 子图5: 质量历史
            subplot(3, 2, 5);
            hold on;
            if isfield(results, 'mass_mg')
                plot(results.time, results.mass_mg * 1e9, 'LineWidth', 1.5, 'DisplayName', 'Mg');
            end
            if isfield(results, 'mass_mgo')
                plot(results.time, results.mass_mgo * 1e9, 'LineWidth', 1.5, 'DisplayName', 'MgO');
            end
            if isfield(results, 'mass_c')
                plot(results.time, results.mass_c * 1e9, 'LineWidth', 1.5, 'DisplayName', 'C');
            end
            title('组分质量历史');
            xlabel('时间 (s)');
            ylabel('质量 (ng)');
            legend('Location', 'best');
            grid on;
            hold off;
            
            % 子图6: 总质量和颗粒直径
            subplot(3, 2, 6);
            yyaxis left;
            if isfield(results, 'mass_mg') && isfield(results, 'mass_mgo')
                total_mass = results.mass_mg + results.mass_mgo;
                if isfield(results, 'mass_c')
                    total_mass = total_mass + results.mass_c;
                end
                plot(results.time, total_mass * 1e9, 'LineWidth', 1.5);
                ylabel('总质量 (ng)');
            end
            
            yyaxis right;
            if isfield(results, 'oxide_thickness')
                % 假设有metal_radius数据或者可以从其他数据计算
                % 这里简化，假设metal_radius可以从质量计算
                if isfield(results, 'mass_mg')
                    % 计算金属核心半径 (m)
                    density_mg = 1740;  % Mg的密度 (kg/m^3)
                    metal_radius = (3 * results.mass_mg ./ (4 * pi * density_mg)).^(1/3);
                    
                    % 计算总半径 (m)
                    total_radius = metal_radius + results.oxide_thickness;
                    
                    % 绘制总直径 (µm)
                    plot(results.time, 2 * total_radius * 1e6, 'LineWidth', 1.5, 'Color', 'r');
                    ylabel('颗粒直径 (µm)');
                end
            end
            
            title('总质量和颗粒直径');
            xlabel('时间 (s)');
            grid on;
            
            % 调整子图间距
            set(gcf, 'Color', 'w');
            tight_spacing = 0.04;
            loose_spacing = 0.12;
            
            % 应用紧凑布局
            % subplot_tight(3, 2, 1, [tight_spacing, loose_spacing]);
        end
        
        function saveResults(results, filename)
            % 保存仿真结果到MAT文件
            %
            % 输入:
            %   results: 仿真结果结构体
            %   filename: 保存文件名
            
            save(filename, 'results');
            fprintf('结果已保存到文件: %s\n', filename);
        end
        
        function exportToCSV(results, filename)
            % 导出仿真结果到CSV文件
            %
            % 输入:
            %   results: 仿真结果结构体
            %   filename: 保存文件名
            
            % 准备数据表
            data_table = table(results.time, results.temperature, results.melted_fraction, ...
                'VariableNames', {'Time', 'Temperature', 'MeltedFraction'});
            
            % 添加阶段信息
            data_table.Stage = categorical(results.stage);
            
            % 添加氧化层厚度（如果存在）
            if isfield(results, 'oxide_thickness')
                data_table.OxideThickness = results.oxide_thickness;
            end
            
            % 添加质量信息（如果存在）
            if isfield(results, 'mass_mg')
                data_table.MassMg = results.mass_mg;
            end
            
            if isfield(results, 'mass_mgo')
                data_table.MassMgO = results.mass_mgo;
            end
            
            if isfield(results, 'mass_c')
                data_table.MassC = results.mass_c;
            end
            
            % 写入CSV文件
            writetable(data_table, filename);
            fprintf('结果已导出到CSV文件: %s\n', filename);
        end
    end
end 