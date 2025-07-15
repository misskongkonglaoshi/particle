classdef ParticleModel < handle
    % ParticleModel 颗粒模型
    % 作为颗粒状态的唯一权威，包含并管理一个ParticleState对象
    
    properties
        particleState  % 颗粒状态对象 (ParticleState)
        params         % 参数对象
    end
    
    methods
        function obj = ParticleModel(params)
            % 构造函数
            %
            % 输入:
            %   params: 参数对象
            if nargin > 0
                obj.params = params;
                % 创建并持有ParticleState实例
                obj.particleState = ParticleState(params);
            end
        end
        
        function stage = determineStage(obj)
            % 根据当前颗粒状态确定所处的阶段 (引入容差以提高鲁棒性)
            
            TOL = 1e-6; % 定义一个小的容差来处理浮点数比较问题
            T_p = obj.particleState.T_p;
            T_melt = obj.params.materials.Mg.melting_point;
            T_ignite = obj.params.materials.Mg.ignition_temp;
            ps = obj.particleState;
            
            if T_p < T_melt - TOL
                stage = 'preheating';
            elseif abs(T_p - T_melt) < TOL && ps.melted_fraction < 1
                stage = 'melting';
            elseif T_p > T_melt - TOL && ps.melted_fraction >= 1 && T_p < T_ignite
                stage = 'liquid_heating';
            elseif T_p >= T_ignite
                stage = 'vaporization';
            else
                % 默认或未知状态处理, 增加诊断信息
                fprintf('警告: 无法确定颗粒阶段。 T_p=%.4f, melted_fraction=%.4f\n', T_p, ps.melted_fraction);
                stage = 'unknown';
            end
        end
        
        function update_oxide_layer(obj, consumed_mg_mass)
            % 更新氧化层厚度
            %
            % 输入:
            %   consumed_mg_mass: 消耗的Mg质量 (kg)
            
            % 更新颗粒状态
            obj.particleState.update_state_from_reaction(-consumed_mg_mass);
        end
        
        function cp = getHeatCapacity(obj, physicalModel)
            % 获取颗粒总热容
            cp = physicalModel.calculate_total_heat_capacity(obj.particleState);
        end
        
        function mass = getMass(obj)
            % 获取颗粒总质量
            mass = obj.particleState.m_mg + obj.particleState.m_mgo + obj.particleState.m_c;
        end
        
        function r = getRadius(obj)
            % 获取颗粒总半径
            r = obj.particleState.r_p;
        end
    end
end 