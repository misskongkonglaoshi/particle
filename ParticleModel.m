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
            % 根据当前颗粒状态确定所处的阶段
            
            T = obj.particleState.T_p;
            T_melt = obj.params.materials.Mg.melting_point;
            T_ignite = obj.params.materials.Mg.ignition_temp;
            melt_frac = obj.particleState.melted_fraction;

            if T >= T_ignite
                stage = 'vaporization';
            elseif T < T_melt
                stage = 'preheating';
            elseif T >= T_melt && T < T_ignite
                if melt_frac < 1.0
                    stage = 'melting';
                else % melt_frac == 1.0
                    stage = 'liquid_heating';
                end
            else
                % 兜底情况，理论上不应到达
                stage = 'preheating';
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