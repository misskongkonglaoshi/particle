%% ===============================================================
%  FlameParticleBVP_Demo.m
%  双区域 + 4 组分 + 未知火焰面
%  运行环境：MATLAB R2018b 及以上
%% ===============================================================
clear; clc; close all

%% ---------- 1. 物理/化学参数 ----------
Rp   = 1e-3;                 % 颗粒半径 [m]
Tp   = 400;                  % 颗粒温度 [K]
Tinf = 300;                  % 远场温度 [K]
rInf = 10*Rp;                % 远场半径
k    = 0.025;                % 导热系数 [W/m/K]
rhoD = 2e-5;                 % ρD [kg/m/s]
hfg  = 2.3e6;                % Mg 蒸发潜热 [J/kg]
Qf   = 50;                   % 火焰面放热 [W]（示例）
MW   = struct('Mg',24.31,'CO2',44.01,'CO',28.01,'MgO',40.31);
sigma  = 5.67e-8;            % Stefan-Boltzmann
eps    = 0.9;                % 发射率
DH_surf = -2.4e6;            % CO+Mg→MgO+C 反应热 [J/kg-CO]（示例）
DH_dep  = -1.8e6;            % MgO 沉积放热 [J/kg-MgO]（示例）

%% ---------- 2. 区域与变量 ----------
% 状态变量 z = [T, Y_Mg, Y_CO2, Y_CO, Y_MgO, dTdr, dY_Mgdr, dY_CO2dr, dY_COdr, dY_MgOdr]  (10×1)
nVar = 10;

%% ---------- 3. ODE ----------
ode = @(r,z,region,p) odefun(r,z,region,p,k,rhoD);

%% ---------- 4. 边界/接口 ----------
bc = @(yl,yr,p) bcfun(yl,yr,p,Rp,Tp,Tinf,k,rhoD,hfg,eps,sigma,Qf,DH_surf,DH_dep,MW);

%% ---------- 5. 网格与初值 ----------
rf0   = 3*Rp;                     % 火焰面初猜
mesh1 = linspace(Rp,rf0,20);
mesh2 = linspace(rf0,rInf,20);
solinit = bvpinit([mesh1 mesh2(2:end)], @guess, [rf0, 1e-6]);  % [r_f, m_dot]

%% ---------- 6. 求解 ----------
opts = bvpset('RelTol',1e-6,'AbsTol',1e-8,'Stats','on');
sol  = bvp4c(ode, bc, solinit, opts);

%% ---------- 7. 结果 ----------
rf = sol.parameters(1);
md = sol.parameters(2);
fprintf('收敛火焰半径  rf = %.4f mm\n', rf*1e3);
fprintf('总蒸发/燃烧流 m_dot = %.4e kg/s\n', md);

%% ---------- 8. 绘图 ----------
r = linspace(Rp,rInf,400);
Z = deval(sol,r);
figure; plot(r*1e3, Z(1,:), 'LineWidth',2); grid on
xlabel('r [mm]'); ylabel('T [K]'); title('温度分布');
figure; plot(r*1e3, Z(2:5,:), 'LineWidth',2); grid on
xlabel('r [mm]'); ylabel('Y_i'); legend('Mg','CO2','CO','MgO');

%% ===============================================================
% 子函数：ODE
function dz = odefun(r,z,region,k,rhoD)
    % z = [T, Y_Mg, Y_CO2, Y_CO, Y_MgO, dTdr, dY_Mgdr, dY_CO2dr, dY_COdr, dY_MgOdr]
    dz = zeros(10,1);
    % 通用对流-扩散算子
    dz(1) = z(6);  dz(2:5) = z(7:10);
    dz(6) = -2/r*z(6);        % 能量方程 (简化，无反应)
    for i = 1:4
        dz(6+i) = -2/r*z(6+i);  % 组分方程 (简化，无反应)
    end
end

%% ===============================================================
% 子函数：边界条件
function res = bcfun(yl,yr,p,Rp,Tp,Tinf,k,rhoD,hfg,eps,sigma,Qf,DH_surf,DH_dep,MW)
    r_f = p(1);  m_dot = p(2);
    % 区域 1 左/右：r=Rp, r=r_f
    y_rp   = yl{1};   y_rf_L = yr{1};
    % 区域 2 左/右：r=r_f, r=rInf
    y_rf_R = yl{2};   y_rInf = yr{2};

    res = zeros(12,1);

    % 1-5 颗粒表面 r=Rp
    res(1) = y_rp(1) - Tp;                                     % T
    % 蒸发 Spalding 关系 → Y_Mg_s
    B = exp(m_dot/(4*pi*Rp*rhoD)) - 1;
    Y_Mg_s = B/(1+B);
    res(2) = y_rp(2) - Y_Mg_s;                                 % Y_Mg
    res(3) = y_rp(8);                                          % dY_CO2/dr = 0
    res(4) = y_rp(5);                                          % Y_MgO = 0 (瞬时沉积)
    % 能量平衡：导热 + 反应放热 + 沉积放热 = 辐射 + 蒸发
    dT = y_rp(6);   dY_CO = y_rp(9);  dY_MgO = y_rp(10);
    Y_CO_s = y_rp(4);
    J_CO   = -4*pi*Rp^2*rhoD*dY_CO + m_dot*Y_CO_s;
    J_MgO  = -4*pi*Rp^2*rhoD*dY_MgO;          % 因 Y_MgO_s = 0
    Q_surf = J_CO  * DH_surf;
    Q_dep  = J_MgO * DH_dep;
    Q_rad  = 4*pi*Rp^2 * eps * sigma * (Tp^4 - Tinf^4);
    Q_evap = m_dot * hfg;
    res(5) = (-4*pi*Rp^2*k*dT) + Q_surf + Q_dep - Q_rad - Q_evap;

    % 6-7 远场 r=rInf
    res(6) = y_rInf(1) - Tinf;      % T
    res(7) = y_rInf(3) - 1.0;       % CO2 = 1

    % 8-12 火焰面 r=r_f
    res(8)  = y_rf_L(1) - y_rf_R(1);            % 温度连续
    res(9)  = y_rf_L(2);                        % Y_Mg = 0
    res(10) = y_rf_R(3);                        % Y_CO2 = 0
    % 化学计量 1:1 摩尔比
    J_Mg  = -4*pi*r_f^2*rhoD*y_rf_L(7) + m_dot*y_rf_L(2);
    J_CO2 = -4*pi*r_f^2*rhoD*y_rf_R(8) + m_dot*y_rf_R(3);
    res(11) = J_Mg/MW.Mg + J_CO2/MW.CO2;
    % 能量跃变
    dT_L = y_rf_L(6);  dT_R = y_rf_R(6);
    res(12) = (dT_R - dT_L)*k + Qf/(4*pi*r_f^2);
end

%% ===============================================================
% 子函数：初值猜测
function z = guess(r,region)
    if region==1
        T  = 400 + 1000*(r-1e-3)/(3e-3-1e-3);
        Y  = [0.2*(1-(r-1e-3)/(3e-3-1e-3)); 0; 0; 0];
    else
        T  = 1400 - 1100*(r-3e-3)/(10e-3-3e-3);
        Y  = [0; (r-3e-3)/(10e-3-3e-3); 1-(r-3e-3)/(10e-3-3e-3); 0];
    end
    z = [T; Y; zeros(5,1)];
end