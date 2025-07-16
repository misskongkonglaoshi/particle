function rate_info = solve_vaporization_bvp(bvp_deps)
% SOLVE_VAPORIZATION_BVP - V28 最终修正: 使用全局变量
%
% 结论: 在深层调用栈中, bvp4c/bvparguments 无法处理任何形式的
% 闭包或嵌套上下文。唯一的解决方法是使用全局变量传递依赖，并向
% bvp4c 传递纯函数句柄。
global GLOBAL_BVP_DEPS;

    % --- 1. 从依赖项中解包变量 ---
    pState = bvp_deps.pState;
    params = bvp_deps.params;
    
    r_p = pState.r_p;
    r_inf = r_p + params.flam_thickness;

    % --- 2. 对未知参数提供合理的初始猜测 ---
    r_f_guess = r_p * 1.5;
    m_dot_guess = 4 * pi * r_p * 1e-4;
    p_guess = [r_f_guess, m_dot_guess];

    % --- 3. 准备BVP的初始结构 (手动构建, 绕过bvpinit) ---
    xinit = [linspace(r_p, r_f_guess, 15), linspace(r_f_guess, r_inf, 15)];

    % 手动计算初始猜测矩阵 yinit
    y_guess_matrix = zeros(10, numel(xinit));
    for i = 1:numel(xinit)
        current_r = xinit(i);
        if current_r <= r_f_guess
            region = 1;
        else
            region = 2;
        end
        y_guess_matrix(:, i) = initial_guess_func(current_r, region, pState, params, r_f_guess, r_inf);
    end
    
    % 手动构建 solinit 结构体 (绕过 bvpinit)
    solinit.x = xinit;
    solinit.y = y_guess_matrix;
    solinit.parameters = p_guess; % <--- 关键补充：将未知参数的猜测提供给求解器
    
    % --- 4. 准备ODE和BC函数的依赖项 ---
    % 将依赖存入全局变量
    GLOBAL_BVP_DEPS = bvp_deps;

    % --- 5. 创建纯函数句柄 ---
    ode_handle = @ode_func;
    bc_handle = @bc_func;

    % --- 6. 求解 ---
    options = bvpset('NMax', 6000, 'RelTol', 1e-4, 'AbsTol', 1e-7);
    try
        sol = bvp4c(ode_handle, bc_handle, solinit, options);
    catch ME
        % 确保在出错时也清理全局变量
        clear global GLOBAL_BVP_DEPS;
        rethrow(ME);
    end
    
    % 清理全局变量
    clear global GLOBAL_BVP_DEPS;

    % --- 7. 从解中提取结果并打包 ---
    r_f_sol = sol.parameters(1);
    m_dot_sol = sol.parameters(2);
    
    y_surf = deval(sol, r_p);
    dYco_dr_surf = y_surf(9);
    Yco_surf = y_surf(4);
    rho_D_avg = params.rho_D_gas;
    A_p = 4 * pi * r_p^2;
    m_dot_co_surf = -(A_p * rho_D_avg) * dYco_dr_surf + m_dot_sol * Yco_surf;
    
    rate_info = package_rate_info_func(m_dot_sol, m_dot_co_surf, sol, r_f_sol, params);

    % =========================================================================
    %  BVP核心方程 (作为此文件内的嵌套函数)
    % =========================================================================

end

% =========================================================================
%  独立的局部函数 (不嵌套)
% =========================================================================

function dzdr = ode_func(r, z, region, p)
    % BVP的常微分方程 (使用全局变量)
    global GLOBAL_BVP_DEPS;
    deps = GLOBAL_BVP_DEPS;

    m_dot_local = p(2);
    
    T = z(1,:);
    props = deps.physicalModel.get_gas_properties(T);
    k_gas = props.k_gas;
    cp_gas = props.Cp_gas;
    rho_D = deps.params.rho_D_gas;
    
    dzdr = zeros(size(z));
    
    dzdr(1:5,:) = z(6:10,:);
    dzdr(6,:) = (m_dot_local * cp_gas ./ (4 * pi * r.^2 .* k_gas)) .* z(6,:) - (2./r) .* z(6,:);
    common_term = m_dot_local ./ (4 * pi * r.^2 * rho_D);
    dzdr(7:10,:) = common_term .* z(7:10,:) - (2./r) .* z(7:10,:);
end

function res = bc_func(yl, yr, p)
    % BVP的边界/界面条件 (使用全局变量)
    % 版本4 (最终版): 基于完整的物理和数学分析构建的所有22个条件
    global GLOBAL_BVP_DEPS;
    deps = GLOBAL_BVP_DEPS;

    % --- 解包参数和状态 ---
    r_f_local = p(1);
    m_dot_local = p(2);
    pState_local = deps.pState;
    r_p_local = pState_local.r_p;

    y_at_rp       = yl(:,1);
    y_at_rf_left  = yr(:,1);
    y_at_rf_right = yl(:,2);
    y_at_rinf     = yr(:,2);

    res = zeros(22,1);
    
    % =========================================================================
    % 组1: 外边界条件 (共9个)
    % =========================================================================

    % --- 条件 1-4: 颗粒表面 r = r_p (4个) ---
    % 1. 温度等于颗粒沸点
    res(1) = y_at_rp(1) - pState_local.T_p;
    % 2. Mg蒸气浓度由蒸发速率决定
    B_m = exp(m_dot_local / (4*pi*r_p_local*deps.params.rho_D_gas)) - 1;
    res(2) = y_at_rp(2) - B_m/(1+B_m);
    % 3. CO2在颗粒表面无通量梯度
    res(3) = y_at_rp(8); % dY_co2/dr = 0
    % 4. 颗粒表面能量守恒
    T_surf = y_at_rp(1); dTdr_surf = y_at_rp(6); dYco_dr_s = y_at_rp(9); Yco_surf = y_at_rp(4);
    k_gas_surf = deps.physicalModel.get_gas_properties(T_surf).k_gas; A_p_local = 4*pi*r_p_local^2;
    Q_cond = -k_gas_surf * dTdr_surf * A_p_local;
    m_dot_co_s = -(A_p_local*deps.params.rho_D_gas)*dYco_dr_s + m_dot_local*Yco_surf;
    H_reac_s = deps.thermo.get_reaction_enthalpy('Mg(g)+CO(g)=MgO(s)+C(s)',T_surf);
    Q_reac_s = m_dot_co_s * (-H_reac_s / deps.params.materials.CO.molar_mass);
    Q_rad = deps.params.emissivity*deps.params.sigma*(T_surf^4-deps.params.ambient_temperature^4)*A_p_local;
    Q_evap = m_dot_local*deps.params.materials.Mg.latent_heat;
    res(4) = (Q_cond + Q_reac_s) - (Q_rad + Q_evap);

    % --- 条件 5-9: 远场 r = r_inf (5个) ---
    % 5. 温度等于环境温度
    res(5) = y_at_rinf(1) - deps.params.ambient_temperature;
    % 6. Y_co2 = 1 (纯CO2环境)
    res(6) = y_at_rinf(3) - 1.0;
    % 7. Y_mg = 0 
    res(7) = y_at_rinf(2);
    % 8. Y_co = 0
    res(8) = y_at_rinf(4);
    % 9. Y_mgo = 0
    res(9) = y_at_rinf(5);

    % =========================================================================
    % 组2: 火焰界面条件 r = r_f (共13个)
    % =========================================================================
    
    % --- 条件 10-14: 状态变量连续性 (5个) ---
    res(10) = y_at_rf_left(1) - y_at_rf_right(1); % T
    res(11) = y_at_rf_left(2) - y_at_rf_right(2); % Y_mg
    res(12) = y_at_rf_left(3) - y_at_rf_right(3); % Y_co2
    res(13) = y_at_rf_left(4) - y_at_rf_right(4); % Y_co
    res(14) = y_at_rf_left(5) - y_at_rf_right(5); % Y_mgo
    
    % --- 条件 15-16: 反应物耗尽 (2个) ---
    % 结合连续性条件，这确保了Mg和CO2在整个界面上浓度为0
    res(15) = y_at_rf_right(2); % Mg在火焰面右侧为0
    res(16) = y_at_rf_left(3);  % CO2在火焰面左侧为0
    
    % --- 条件 17-18: 单边导数约束 (2个) ---
    % 描述反应物在火焰面被完全消耗，浓度曲线平滑触底，因此在该侧导数为0
    res(17) = y_at_rf_left(7); % dY_mg/dr @ r_f_left = 0
    res(18) = y_at_rf_right(8);  % dY_co2/dr @ r_f_right = 0
    
    % --- 条件 19: 能量通量跳变 (1个) ---
    mw = deps.params.materials; rho_D_avg_local = deps.params.rho_D_gas; A_f_local = 4*pi*r_f_local^2;
    calc_flux = @(Y, dYdr) -(A_f_local * rho_D_avg_local) * dYdr + m_dot_local * Y;
    J_mg_left  = calc_flux(y_at_rf_left(2), y_at_rf_left(7));
    T_f = y_at_rf_left(1); H_reac_f = deps.thermo.get_reaction_enthalpy('Mg(g)+CO2(g)=MgO(g)+CO(g)',T_f);
    Q_flame = J_mg_left*(-H_reac_f/mw.Mg.molar_mass);
    k_gas_f = deps.physicalModel.get_gas_properties(T_f).k_gas;
    res(19) = (y_at_rf_right(6)-y_at_rf_left(6))*k_gas_f + Q_flame/A_f_local;

    % --- 条件 20-22: 化学计量关系 (3个) ---
    J_co2_right= calc_flux(y_at_rf_right(3), y_at_rf_right(8));
    J_co_jump  = calc_flux(y_at_rf_right(4), y_at_rf_right(9)) - calc_flux(y_at_rf_left(4), y_at_rf_left(9));
    J_mgo_jump = calc_flux(y_at_rf_right(5), y_at_rf_right(10)) - calc_flux(y_at_rf_left(5), y_at_rf_left(10));
    res(20) = J_mg_left/mw.Mg.molar_mass + J_co2_right/mw.CO2.molar_mass;
    res(21) = J_mg_left/mw.Mg.molar_mass - J_co_jump/mw.CO.molar_mass;
    res(22) = J_mg_left/mw.Mg.molar_mass - J_mgo_jump/mw.MgO.molar_mass;
end

function y = initial_guess_func(r, region, pState, params, r_f_guess, r_inf)
    % BVP的初始猜测值 (版本B.1: 修正CO分布)
    % 目标: 提供一个物理上更平滑、形状更真实的初始猜测，以提高收敛性。
    y = zeros(10, 1);
    
    T_p = pState.T_p;
    T_amb = params.ambient_temperature;
    T_f_guess = (T_p + T_amb) / 2 + 800; % 更合理的火焰温度猜测
    Y_mg_s_guess = 0.1; % Mg在表面的质量分数猜测
    Y_co_peak_guess = 0.2; % CO在火焰面的峰值猜测
    
    % 定义火焰过渡带的“厚度”，delta越小，过渡越陡峭
    delta1 = (r_f_guess - pState.r_p) / 4; % 区域1的厚度参数
    delta2 = (r_inf - r_f_guess) / 4; % 区域2的厚度参数
    
    % 使用tanh创建0到1的平滑过渡因子
    % frac1: 从 r_p (frac~0) 平滑过渡到 r_f (frac~1)
    frac1 = (tanh((r - (pState.r_p + r_f_guess)/2) / (delta1+eps)) + 1) / 2;
    % frac2: 从 r_f (frac~0) 平滑过渡到 r_inf (frac~1)
    frac2 = (tanh((r - (r_f_guess + r_inf)/2) / (delta2+eps)) + 1) / 2;
    
    % --- 区域1: r_p -> r_f ---
    if region == 1
        % 温度: 从 T_p 平滑上升到 T_f_guess
        y(1) = T_p + (T_f_guess - T_p) * frac1;
        % Y_mg: 从 Y_mg_s_guess 平滑下降到 0
        y(2) = Y_mg_s_guess * (1 - frac1);
        % Y_co2: 始终为 0
        y(3) = 0;
        % Y_co: 从0 (r_p处) 平滑上升到 Y_co_peak_guess (r_f处)
        y(4) = Y_co_peak_guess * frac1;
        % Y_mgo: 从0平滑上升到某个峰值
        y(5) = 0.1 * frac1;
        
    % --- 区域2: r_f -> r_inf ---
    else % region == 2
        % 温度: 从 T_f_guess 平滑下降到 T_amb
        y(1) = T_f_guess - (T_f_guess - T_amb) * frac2;
        % Y_mg: 始终为0
        y(2) = 0;
        % Y_co2: 从0平滑上升到1
        y(3) = 1.0 * frac2;
        % Y_co: 从 Y_co_peak_guess (r_f处) 平滑下降到 0 (r_inf处)
        y(4) = Y_co_peak_guess * (1 - frac2);
        % Y_mgo: 从某个值平滑下降到0
        y(5) = 0.1 * (1 - frac2);
    end
end

function rate_info_out = package_rate_info_func(m_dot_mg_total, m_dot_co_surf, sol, r_f, params)
    % 结果打包的独立版本
    mw_mg = params.materials.Mg.molar_mass;
    mw_mgo = params.materials.MgO.molar_mass;
    mw_c = params.materials.C.molar_mass;
    mw_co = params.materials.CO.molar_mass;

    rate_info_out.dmdt_mg = -m_dot_mg_total;
    
    m_dot_mgo_surf = m_dot_co_surf * (mw_mgo / mw_co);
    rate_info_out.dmdt_c = m_dot_co_surf * (mw_c / mw_co);
    
    m_dot_mg_surf_equiv = m_dot_co_surf * (mw_mg/mw_co);
    m_dot_mg_flame = m_dot_mg_total - m_dot_mg_surf_equiv;
    m_dot_mgo_flame = m_dot_mg_flame * (mw_mgo / mw_mg);

    rate_info_out.dmdt_mgo = m_dot_mgo_surf + m_dot_mgo_flame;
    
    rate_info_out.r_f = r_f;
    sol_at_rf = deval(sol, r_f);
    rate_info_out.T_f = sol_at_rf(1);
end 