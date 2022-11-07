function Trp_RES_CO2_pre(obj, varargin)
    % preparation for the crossover RES model of CO2
    % format:
    %  - Trp_RES_CO2_pre(point)
    %  - Trp_RES_CO2_pre(point, prop_req)

    point = varargin{1};
    switch nargin
    case 2
        prop_req = 'BOTH';
    case 3
        prop_req = varargin{2};
    end
    
    kexi = cr_kexi(obj, point);
    if point.mc_num == 2
        switch prop_req
        case {'TCX', 'BOTH'}
            cr_chempot_T = cr_ChemPotD21(obj, point, 'T_P');
            cr_chempot_n = cr_ChemPotD21(obj, point, 'n_P');
            point.tcx.chempot_T = cr_chempot_T;
            point.tcx.chempot_n = cr_chempot_n;
            [alpha_bg, beta_bg, d_gamma_bg] = nec_coeff_bg(obj, point);
            point.tcx.alpha_bg = alpha_bg;
            point.tcx.beta_bg = beta_bg;
            point.tcx.d_gamma_bg = d_gamma_bg;
        end
    end
    
    switch prop_req
    case 'VIS'
        point.vis.cr_kexi = kexi;
    case 'TCX'
        point.tcx.cr_kexi = kexi;
    case 'BOTH'
        point.vis.cr_kexi = kexi;
        point.tcx.cr_kexi = kexi;
    end
end

function val = cr_kexi(obj, point)
    % kexi, correlation length, [nm]
    % ref: Int J Thremophys (2013) 34: 191-212
    
    eps = 1E-10;

    % eq.(36)
    nu = 0.630;
    gamma = 1.239;

    %Tcr = obj.eos.Tcr_mix(point.X);
    Dcr = point.Dcr; %1 / obj.eos.Vcr_mix(point.X);
    %Pcr = obj.eos.P('0', Tcr, Dcr, point.X);

    if ismethod(obj.eos, 'checkCriticalExponent')
        Gamma0 = cr_Gamma0(obj, point);
        kexi0 = cr_kexi0(obj, point); 
        Tref = cr_Tref(obj, point);
    else
        Gamma0 = point.X * obj.tcx.Gamma0; % debug here.
        kexi0 = point.X * obj.tcx.kexi0;
        Tref = point.Tcr * obj.tcx.Trefr; % [K]
    end
    
    if point.mc_num == 1
    %    XT = point.D / obj.eos.P(point, 'D');
    %    XTref = point.D / obj.eos.P('D', Tref, point.D, point.X);
        
        
        XT = point.D / obj.tcx.eos.P(point, 'D');
        XTref = point.D / obj.tcx.eos.P('D', Tref, point.D, point.X);
        
        
    else
        XT = point.D / cr_P_D_constChemPot(obj, point);
        XTref = point.D / cr_P_D_constChemPot(obj, point, Tref);
    end
    DXT = XT - XTref * Tref / point.T;
    if DXT < eps
        DXT = eps;
    end
    val = kexi0 * (point.Pcr / Dcr^2 * DXT / Gamma0)^(nu / gamma); % [nm] original
    % val = kexi0 * (Phys.Rg * Tcr / Dcr * DXT / Gamma0)^(nu / gamma); % [nm] proposed by Sengers, Chap 10 of Advances in ...

end

function val = cr_qC1(obj)
    point_cr = Point('MTRCODE',obj,'MC',{obj.name},'PHASE','FL','T',obj.fluid.Tcr,'D',obj.fluid.Dcr);
    obj.calProp({'CP','Cratio','SRES','VIS0','VIS','TCX0','TCX'}, point_cr);
    Gamma0 = cr_Gamma0(obj, point_cr);
    kexi0 = cr_kexi0(obj, point_cr); 
    P_T = obj.eos.P('T',obj.fluid.Tcr,obj.fluid.Dcr,1);
    % carbon dioxide: P_T = 1.593395912963463e+05
    qC = Phys.kB * obj.fluid.Tcr^2 / 16 / point_cr.vis.visbg / point_cr.tcx.tcxbg / obj.fluid.Pcr * Gamma0*1e9 / kexi0^2 * P_T^2;
    val = 1 / qC;
    % carbon dioxide: qC1 = 3.131177361390885
end

function val = cr_P_D_constChemPot_USELESS(obj, point, T)
    switch nargin
    case 2
        pt = point;
    case 3
        pt = obj.copyCleanPoint(point, T, [], -1, []);
    end

    chempotd21_P_D = cr_ChemPotD21(obj, pt, 'P_D');
    chempotd21_D_P = cr_ChemPotD21(obj, pt, 'D_P');

    % P_n = obj.eos.P(pt, 'n');
    val = - chempotd21_P_D / chempotd21_D_P;
end

function val = cr_P_D_constChemPot(obj, point, T)
    switch nargin
    case 2
        pt = point;
    case 3
        pt = obj.copyCleanPoint(point, T, [], -1, []);
    end

    chempotd21_n_D = cr_ChemPotD21(obj, pt, 'n_D');
    chempotd21_D_n = cr_ChemPotD21(obj, pt, 'D_n');

    P_n = obj.eos.P(pt, 'n');
    val = obj.eos.P(pt, 'D') - P_n(2) * chempotd21_D_n / chempotd21_n_D;
end

function val = cr_ChemPotD21(obj, varargin)
    % To calculate the difference of chemical potentials of component 2 and component 1 [J/mol] and its partial derivatives.
    % format:
    %  - cr_ChemPotD21(point)
    %  - cr_ChemPotD21(point, pd)

    point = varargin{1};
    if point.mc_num ~= 2
        errordlg('Material.cr_ChemPotD21: Only 2-component system is supported.');
    end

    switch nargin
    case 2
        pd = '0';
    case 3
        pd = varargin{2};
    end

    switch pd
    case '0'
        chempot = obj.ChemPot(point);
        val = chempot(2) - chempot(1);
        return
    case 'T_P'
        step = point.T * 1e-8;
        point_p = obj.copyCleanPoint(point, point.T + step, -1, [], []);
        point_m = obj.copyCleanPoint(point, point.T - step, -1, [], []);
    case 'P_D'
        step = point.P * 1e-4;
        point_p = obj.copyCleanPoint(point, [], [], point.P + step, -1);
        point_m = obj.copyCleanPoint(point, [], [], point.P - step, -1);
    case 'D_P'
        step = point.D * 1e-4;
        point_p = obj.copyCleanPoint(point, [], point.D + step, [], -1);
        point_m = obj.copyCleanPoint(point, [], point.D - step, [], -1);
    case 'D_n'
        step = point.D * 1e-6;
        point_p = obj.copyCleanPoint(point, [], point.D + step, -1, []);
        point_m = obj.copyCleanPoint(point, [], point.D - step, -1, []);
    case 'n_P'
        [x_p, x_m, n_p, n_m] = Bsc.PartialMolar(point.X, 2);
        step = n_p - n_m;
        point_p = obj.copyCleanPoint(point, [], -1, [], x_p);
        point_m = obj.copyCleanPoint(point, [], -1, [], x_m);
    case 'n_D'
        [x_p, x_m, n_p, n_m] = Bsc.PartialMolar(point.X, 2);
        step = n_p - n_m;
        point_p = obj.copyCleanPoint(point, [], [], -1, x_p);
        point_m = obj.copyCleanPoint(point, [], [], -1, x_m);
    end

    chempot_p = cr_ChemPotD21(obj, point_p);
    chempot_m = cr_ChemPotD21(obj, point_m);
    %d_chempot = (chempot_p - chempot_m) / step;
    %val = d_chempot(2) - d_chempot(1);
    val = (chempot_p - chempot_m) / step;
end

function val = cr_Gamma0(obj, point)
    gamma = 1.239;

    Tcr = point.Tcr; %obj.eos.Tcr_mix(point.X);
    Dcr = point.Dcr; %1 / obj.eos.Vcr_mix(point.X);
    %Pcr = obj.eos.P('0', Tcr, Dcr, point.X);

    kapa = 1e-12;
    Tcrp = Tcr * (1 + kapa);
    if point.mc_num == 1
        XTr = point.Pcr / Dcr / obj.eos.P('D', Tcrp, Dcr, point.X); % original
        % XTr = Phys.Rg * Tcr / obj.eos.P('D', Tcrp, Dcr, point.X); % Sengers
    else
        point_Dcr = obj.copyCleanPoint(point, Tcrp, Dcr, -1, []);
        XTr = point.Pcr / Dcr / cr_P_D_constChemPot(obj, point_Dcr); % original
        % XTr = Phys.Rg * Tcr / cr_P_D_constChemPot(obj, point_Dcr); % Sengers
    end
    val = XTr * kapa^gamma;

    if val < 0
        errordlg('Trp_RES_CO2_pre.cr_Gamma0: Gamma0 < 0.');
    end
end

function val = cr_kexi0(obj, point)
    %kexi0 = point.X * obj.vis.kexi0;
    %kexi0 = point.X * obj.tcx.kexi0;
    
    alpha = 0.110;

    Tcr = point.Tcr;
    Dcr = point.Dcr;
    tt = 1e-12;
    T = Tcr / (1 - tt);
    pt = obj.copyCleanPoint(point, T, Dcr, -1, []);
    CV = obj.CV(pt);
    
    val = ( 0.0188 * Phys.kB / alpha / Dcr / CV * tt^(-alpha) )^(1/3) * 1e9; % [nm]
    % val = point.X * obj.vis.kexi0;

    if val < 0
        errordlg('Trp_RES_CO2_pre.cr_kexi0: kexi0 < 0.');
    end
end

function val = cr_Tref(obj, point)
    %val = Tcr * obj.vis.Trefr; % [K]
    %val = Tcr * obj.tcx.Trefr; % [K]
    Tcr = point.Tcr;
    val = max(min([obj.fluid.Tupper - 50; 2*Tcr]), 1.5*Tcr);
end

function [alpha_bg, beta_bg, d_gamma_bg] = nec_coeff_bg(obj, point)
    % ref: Kiselev. Fluid Phase Equilibria 142 (1998) 253-280.

    % eq.(42): atomic diffusion volumn
    sumv = zeros(point.mc_num, 1);
    for i_mc = 1 : point.mc_num
        switch point.mc_name{i_mc}
        case 'CO2'
            sumv(i_mc) = 26.90; 
        case 'ethane'
            sumv(i_mc) = 44.88; 
        end
    end

    % Table 5
    coeff_alpha = zeros(19, 1);
    coeff_beta  = zeros(19, 1);
    coeff_gamma = zeros(12, 1);
    coeff_alpha(6) = 1.33315e-6;
    coeff_alpha(7) = -2.34664e-6;
    coeff_alpha(12) = -1.43428e-6;
    coeff_alpha(13) = 1.70776e-6;
    coeff_alpha(18) = 1.96473e-7;
    coeff_alpha(19) = -2.44598e-7;
    coeff_beta(1) = 2.40238;
    coeff_beta(3) = -7.43324e-6;
    coeff_beta(6) = 6.39961e-6;
    coeff_beta(12) = -1.15818e-6;
    coeff_beta(18) = 9.40929e-8;
    coeff_gamma(1) = -1.17562e-6;
    coeff_gamma(2) = -5.14420e-7;
    coeff_gamma(3) = 5.54845e-6;
    coeff_gamma(5) = -6.72675e-6;
    coeff_gamma(6) = 2.60059e-6;
    coeff_gamma(9) = 1.74581e-6;
    coeff_gamma(10) = -1.24763e-6;

    x = point.X(2);

    % eq.(39): binary diffusion coefficient [m2/s]
    D0 = 1.01325e-2 * point.T^1.75 * (1e-3 * (obj.fluid.MW0(1) + obj.fluid.MW0(2))/(obj.fluid.MW0(1) * obj.fluid.MW0(2)))^(1/2) / point.P / (sum(sumv.^(1/3)))^2;

    % eq.(36~37): dilute gas kinetic coefficient
    alpha_0 = (point.D*1e-3) * D0 * (point.MW*1e3)^2 / (Phys.Rg*1e3) / point.T * prod(point.X); % [kg s /m3]
    beta_0 = (Phys.Rg*1e3) * alpha_0 * (coeff_beta(1) + x * coeff_beta(2) - log(x / (1 - x))); % [kg /m /s /K]
    
    % eq.(33): composition-dependent coefficient
    Gamma_cx = point.Tcr^(1/6) * (point.Pcr/point.Dcr/Phys.Rg/point.Tcr)^5 * sqrt(point.MW*1e3) / (point.Pcr*1e-6)^(2/3);

    Dr = point.Dr;

    % eq.(31~32, 34)
    k_alpha = zeros(6, 1);
    k_beta = zeros(6, 1);
    k_gamma = zeros(6, 1);
    for kk = 1 : 6
        k_alpha(kk) = (coeff_alpha(3*kk) + coeff_alpha(3*kk+1) * x) * Dr^(kk+1);
        k_beta(kk) = (coeff_beta(3*kk) + coeff_beta(3*kk+1) * x) * Dr^(kk+1);
        k_gamma(kk) = (coeff_gamma(2*kk-1) + coeff_gamma(2*kk) * x) * Dr^kk;
    end

    alpha_ex =prod(point.X) * (1000*Phys.Rg)^(-7/6) / point.Tcr / Gamma_cx * sum(k_alpha);
    beta_ex = prod(point.X) * (1000*Phys.Rg)^(-1/6) / point.Tcr / Gamma_cx * sum(k_beta);

    % eq.(29~30)
    alpha_bg = alpha_0 + alpha_ex;
    beta_bg = beta_0 + beta_ex;

    % eq.(34)
    d_gamma_bg = point.T * beta_bg^2 / alpha_bg + prod(point.X) * (1000*Phys.Rg)^(5/6) / Gamma_cx * sum(k_gamma);
end