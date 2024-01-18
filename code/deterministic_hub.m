function [f, A, b, Aeq, beq, lb, ub, N, Cost, Emiss, variables] = deterministic_hub(time, Ns, ws, plt, comps, a, R, C_l, C_t, C_g, C_c, irr, obj, f1f2)
    
    S_pannello = 1.768;     % m2/pv panel
    N_pv_pan = 9541;        % num. panels pv
    S_stc_pannello = 2.2;   % m2/panel stc
    N_stc_pan = 7667;       % num. solar thermal panels
    S_roof = 16868;         % m2 max available
    
    dT = 1;                 % time interval of 1 h
    
    emis_pow = 493.8e-6;    % tCO2/kWh
    emis_gas = 1.975e-3;    % tCO2/m^3
    emiss2euro = 30;        % euro/tCO2
    
    c_ess_p = 160;          % euros/kW
    c_ess_e = 240;          % euros/kWh
    c_tes_p = 75;           % euros/kW
    c_tes_e = 125;          % euros/kWh
    c_chp = 800;            % euros/kW
    c_boiler = 200;         % euros/kW
    % boiler a condensazione
    % c_cboiler = 180;      % euros/kW
    c_hp = 900;             % euros/kW
    c_Spv = 226;            % euro/m2
    c_Sstc = 400;           % euro/m2
    c_purch = .2;           % euros/kWh
    c_sold = .05;           % euros/kWh
    c_gas = .8;             % euro/m3
    
    sigma_min_ess = .1;     % min SoC ESS
    sigma_max_ess = .9;     % max SoC ESS
    sigma_min_tes = 0;      % min SoC TES
    sigma_max_tes = .9;     % max SoC TES
    eta = 0.93;             % pv yield ESS and TES
    eta_pv = .23;           % rend. pv
    eta_bos = .91;          % rend. inverter for pv
    eta_stc = .65;          % rend. thermal panels
    eta_e_chp = .4;         % rend. electric chp
    eta_t_chp = .5;         % rend. thermal chp
    eta_boiler = .9;        % rend. boiler
    % condensing boiler
    % eta_cboiler = 1.06;   % boiler yield in condensatezione
    COP = 4;                % performance coeff.
    EER = 4;                % energy efficiency ratio
    PCI_gas = 9.806;        % kWh/m3
    
    % Define variables and variable's indexes
    ro = 0;
    V = var_indexing(time, comps, ro, plt); variables = V{1};
    [Pn_hp, Pn_ess, Pn_tes, Pn_boiler, Pn_chp, S_pv, S_stc, E_ess,...
        E_tes, P_hp, P_ess, P_tes, P_boiler, P_chp, P_pv, P_stc, P_sold,...
        P_purch, SoC_ess, SoC_tes, G_grid, SoC_ess_0, SoC_tes_0, P_ess_0,...
        P_tes_0, Cost, Emiss] = idx_variables(time, comps);

    N = length(variables);
    % CREATE LOWER AND UPPER BOUNDS
    lb = zeros(N,1);
    lb(P_ess) = -Inf;
    lb(P_tes) = -Inf;
    
    ub = Inf(size(variables));
    ub(S_pv) = S_pannello*N_pv_pan;
    ub(S_stc) = S_stc_pannello*N_stc_pan;

    G_dis = 1+1+time*12+4+1; % 4 eq. soc + boiler + 2 eq. chp + hp + 4 eq.
    % P batteries + [4 eq. ratio P/E batt. (no t)] + 1 per function
    % target (costs) in inequality constraints + 1 per function
    % (emissions) in inequality constraints + 1 per constraint
    % max roof space occupancy (S_stc + S_pv = S_roof)
    G_eq = time*8+2+2;  % 4 eq. budget + 2 eq. batteries + 2 eq.
    % allocationz Spv-Ppv and Sstc-Pstc + 2 eq. allocationz SoC0
    A = zeros(G_dis,N);
    b = zeros(G_dis,1);
    Aeq = zeros(G_eq,N);
    beq = zeros(G_eq,1);
    q = 1;
    
    %--- Inequality Constraints ---
    for t=1:G_dis
        
        if t<=time
           % Limit on ESS Minimum SoC
            A(t,SoC_ess(q)) = -1;
            A(t,E_ess) = sigma_min_ess;
            q = q + 1;
            
        elseif t>time && t<=2*time
            % Limit on SoC maximum ESS
            A(t,SoC_ess(q)) = 1;
            A(t,E_ess) = -sigma_max_ess;
            q = q + 1;
            
        elseif t>2*time && t<=3*time
            % Limit on SoC minimal TES.
            A(t,SoC_tes(q)) = -1;
            A(t,E_tes) = sigma_min_tes;
            q = q + 1;
            
        elseif t>3*time && t<=4*time
            % Limit on maximum SoC TES
            A(t,SoC_tes(q)) = 1;
            A(t,E_tes) = -sigma_max_tes;
            q = q + 1;
            
        elseif t>4*time && t<=5*time
            % Limit max boiler
            A(t,P_boiler(q)) = 1;
            A(t,Pn_boiler) = -1;
            q = q + 1;
            
        elseif t>5*time && t<=6*time
            % Limit max hp
            A(t,P_hp(q)) = 1;
            A(t,Pn_hp) = -1;
            q = q + 1;
            
        elseif t>6*time && t<=7*time
            % Limit min chp
            A(t,P_chp(q)) = -1;
            A(t,Pn_chp) = .5;
            q = q + 1;
            
        elseif t>7*time && t<=8*time
            % Limit max chp
            A(t,P_chp(q)) = 1;
            A(t,Pn_chp) = -1;
            q = q + 1;
            
        elseif t>8*time && t<=9*time
            % Limit potenza max ESS
            A(t,P_ess(q)) = -1;
            if (q-1) == 0
                A(t,P_ess_0) = 1;
            else
                A(t,P_ess(q-1)) = 1;
            end
            A(t,Pn_ess) = -1;
            q = q + 1;
            
        elseif t>9*time && t<=10*time
            % Limit potenza min ESS 
            A(t,P_ess(q)) = 1;
            if (q-1) == 0
                A(t,P_ess_0) = -1;
            else
                A(t,P_ess(q-1)) = -1;
            end
            A(t,Pn_ess) = -1;
            q = q + 1;
            
        elseif t>10*time && t<=11*time
            % Limit potenza max TES
            A(t,P_tes(q)) = -1;
            if (q-1) == 0
                A(t,P_tes_0) = 1;
            else
                A(t,P_tes(q-1)) = 1;
            end
            A(t,Pn_tes) = -1;
            q = q + 1;
            
        elseif t>11*time && t<=12*time
            % Limit potenza min TES
            A(t,P_tes(q)) = 1;
            if (q-1) == 0
                A(t,P_tes_0) = 1;
            else
                A(t,P_tes(q-1)) = -1;
            end
            A(t,Pn_tes) = -1;
            q = q + 1;
            
        elseif t>12*time && t<=12*time+1
            % Limit ratio min ESS
            A(t,Pn_ess) = 1;
            A(t,E_ess) = -1;
            
        elseif t>12*time+1 && t<=12*time+2
            % Limit ratio max ESS
            A(t,Pn_ess) = -8; % ratio E/P = 8 (ess)
            A(t,E_ess) = 1;
        
        elseif t>12*time+2 && t<=12*time+3
            % Limit ratio min TES
            A(t,Pn_tes) = 1;
            A(t,E_tes) = -1;
            
        elseif t>12*time+3 && t<=12*time+4
            % Limit ratio max TES
            A(t,Pn_tes) = -5; % ratio E/P = 5 (tes)
            A(t,E_tes) = 1;
            
        elseif t>12*time+4 && t<=12*time+5
            A(t,Cost) = -1;
            A(t,Emiss) = emiss2euro;
            % Investment costs
            A(t,Pn_hp) = c_hp * a(1) * (1-R(1));
            A(t,Pn_ess) = c_ess_p * a(2) * (1-R(2));
            A(t,E_ess) = c_ess_e * a(2) * (1-R(2));
            A(t,Pn_tes) = c_tes_p * a(3) * (1-R(3));
            A(t,E_tes) = c_tes_e * a(3) * (1-R(3));
            A(t,Pn_boiler) = c_boiler * a(4) * (1-R(4));
            A(t,Pn_chp) = c_chp * a(5) * (1-R(5));
            A(t,S_pv) = c_Spv * a(6) * (1-R(6));
            A(t,S_stc) = c_Sstc * a(7) * (1-R(7));

            % Operational costs
            for h=1:Ns
                A(t,P_purch) = c_purch * dT * 365 * ws(h);
                A(t,P_sold) = c_sold * dT * 365 * ws(h);
                A(t,G_grid) = c_gas * dT * 365 * ws(h);
            end
            
        elseif t>12*time+5 && t<=12*time+6
            A(t,Emiss) = -1;
            % Hourly emissions
            A(t,P_purch) = emis_pow;
            A(t,G_grid) = emis_gas;
            
        else
            % Max roof area limit
            A(t,S_pv) = 1;
            A(t,S_stc) = 1;
            b(t) = S_roof;

        end
        
        if q > time
            q = 1;
        end
    end
    
    
    %--- Equality Constraints ---
    q = 1;
    for t=1:G_eq
        if t<=time 
            % Equations SoC ESS
            Aeq(t,SoC_ess(q)) = -1/(eta*dT);
            if (q-1) == 0
                Aeq(t,SoC_ess_0) = 1/(eta*dT);
            else
                Aeq(t,SoC_ess(q-1)) = 1/(eta*dT);
            end
            Aeq(t,Pn_ess) = 1;
            q = q + 1;
            
        elseif t>time && t<=2*time
            % Equations SoC TES
            Aeq(t,SoC_tes(q)) = -1/(eta*dT);
            if (q-1) == 0
                Aeq(t,SoC_tes_0) = 1/(eta*dT);
            else
                Aeq(t,SoC_tes(q-1)) = 1/(eta*dT);
            end
            Aeq(t,Pn_tes) = 1;
            q = q + 1;
            
        elseif t>2*time && t<=3*time
            % Equations photovoltaic PV panels
            Aeq(t,P_pv(q)) = 1;
            Aeq(t,S_pv) = -eta_pv*eta_bos*irr(q);
            q = q + 1;
            
        elseif t>3*time && t<=4*time
            % Equations photovoltaic panels STC
            Aeq(t,P_stc(q)) = 1;
            Aeq(t,S_stc) = -eta_stc*irr(q);
            q = q + 1;
            
        elseif t>4*time && t<=5*time
            % ELECTRIC balance equations
            Aeq(t,P_sold(q)) = -1;
            Aeq(t,P_purch(q)) = 1;
            Aeq(t,P_chp(q)) = 1;
            Aeq(t,P_pv(q)) = 1;
            Aeq(t,P_hp(q)) = -1;
            Aeq(t,P_ess(q)) = -1;
            
            beq(t) = C_l(q);
            q = q + 1;
            
        elseif t>5*time && t<=6*time
            % THERMAL balance equation
            Aeq(t,P_chp(q)) = eta_t_chp/eta_e_chp;
            Aeq(t,P_stc(q)) = 1;
            Aeq(t,P_hp(q)) = COP;
            Aeq(t,P_tes(q)) = -1;
            Aeq(t,P_boiler(q)) = eta_boiler;
            
            beq(t) = C_t(q);
            q = q + 1;
            
        elseif t>6*time && t<=7*time
            % GAS budget equation
            Aeq(t,P_chp(q)) = -eta_t_chp/(PCI_gas*eta_e_chp);
            Aeq(t,P_boiler(q)) = -1/PCI_gas;
            Aeq(t,G_grid(q)) = 1;
            
            beq(t) = C_g(q);
            q = q + 1;
            
        elseif t>7*time && t<=8*time
            % COOLING budget equation
            Aeq(t,P_hp(q)) = EER;
            
            beq(t) = C_c(q);
            q = q + 1;
            
        elseif t>8*time && t<=8*time+1
            Aeq(t,SoC_ess_0) = 1;
            Aeq(t,E_ess) = -.05;
            
        elseif t>8*time+1 && t<=8*time+2
            Aeq(t,SoC_tes_0) = 1;
            Aeq(t,E_tes) = -.05;
            
        elseif t>8*time+2 && t<=8*time+3
            % Initial and final ESS balance equation
            Aeq(t,SoC_ess_0) = 1;
            Aeq(t,SoC_ess(end)) = -1;
            
        elseif t>8*time+3 && t<=8*time+4
            % Initial and final TES budget equation
            Aeq(t,SoC_tes_0) = 1;
            Aeq(t,SoC_tes(end)) = -1;
            
        end
        
        if q > time
            q = 1;
        end
    end
    
    %---- Objective function ----
    % f = @(x) [ x(Cost), x(Emiss)];
    f = zeros(N,1);
    f(Cost) = 1;
    
    if plt
        variables_sym = V{2}; clear V;
        
        disp('Objective function')
        f' * variables_sym

        disp('Inequality Constraints')
        A * variables_sym

        disp('Equality Constraints')
        Aeq * variables_sym
    end
    
end

