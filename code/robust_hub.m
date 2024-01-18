function [f, A, b, Aeq, beq, lb, ub, N, Cost, Emiss, variables] = robust_hub(time, Ns, ws, plt, comps, a, R, C_l, C_t, C_g, C_c, irr, symmetric)
    
    S_pannello = 1.768;     % m2/pv panel
    N_pv_pan = 9541;        % num. panels pv
    S_stc_pannello = 2.2;   % m2/panel stc
    N_stc_pan = 7667;       % num. solar thermal panels
    S_roof = 16868;         % m2 max available
    
    dT = 1;                 % 1 h time interval
    
    emis_pow = 493.8e-6;    % tCO2/kWh
    emis_gas = 1.975e-3;    % tCO2/m^3
    emiss2euro = 30;        % euro/tCO2
    
    c_ess_p = 160;          % â'¬/kW
    c_ess_e = 240;          % â'¬/kWh
    c_tes_p = 75;           % â'¬/kW
    c_tes_e = 125;          % â'¬/kWh
    c_chp = 800;            % â'¬/kW
    c_boiler = 200;         % â'¬/kW
    % boiler a condensazione
    % c_cboiler = 180;      % â'¬/kW
    c_hp = 900;             % â'¬/kW
    c_Spv = 226;            % â'¬/m2
    c_Sstc = 400;           % â'¬/m2
    c_purch = .2;           % â'¬/kWh
    c_sold = .05;           % â'¬/kWh
    c_gas = .8;             % â'¬/m3
    
    sigma_min_ess = .1;     % min SoC ESS
    sigma_max_ess = .9;     % max SoC ESS
    sigma_min_tes = 0;      % min SoC TES
    sigma_max_tes = .9;     % max SoC TES
    eta = 0.93;             % batt. yield ESS and TES
    eta_pv = .23;           % rend. pv
    eta_bos = .91;          % rend. inverter for pv
    eta_stc = .65;          % rend. thermal panels
    eta_e_chp = .4;         % rend. electric chp
    eta_t_chp = .5;         % rend. thermal chp
    eta_boiler = .9;        % rend. boiler
    % condensing boiler
    % eta_cboiler = 1.06;   % rend. condensing boiler
    COP = 4;                % performance coeff.
    EER = 4;                % energy efficiency ratio
    PCI_gas = 9.806;        % kWh/m3
    
    
    epsilon = .1;           % Uncertainty of variables (or parameters).
    
    ro = 1;
    V = var_indexing(time, comps, ro, plt, symmetric); variables = V{1};
    [Pn_hp, Pn_ess, Pn_tes, Pn_boiler, Pn_chp, S_pv, S_stc, E_ess,...
        E_tes, P_hp, P_ess, P_tes, P_boiler, P_chp, P_pv, P_stc, P_sold,...
        P_purch, SoC_ess, SoC_tes, G_grid, SoC_ess_0, SoC_tes_0, P_ess_0,...
        P_tes_0, Cost, Emiss, z1, z2, z3, z4, z5, z6, z7, ...
        p1, p2, p3, p4, p5, p6, p7] = idx_variables(time, comps);
        
    % CREATE VARIABLES FOR INDEXING
    N = length(variables); 
    for v = 1:N 
        eval([variables{v},' = ', num2str(v),';'])
    end

    
    % CREATE LOWER AND UPPER BOUNDS
    lb = zeros(N,1);
    lb(P_ess) = -Inf;
    lb(P_tes) = -Inf;
    
    ub = Inf(size(variables));
    ub(S_pv) = S_pannello*N_pv_pan;
    ub(S_stc) = S_stc_pannello*N_stc_pan;
    
    
    % #### Costruiamo i parametri per il modello robusto ####
    % NOTA = length(comps) = 7 (+2 per energia ess e tes); numero risorse
    % considerate nello studio di pianificazione energy hub + 1 per il
    % costo delle emissioni
    
    G_dis_ro = 3*time+(length(comps)+2)+1+2*4*time+2*2*time; % time*3 
    % poiché 3 variabili incognite nella prima equazione variabili nel tempo
    % + 10 variabili incognite non variabili + 4*time per le quazioni di
    % bilancio ma siccome sono eq. equality - disequality sono 
    % moltiplicate per 2 + 2*time per le equazioni dell'irraggiamento,
    % anche queste moltiplicate per 2 poiché eq. di equality
    % originariamente
    
    % Numero di equazioni (senza considerare le variazioni nel tempo, 
    % tanto il numero di variabiliincerte rimane sempre uguale al
    % variare del tempo) in cui compaiono le variabili incerte. 
    % Queste compaiono nelle equazioni di equality 
    % (24 * 4, 4 poichè sono 4 i profili di domanda) e in più 1 nella 
    % funzione obiettivo e 2 per l'irraggiamento
    L_eqRobuste = 7;
    Gamma = zeros(L_eqRobuste,1);
    % num. variabili incerte per primo vincolo (objective fcn)
    % Il numero sarà  pari al numero di componenti + 3 per considerare i
    % costi di vendita e acquisto energia elettrica e acquisto gas + 1 per
    % considerare il fattore di trasformazione da tCO2 ad euro delle
    % emissioni
    J = (length(comps)+2) + 3 + 1;
    for l=1:L_eqRobuste
        if l > 1
            J = 1;
        end
        Gamma(l,1) = gamma_def(J,epsilon,0,0); 
        % Il 3° elemento indica il metodo per il calcolo di Gamma 
        % (0 = epsilon ; 1 = combinatoriale sempl. ;
        % 2 = metodo semplif. combinatoriale), e lo 0 al quarto indica che 
        % non si vuole il plot del Gamma al variare di J
    end
    
    
    G_dis = 1+1+1+time*12+4+time*8+time*4+G_dis_ro; % 4 eq. soc + boiler + 2 eq. chp
    % + hp + 4 eq. P batterie + [4 eq. ratio P/E batt. (no t)] 
    % + 1 per funzione obiettivo nei vincoli di disuguaglianza + 1 per funzione
    % (emissioni) nei vincoli di disuguaglianza
    % + 1 per vincolo max occupazione spazio dei tetti
    % (S_stc + S_pv = S_roof) + 4 eq. bilancio 
    % (essendo scritte come less equal and upper equal sono 4*2=8) + 2 eq.
    % dell'irraggiamente (essendo scritte come less equal and upper equal
    % sono 2*2=4)
    G_eq = time*2+2+2;  % 2 eq. batterie (var. temp) + 2 eq. assegnaz Spv-Ppv 
    % (var. temp) e Sstc-Pstc + 2 eq. assegnaz SoC0 (no t) + 2 eq. per
    % SoC_0 == SoC_T
    
    A = zeros(G_dis,N);
    b = zeros(G_dis,1);
    Aeq = zeros(G_eq,N);
    beq = zeros(G_eq,1);
    q = 1;    
    
    %---- Inequality Constraints ----
    for t=1:G_dis
        
        if t<=time
            % Limite sull'SoC minimo ESS
            A(t,SoC_ess(q)) = -1;
            A(t,E_ess) = sigma_min_ess;
            q = q + 1;
            
        elseif t>time && t<=2*time
            % Limite sull'SoC massimo ESS
            A(t,SoC_ess(q)) = 1;
            A(t,E_ess) = -sigma_max_ess;
            q = q + 1;
            
        elseif t>2*time && t<=3*time
            % Limite sull'SoC minimio TES
            A(t,SoC_tes(q)) = -1;
            A(t,E_tes) = sigma_min_tes;
            q = q + 1;
            
        elseif t>3*time && t<=4*time
            % Limite sull'SoC massimo TES
            A(t,SoC_tes(q)) = 1;
            A(t,E_tes) = -sigma_max_tes;
            q = q + 1;
            
        elseif t>4*time && t<=5*time
            % Limite max boiler
            A(t,P_boiler(q)) = 1;
            A(t,Pn_boiler) = -1;
            q = q + 1;
            
        elseif t>5*time && t<=6*time
            % Limite max hp
            A(t,P_hp(q)) = 1;
            A(t,Pn_hp) = -1;
            q = q + 1;
            
        elseif t>6*time && t<=7*time
            % Limite min chp
            A(t,P_chp(q)) = -1;
            A(t,Pn_chp) = .5;
            q = q + 1;
            
        elseif t>7*time && t<=8*time
            % Limite max chp
            A(t,P_chp(q)) = 1;
            A(t,Pn_chp) = -1;
            q = q + 1;
            
        elseif t>8*time && t<=9*time
            % Limite potenza max ESS
            A(t,P_ess(q)) = -1;
            if (q-1) == 0
                A(t,P_ess_0) = 1;
            else
                A(t,P_ess(q-1)) = 1;
            end
            A(t,Pn_ess) = -1;
            q = q + 1;
            
        elseif t>9*time && t<=10*time
            % Limite potenza min ESS 
            A(t,P_ess(q)) = 1;
            if (q-1) == 0
                A(t,P_ess_0) = -1;
            else
                A(t,P_ess(q-1)) = -1;
            end
            A(t,Pn_ess) = -1;
            q = q + 1;
            
        elseif t>10*time && t<=11*time
            % Limite potenza max TES
            A(t,P_tes(q)) = -1;
            if (q-1) == 0
                A(t,P_tes_0) = 1;
            else
                A(t,P_tes(q-1)) = 1;
            end
            A(t,Pn_tes) = -1;
            q = q + 1;
            
        elseif t>11*time && t<=12*time
            % Limite potenza min TES
            A(t,P_tes(q)) = 1;
            if (q-1) == 0
                A(t,P_tes_0) = 1;
            else
                A(t,P_tes(q-1)) = -1;
            end
            A(t,Pn_tes) = -1;
            q = q + 1;
            
        elseif t>12*time && t<=12*time+1
            % Limite ratio min ESS
            A(t,Pn_ess) = 1;
            A(t,E_ess) = -1;
            
        elseif t>12*time+1 && t<=12*time+2
            % Limite ratio max ESS
            A(t,Pn_ess) = -8; % ratio E/P = 8 (ess)
            A(t,E_ess) = 1;
        
        elseif t>12*time+2 && t<=12*time+3
            % Limite ratio min TES
            A(t,Pn_tes) = 1;
            A(t,E_tes) = -1;
            
        elseif t>12*time+3 && t<=12*time+4
            % Limite ratio max TES
            A(t,Pn_tes) = -5; % ratio E/P = 5 (tes)
            A(t,E_tes) = 1;
            
        elseif t>12*time+4 && t<=12*time+5
            % PRIMA FUNZIONE OBIETTIVO: Costi
            A(t,Cost) = -1;
            A(t,Emiss) = emiss2euro;
            % Costo di investimento
            A(t,Pn_hp) = c_hp * a(1) * (1-R(1));
            A(t,Pn_ess) = c_ess_p * a(2) * (1-R(2));
            A(t,E_ess) = c_ess_e * a(2) * (1-R(2));
            A(t,Pn_tes) = c_tes_p * a(3) * (1-R(3));
            A(t,E_tes) = c_tes_e * a(3) * (1-R(3));
            A(t,Pn_boiler) = c_boiler * a(4) * (1-R(4));
            A(t,Pn_chp) = c_chp * a(5) * (1-R(5));
            A(t,S_pv) = c_Spv * a(6) * (1-R(6));
            A(t,S_stc) = c_Sstc * a(7) * (1-R(7));
            % Variabile polyhedral set
            A(t,z1) = Gamma(1,1);
            % Variabili box set
            A(t,p1) = 1;
            
            % Costi di operation
            for h=1:Ns
                A(t,P_purch) = c_purch * dT * 365 * ws(h);
                A(t,P_sold) = c_sold * dT * 365 * ws(h);
                A(t,G_grid) = c_gas * dT * 365 * ws(h);
            end
            
        elseif t>12*time+5 && t<=12*time+5+time
            % Vincolo addizionale Costi di acquisto
            A(t,P_purch(q)) = epsilon*c_purch;
            A(t,z1) = -1;
            A(t,p1(q)) = -1; % A(t,p1(1:time)) = -1;
            q = q + 1;
        
        elseif t>12*time+5+time && t<=12*time+5+2*time
            % Vincolo addizionale Costi di vendita
            A(t,P_sold(q)) = epsilon*c_sold;
            A(t,z1) = -1;
            A(t,p1(time+q)) = -1; % A(t,p1(time+1:2*time)) = -1;
            q = q + 1;
            
        elseif t>12*time+5+2*time && t<=12*time+5+3*time
            % Vincolo addizionale Costi di acquisto gas dalla rete
            A(t,G_grid(q)) = epsilon*c_gas;
            A(t,z1) = -1;
            A(t,p1(2*time+q)) = -1; % A(t,p1(2*time+1:3*time)) = -1;
            q = q + 1;
            
        elseif t>12*time+5+3*time && t<=12*time+5+3*time+1
            % Vincolo addizionale Costi di investimento Potenza ESS
            A(t,Pn_ess) = epsilon*c_ess_p;
            A(t,z1) = -1;
            A(t,p1(3*time+1)) = -1;
            
        elseif t>12*time+5+3*time+1 && t<=12*time+5+3*time+2
            % Vincolo addizionale Costi di investimento Energia ESS
            A(t,E_ess) = epsilon*c_ess_e;
            A(t,z1) = -1;
            A(t,p1(3*time+2)) = -1;
            
        elseif t>12*time+5+3*time+2 && t<=12*time+5+3*time+3
            % Vincolo addizionale Costi di investimento Potenza TES
            A(t,Pn_tes) = epsilon*c_tes_p;
            A(t,z1) = -1;
            A(t,p1(3*time+3)) = -1;
            
        elseif t>12*time+5+3*time+3 && t<=12*time+5+3*time+4
            % Vincolo addizionale Costi di investimento Energia TES
            A(t,E_tes) = epsilon*c_tes_e;
            A(t,z1) = -1;
            A(t,p1(3*time+4)) = -1;
        
        elseif t>12*time+5+3*time+4 && t<=12*time+5+3*time+5
            % Vincolo addizionale Costi di investimento CHP
            A(t,Pn_chp) = epsilon*c_chp;
            A(t,z1) = -1;
            A(t,p1(3*time+5)) = -1;
            
        elseif t>12*time+5+3*time+5 && t<=12*time+5+3*time+6
            % Vincolo addizionale Costi di investimento HP
            A(t,Pn_hp) = epsilon*c_hp;
            A(t,z1) = -1;
            A(t,p1(3*time+6)) = -1;
            
        elseif t>12*time+5+3*time+6 && t<=12*time+5+3*time+7
            % Vincolo addizionale Costi di investimento BOILER
            A(t,Pn_boiler) = epsilon*c_boiler;
            A(t,z1) = -1;
            A(t,p1(3*time+7)) = -1;
            
        elseif t>12*time+5+3*time+7 && t<=12*time+5+3*time+8
            % Vincolo addizionale Costi di investimento PV
            A(t,S_pv) = epsilon*c_Spv;
            A(t,z1) = -1;
            A(t,p1(3*time+8)) = -1;
            
        elseif t>12*time+5+3*time+8 && t<=12*time+5+3*time+9
            % Vincolo addizionale Costi di investimento STC
            A(t,S_stc) = epsilon*c_Sstc;
            A(t,z1) = -1;
            A(t,p1(3*time+9)) = -1;
        
        elseif t>12*time+5+3*time+9 && t<=12*time+5+3*time+10
            % Vincolo addizionale Costi di Emissione
            A(t,Emiss) = epsilon*emiss2euro;
            A(t,z1) = -1;
            A(t,p1(3*time+10)) = -1;
            
        elseif t>12*time+5+3*time+10 && t<=12*time+5+3*time+11
            % SECONDA FUNZIONE OBIETTIVO: Emissioni
            A(t,Emiss) = -1;
            
            % Emissioni orarie
            A(t,P_purch) = emis_pow;
            A(t,G_grid) = emis_gas;
            
        elseif t>12*time+5+3*time+11 && t<=12*time+6+3*time+11
            A(t,S_pv) = 1;
            A(t,S_stc) = 1;
            b(t) = S_roof;
            
        elseif t>12*time+6+3*time+11 && t<=12*time+6+time+3*time+11
            % Equazione di bilancio ELETTRICO (segno disequazione: <)
            A(t,P_sold(q)) = -1;
            A(t,P_purch(q)) = 1;
            A(t,P_chp(q)) = 1;
            A(t,P_pv(q)) = 1;
            A(t,P_hp(q)) = -1;
            A(t,P_ess(q)) = 1;
            A(t,p2(q)) = 1;
            A(t,z2(q)) = min(Gamma(2,1),1);
            
            
            b(t) = C_l(q);
            q = q + 1;
            
        elseif t>15*time+6+time+11 && t<=15*time+6+time+11+time
            % Vincolo addizionale Equazione di Bilancio ELETTRICO (<)
            A(t,p2(q)) = -1;
            A(t,z2(q)) = -1;
            b(t) = -epsilon*C_l(q);
            q = q + 1;
            
        elseif t>15*time+6+time+11+time && t<=15*time+6+2*time+11+time
            % Equazione di bilancio ELETTRICO (segno disequazione: >)
            A(t,P_sold(q)) = 1;
            A(t,P_purch(q)) = -1;
            A(t,P_chp(q)) = -1;
            A(t,P_pv(q)) = -1;
            A(t,P_hp(q)) = 1;
            A(t,P_ess(q)) = -1;
            A(t,p2(q)) = -1;
            A(t,z2(q)) = -min(Gamma(2,1),1);
            
            b(t) = -C_l(q);
            q = q + 1;
        
        elseif t>15*time+6+2*time+11+time && t<=15*time+6+2*time+11+2*time
            % Vincolo addizionale Equazione di Bilancio ELETTRICO (>)
            A(t,p2(q)) = 1;
            A(t,z2(q)) = 1;
            b(t) = epsilon*C_l(q);
            q = q + 1;
            
        elseif t>15*time+6+2*time+11+2*time && t<=15*time+6+3*time+11+2*time
            % Equazione di bilancio TERMICO (segno disequazione: <)
            A(t,P_chp(q)) = eta_t_chp/eta_e_chp;
            A(t,P_stc(q)) = 1;
            A(t,P_hp(q)) = COP;
            A(t,P_tes(q)) = 1;
            A(t,P_boiler(q)) = eta_boiler;
            A(t,p3(q)) = 1;
            A(t,z3(q)) = min(Gamma(3,1),1);
            
            
            b(t) = C_t(q);
            q = q + 1;
            
        elseif t>15*time+6+3*time+11+2*time && t<=15*time+6+3*time+11+3*time
            % Vincolo addizionale Equazione di Bilancio TERMICO (<)
            A(t,p3(q)) = -1;
            A(t,z3(q)) = -1;
            b(t) = -epsilon*C_t(q);
            q = q + 1;
            
        elseif t>15*time+6+3*time+11+3*time && t<=15*time+6+4*time+11+3*time
            % Equazione di bilancio TERMICO (segno disequazione: >)
            A(t,P_chp(q)) = -eta_t_chp/eta_e_chp;
            A(t,P_stc(q)) = -1;
            A(t,P_hp(q)) = -COP;
            A(t,P_tes(q)) = -1;
            A(t,P_boiler(q)) = -eta_boiler;
            A(t,p3(q)) = -1;
            A(t,z3(q)) = -min(Gamma(3,1),1);
            
            b(t) = -C_t(q);
            q = q + 1;
            
        elseif t>15*time+6+4*time+11+3*time && t<=15*time+6+4*time+11+4*time
            % Vincolo addizionale Equazione di Bilancio TERMICO (>)
            A(t,p3(q)) = 1;
            A(t,z3(q)) = 1;
            b(t) = epsilon*C_t(q);
            q = q + 1;
            
        elseif t>15*time+6+4*time+11+4*time && t<=15*time+6+5*time+11+4*time
            % Equazione di bilancio GAS (segno disequazione: <)
            A(t,P_chp(q)) = -eta_t_chp/(PCI_gas*eta_e_chp);
            A(t,P_boiler(q)) = -1/PCI_gas;
            A(t,G_grid(q)) = 1;
            A(t,p4(q)) = 1;
            A(t,z4(q)) = min(Gamma(4,1),1);
            
            b(t) = C_g(q);
            q = q + 1;
            
        elseif t>15*time+6+5*time+11+4*time && t<=15*time+6+5*time+11+5*time
            % Vincolo addizionale Equazione di Bilancio GAS (<)
            A(t,p4(q)) = -1;
            A(t,z4(q)) = -1;
            b(t) = -epsilon*C_g(q);
            q = q + 1;
            
        elseif t>15*time+6+5*time+11+5*time && t<=15*time+6+6*time+11+5*time
            % Equazione di bilancio GAS (segno disequazione: >)
            A(t,P_chp(q)) = eta_t_chp/(PCI_gas*eta_e_chp);
            A(t,P_boiler(q)) = 1/PCI_gas;
            A(t,G_grid(q)) = -1;
            A(t,p4(q)) = -1;
            A(t,z4(q)) = -min(Gamma(4,1),1);
            
            b(t) = -C_g(q);
            q = q + 1;
            
        elseif t>15*time+6+6*time+11+5*time && t<=15*time+6+6*time+11+6*time
            % Vincolo addizionale Equazione di Bilancio GAS (>)
            A(t,p4(q)) = 1;
            A(t,z4(q)) = 1;
            b(t) = epsilon*C_g(q);
            q = q + 1;

        elseif t>15*time+6+6*time+11+6*time && t<=15*time+6+7*time+11+6*time
            % Equazione di bilancio COOLING (segno disequazione: <)
            A(t,P_hp(q)) = EER;
            A(t,p5(q)) = 1;
            A(t,z5(q)) = min(Gamma(5,1),1);
            
            b(t) = C_c(q);
            q = q + 1;
            
        elseif t>15*time+6+7*time+11+6*time && t<=15*time+6+7*time+11+7*time
            % Vincolo addizionale Equazione di Bilancio COOLING (<)
            A(t,p5(q)) = -1;
            A(t,z5(q)) = -1;
            b(t) = -epsilon*C_c(q);
            q = q + 1;
            
        elseif t>15*time+6+7*time+11+7*time && t<=15*time+6+8*time+11+7*time
            % Equazione di bilancio COOLING (segno disequazione: >)
            A(t,P_hp(q)) = -EER;
            A(t,p5(q)) = -1;
            A(t,z5(q)) = -min(Gamma(5,1),1);
            
            b(t) = -C_c(q);
            q = q + 1;
            
        elseif t>15*time+6+8*time+11+7*time && t<=15*time+6+8*time+11+8*time
            % Vincolo addizionale Equazione di Bilancio COOLING (>)
            A(t,p5(q)) = 1;
            A(t,z5(q)) = 1;
            b(t) = epsilon*C_c(q);
            q = q + 1;
        
        elseif t>15*time+6+8*time+11+8*time && t<=15*time+6+8*time+11+9*time
            % Equazione pannelli fotovoltaici PV (<)
            A(t,P_pv(q)) = 1;
            A(t,S_pv) = -epsilon*eta_pv*eta_bos*irr(q);
            A(t,p6(q)) = 1;
            A(t,z6(q)) = min(Gamma(6,1),1);
            q = q + 1;
        
        elseif t>15*time+6+8*time+11+9*time && t<=15*time+6+8*time+11+10*time
            % Equazione addizionale pannelli fotovoltaici PV (<)
            A(t,p6(q)) = -1;
            A(t,z6(q)) = -1;
            b(t) = -epsilon*irr(q);
            q = q + 1;
        
        elseif t>15*time+6+8*time+11+10*time && t<=15*time+6+8*time+11+11*time
            % Equazione pannelli fotovoltaici PV (>)
            A(t,P_pv(q)) = -1;
            A(t,S_pv) = epsilon*eta_pv*eta_bos*irr(q);
            A(t,p6(q)) = -1;
            A(t,z6(q)) = -min(Gamma(6,1),1);
            q = q + 1;
            
        elseif t>15*time+6+8*time+11+11*time && t<=15*time+6+8*time+11+12*time
            % Equazione addizionale pannelli fotovoltaici PV (>)
            A(t,p6(q)) = 1;
            A(t,z6(q)) = 1;
            b(t) = epsilon*irr(q);
            q = q + 1;
            
       elseif t>15*time+6+8*time+11+12*time && t<=15*time+6+8*time+11+13*time
            % Equazione pannelli fotovoltaici STC (<)
            A(t,P_stc(q)) = 1;
            A(t,S_stc) = -epsilon*eta_stc*irr(q);
            A(t,p7(q)) = 1;
            A(t,z7(q)) = min(Gamma(7,1),1);
            q = q + 1;
        
        elseif t>15*time+6+8*time+11+13*time && t<=15*time+6+8*time+11+14*time
            % Equazione addizionale pannelli fotovoltaici STC (<)
            A(t,p7(q)) = -1;
            A(t,z7(q)) = -1;
            b(t) = -epsilon*irr(q);
            q = q + 1;
        
        elseif t>15*time+6+8*time+11+14*time && t<=15*time+6+8*time+11+15*time
            % Equazione pannelli fotovoltaici STC (>)
            A(t,P_stc(q)) = -1;
            A(t,S_stc) = epsilon*eta_stc*irr(q);
            A(t,p7(q)) = -1;
            A(t,z7(q)) = -min(Gamma(7,1),1);
            q = q + 1;
            
        elseif t>15*time+6+8*time+11+15*time && t<=15*time+6+8*time+11+16*time
            % Equazione addizionale pannelli fotovoltaici STC (>)
            A(t,p7(q)) = 1;
            A(t,z7(q)) = 1;
            b(t) = epsilon*irr(q);
            q = q + 1;
            
        end
        
        if q > time
            q = 1;
        end
    end
    
    
    
    %---- Equality Constraints ----
    q = 1;
    for t=1:G_eq
        if t<=time 
            % Equazioni SoC ESS
            Aeq(t,SoC_ess(q)) = -1/(eta*dT);
            if (q-1) == 0
                Aeq(t,SoC_ess_0) = 1/(eta*dT);
            else
                Aeq(t,SoC_ess(q-1)) = 1/(eta*dT);
            end
            Aeq(t,Pn_ess) = 1;
            q = q + 1;
            
        elseif t>time && t<=2*time
            % Equazione SoC TES
            Aeq(t,SoC_tes(q)) = -1/(eta*dT);
            if (q-1) == 0
                Aeq(t,SoC_tes_0) = 1/(eta*dT);
            else
                Aeq(t,SoC_tes(q-1)) = 1/(eta*dT);
            end
            Aeq(t,Pn_tes) = 1;
            q = q + 1;
            
%         elseif t>2*time && t<=3*time
%             % Equazione pannelli fotovoltaici PV
%             Aeq(t,P_pv(q)) = 1;
%             Aeq(t,S_pv) = -eta_pv*eta_bos*irr(q);
%             q = q + 1;
%             
%         elseif t>3*time && t<=4*time
%             % Equazione pannelli fotovoltaici STC
%             Aeq(t,P_stc(q)) = 1;
%             Aeq(t,S_stc) = -eta_stc*irr(q);
%             q = q + 1;
            
        elseif t>2*time && t<=4*time+1
            % Equazione SoC iniziale ESS
            Aeq(t,SoC_ess_0) = 1;
            Aeq(t,E_ess) = -.1;
            
        elseif t>2*time+1 && t<=2*time+2
            % Equazione SoC iniziale TES
            Aeq(t,SoC_tes_0) = 1;
            Aeq(t,E_tes) = -.1;
            
        elseif t>2*time+2 && t<=2*time+3
            % Equazione di bilancio iniziale e finale ESS
            Aeq(t,SoC_ess_0) = 1;
            Aeq(t,SoC_ess(end)) = -1;
            
        elseif t>2*time+3 && t<=2*time+4
            % Equazione di bilancio iniziale e finale TES
            Aeq(t,SoC_tes_0) = 1;
            Aeq(t,SoC_tes(end)) = -1;
            
        end
        
        if q > time
            q = 1;
        end
    end
    
    %---- Funzione Obiettivo ----
    % f = @(x) [x(Cost), x(Emiss)];
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
