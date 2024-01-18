function [Pn_hp, Pn_ess, Pn_tes, Pn_boiler, Pn_chp, S_pv, S_stc, E_ess,...
    E_tes, P_hp, P_ess, P_tes, P_boiler, P_chp, P_pv, P_stc, P_sold,...
    P_purch, SoC_ess, SoC_tes, G_grid, SoC_ess_0, SoC_tes_0, P_ess_0,...
    P_tes_0, Cost, Emiss, z1, z2, z3, z4, z5, z6, z7,...
    p1, p2, p3, p4, p5, p6, p7, lam1, lam2, lam3, lam4, lam5, lam6, lam7, ...
    mu1, mu2, mu3, mu4, mu5, mu6, mu7] = idx_variables(time, comps)

    % IDX_VARIABLES   Defines constants for variables.
    %   Example:
    %
    %   [Pn_hp, Pn_ess, Pn_tes, Pn_boiler, Pn_chp, S_pv, S_stc, E_ess,...
    %       E_tes, P_hp, P_ess, P_tes, P_boiler, P_chp, P_pv, P_stc, P_sold,...
    %       P_purch, SoC_ess, SoC_tes, G_grid, SoC_ess_0, SoC_tes_0, P_ess_0,...
    %       P_tes_0, Cost, Emiss, z1, z2, z3, z4, z5, ...
    %       p1, p2, p3, p4, p5, lam1, lam2, lam3, lam4, lam5, lam6, lam7, ...
    %       mu1, mu2, mu3, mu4, mu5, mu6, mu7] = idx_variables(time, comps)

    %% define variables indexes DETERMINISTIC ONLY
    Pn_hp = 1;
    Pn_ess = 2;
    Pn_tes = 3;
    Pn_boiler = 4;
    Pn_chp = 5;
    S_pv = 6;
    S_stc = 7;
    E_ess = 8;
    E_tes = 9;

    P_hp = [10:time+10-1]';
    P_ess = [time+10:2*time+10-1]';
    P_tes = [2*time+10:3*time+10-1]';
    P_boiler = [3*time+10:4*time+10-1]';
    P_chp = [4*time+10:5*time+10-1]';
    P_pv = [5*time+10:6*time+10-1]';
    P_stc = [6*time+10:7*time+10-1]';


    P_sold = [7*time+10:8*time+10-1]';
    P_purch = [8*time+10:9*time+10-1]';
    SoC_ess = [9*time+10:10*time+10-1]';
    SoC_tes = [10*time+10:11*time+10-1]';
    G_grid = [11*time+10:12*time+10-1]';

    SoC_ess_0 = 12*time+10;
    SoC_tes_0 = 12*time+10+1;
    P_ess_0 = 12*time+10+2;
    P_tes_0 = 12*time+10+3;

    Cost = 12*time+10+4;
    Emiss = 12*time+10+5;


    %% define variables indexes ROBUST ONLY

    % variables polyhedral set
    z1 = 12*time+10+6;
    z2 = [12*time+10+6+1:12*time+10+6+time];
    z3 = [12*time+10+6+time+1:12*time+10+6+2*time];
    z4 = [12*time+10+6+2*time+1:12*time+10+6+3*time];
    z5 = [12*time+10+6+3*time+1:12*time+10+6+4*time];
    z6 = [12*time+10+6+4*time+1:12*time+10+6+5*time];
    z7 = [12*time+10+6+5*time+1:12*time+10+6+6*time];

    % variables box set
    p1 = [12*time+10+6+6*time+1:12*time+10+6+9*time+(length(comps)+2)+1];
    p2 = [12*time+10+6+9*time+2+(length(comps)+2):12*time+10+6+10*time+(length(comps)+2)+1];
    p3 = [12*time+10+6+10*time+2+(length(comps)+2):12*time+10+6+11*time+(length(comps)+2)+1];
    p4 = [12*time+10+6+11*time+2+(length(comps)+2):12*time+10+6+12*time+(length(comps)+2)+1];
    p5 = [12*time+10+6+12*time+2+(length(comps)+2):12*time+10+6+13*time+(length(comps)+2)+1];
    p6 = [12*time+10+6+13*time+2+(length(comps)+2):12*time+10+6+14*time+(length(comps)+2)+1];
    p7 = [12*time+10+6+14*time+2+(length(comps)+2):12*time+10+6+15*time+(length(comps)+2)+1];
    
    % Variables for asymetrical set
    lam1 = [12*time+10+6+15*time+(length(comps)+2)+2:12*time+10+6+18*time+2*(length(comps)+2)+2];
    lam2 = [12*time+10+6+18*time+2*(length(comps)+2)+2+1:12*time+10+6+19*time+2*(length(comps)+2)+2];
    lam3 = [12*time+10+6+19*time+2*(length(comps)+2)+2+1:12*time+10+6+20*time+2*(length(comps)+2)+2];
    lam4 = [12*time+10+6+20*time+2*(length(comps)+2)+2+1:12*time+10+6+21*time+2*(length(comps)+2)+2];
    lam5 = [12*time+10+6+21*time+2*(length(comps)+2)+2+1:12*time+10+6+22*time+2*(length(comps)+2)+2];
    lam6 = [12*time+10+6+22*time+2*(length(comps)+2)+2+1:12*time+10+6+23*time+2*(length(comps)+2)+2];
    lam7 = [12*time+10+6+23*time+2*(length(comps)+2)+2+1:12*time+10+6+24*time+2*(length(comps)+2)+2];
    
    mu1 = [12*time+10+6+24*time+2*(length(comps)+2)+3:12*time+10+6+27*time+3*(length(comps)+2)+3];
    mu2 = [12*time+10+6+27*time+3*(length(comps)+2)+3+1:12*time+10+6+28*time+3*(length(comps)+2)+3];
    mu3 = [12*time+10+6+28*time+3*(length(comps)+2)+3+1:12*time+10+6+29*time+3*(length(comps)+2)+3];
    mu4 = [12*time+10+6+29*time+3*(length(comps)+2)+3+1:12*time+10+6+30*time+3*(length(comps)+2)+3];
    mu5 = [12*time+10+6+30*time+3*(length(comps)+2)+3+1:12*time+10+6+31*time+3*(length(comps)+2)+3];
    mu6 = [12*time+10+6+31*time+3*(length(comps)+2)+3+1:12*time+10+6+32*time+3*(length(comps)+2)+3];
    mu7 = [12*time+10+6+32*time+3*(length(comps)+2)+3+1:12*time+10+6+33*time+3*(length(comps)+2)+3];
end
