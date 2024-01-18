function [vargout] = var_indexing(T, comps, ro, PLOT, symmetric)
    % CREATE VARIABLES
    % Create the variables Nominal powers
    n = length(comps);
    for i=1:(n-2)
        variables{i,1} = strcat('Pn_',comps{i,1});
    end
    % Creo le variabili S_pv e S_stc
    variables{i+1,1} = strcat('S_',comps{end-1,1});
    variables{i+2,1} = strcat('S_',comps{end,1});
    variables{i+3,1} = strcat('E_',comps{2,1});
    variables{i+4,1} = strcat('E_',comps{3,1});
    
    k = length(variables)+1;
    for i=1:length(comps)
        for j=1:T
            variables{k,1} = strcat('P_',comps{i,1},'_',num2str(j));
            k = k + 1;
        end
    end
    
    for j=1:T
        variables{k,1} = strcat('P_sold_',num2str(j));
        k = k + 1;
    end
    
    for j=1:T
        variables{k,1} = strcat('P_purch_',num2str(j));
        k = k + 1;
    end
    for j=1:T
        variables{k,1} = strcat('SoC_ess_',num2str(j));
        k = k + 1;
    end
    for j=1:T
        variables{k,1} = strcat('SoC_tes_',num2str(j));
        k = k + 1;
    end
    for j=1:T
        variables{k,1} = strcat('G_grid_',num2str(j));
        k = k + 1;
    end
    
    
    variables{k,1} = 'SoC_ess_0';
    k = k + 1;
    variables{k,1} = 'SoC_tes_0';
    k = k + 1;
    variables{k,1} = 'P_ess_0';
    k = k + 1;
    variables{k,1} = 'P_tes_0';
    k = k + 1;
    variables{k,1} = 'Cost';
    k = k + 1;
    variables{k,1} = 'Emiss';
    k = k + 1;
    

    if ro
        % Variabile additional per Uncertainty Set (US) poliedrico vincolo
        % objective function
        variables{k,1} = 'z1';
        k = k + 1;
        % Variabile additional per US poliedrico vincoli ELECTRICITY demand
        for i=1:T
            variables{k,1} = strcat('z2_',num2str(i));
            k = k + 1;
        end
        % Variabile additional per US poliedrico vincoli THERMAL demand
        for i=1:T
            variables{k,1} = strcat('z3_',num2str(i));
            k = k + 1;
        end
        % Variabile additional per US poliedrico vincoli GAS demand
        for i=1:T
            variables{k,1} = strcat('z4_',num2str(i));
            k = k + 1;
        end
        % Variabile additional per US poliedrico vincoli COOLING demand
        for i=1:T
            variables{k,1} = strcat('z5_',num2str(i));
            k = k + 1;
        end
        % Variabile additional per US poliedrico vincoli IRR PV
        for i=1:T
            variables{k,1} = strcat('z6_',num2str(i));
            k = k + 1;
        end
        % Variabile additional per US poliedrico vincoli IRR STC
        for i=1:T
            variables{k,1} = strcat('z7_',num2str(i));
            k = k + 1;
        end
        
        % Variabile additional per US box vincolo funzione obiettivo
        for f=1:(T*3+(length(comps)+2)+1)
            variables{k,1} = strcat('p1_',num2str(f));
            k = k + 1;
        end
        
        % Variabile additional per US box vincoli DEMANDS
        L_eqRobuste = 7;
        for g=2:L_eqRobuste
            for f=1:T
                variables{k,1} = strcat(strcat(strcat('p',num2str(g)),...
                '_'),num2str(f));
                k = k + 1;
            end
        end
        
        if ~symmetric       
            % Variabile additional "lambda" e "mu" per US box vincoli
            % DEMANDS (asymmetric constraints)
            L_eqRobuste = 7;
            for g=2:L_eqRobuste
                for f=1:T
                    % variabile v
                    variables{k,1} = strcat(strcat(strcat('v',num2str(g)),...
                    '_'),num2str(f));
                    k = k + 1;
                    
                    % variabile mu
                    variables{k,1} = strcat(strcat(strcat('mu',num2str(g)),...
                    '_'),num2str(f));
                    k = k + 1;
                end
            end
        end
    end

    vargout{1} = variables;
    
    % Creation of symbolic variables
    if PLOT
        N = length(variables);
        parfor v = 1:N
            variables_sym(v,1) = sym([strcat(variables{v}, '_sym')]);
        end
        vargout{2} = variables_sym;
%     else
%         vargout{2} = 0;
    end


end

