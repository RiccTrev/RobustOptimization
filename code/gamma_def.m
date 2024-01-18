function [Gamma, B] = gamma_def(J,epsilon,method,flag)
    %definisce per l'ottimizzazione robusta il valore di gamma: definisce
    %dimensioni dell'uncertainty set. 
    Tmin = 0.0;
    q = .01;        % Resolution of Ti 
                    % (ATTENZIONE: la risoluzione incide in maniera molto
                    % forte sul costo finale a fronte di tempi di calcolo 
                    % molto più bassi)

    Ti=[Tmin:q:J]';


    t = 1;
    for T=Tmin:q:J
        if method == 1

            v = (T+J)/2;
            v_f = floor(v);
            mu = v - v_f;
            C1 = 0;

            for k=v_f+1:J
                C1 = C1 + nchoosek(J,k);
            end
            %Questo metodo sembra che parta da 1 e non da zero
            if T >= 1
                Bound1(t,J) = (1/(2^J))*((1-mu)*nchoosek(J,v_f) + C1);
            else
                Bound1(t,J) = 1;
            end
            t = t + 1;

        elseif method == 0

            Bound2(t,J) = exp(-(T^2)/(2*J));    
            t = t + 1;

        elseif method == 2

            if mod(ceil(T)+J,2) == 0
                n_i = J+1-sign(T);
            else
                n_i = J-sign(T);
            end
            v = (T+n_i)/2;

            V1 = 0;
            for k=ceil(v):n_i
                V1 = V1 + nchoosek(n_i,k);
            end
            Bound3(t,J) = (1/2^n_i)*V1;
            t = t + 1;
            
        elseif method == 3
            
            Bound4(t,J) = exp(-(T^2)/(2));    
            t = t + 1;
            
        end
    end


    if method == 1
        Bound = Bound1;
    elseif method == 0
        Bound = Bound2;
    elseif method == 2
        Bound = Bound3;
    elseif method == 3
        Bound = Bound4;
    end

    try
        G = find(Bound(:,end) <= epsilon);
        Gamma = (G(1)-1)*q;
    catch
        % disp('No valid Gamma for given epsilon')
        % disp('Return max value')
        Gamma = J;
    end

    B = Bound(:,end);

    if flag
        figure(1)
        hold on; grid on
            plot(Ti,Bound(:,end),'-b','LineWidth',1.5)
%             plot(Ti,Bound1(:,end),'-b','LineWidth',1.5)
%             plot(Ti,Bound2(:,end),'-r','LineWidth',1.5)
%             plot(Ti,Bound3(:,end),'-r','LineWidth',1.5)
            xlabel('\Gamma_{i}')
            ylabel('Probability of constraint violation (\epsilon) [%]')
            legend('B(J_{i}, \Gamma_{i})','e^{{-\Gamma_{i}^2}/_{2J_{i}}}','Metodo Combinatoriale')
    end



end

