clear, clc, close all

n_variables = 2;
syms x1 x2
form = input('insert case under study: ');

switch form
    case 0
        % --- Problema deterministico ---
        ub = Inf(n_variables,1);
        lb = zeros(n_variables,1);
        A = [10 20; 6 8];
        b = [140; 72];

        Aeq = [];
        beq = [];
        f = [-8 -12];

        f*[x1; x2]
        A*[x1; x2]
        b

        %Calcoliamo il problema
        options = optimoptions('linprog','Algorithm','dual-simplex');
        [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,options);
        x
        abs(fval)
        
    case 1
        % --- Problema deterministico opposto 
        % (max con vincoli upper-equal) ---
        ub = Inf(n_variables,1);
        lb = zeros(n_variables,1);
        A = [-10 -20; -6 -8];
        b = [-140; -72];

        Aeq = [];
        beq = [];
        f = [8 12];

        f*[x1; x2]
        A*[x1; x2]
        b

        %Calcoliamo il problema
        options = optimoptions('linprog','Algorithm','dual-simplex');
        [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,options);
        x
        abs(fval)
        
    case 2
        % --- Problema deterministico con vincoli di equality---
        ub = Inf(n_variables,1);
        lb = zeros(n_variables,1);
        A = [];
        b = [];

        Aeq = [10 20; 6 8];
        beq = [140; 72];
        f = [-8 -12];

        f*[x1; x2]
        Aeq*[x1; x2]
        beq

        %Calcoliamo il problema
        options = optimoptions('linprog','Algorithm','dual-simplex');
        [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,options);
        x
        abs(fval)
        
    case 3
        % --- Problema deterministico con vincoli equality splittati in 
        % inequality constraints---
        ub = Inf(n_variables,1);
        lb = zeros(n_variables,1);
        A = [10 20; -10 -20; 6 8; -6 -8];
        b = [140; -140; 72; -72];

        Aeq = [];
        beq = [];
        f = [-8 -12];

        f*[x1; x2]
        A*[x1; x2]
        b

        %Calcoliamo il problema
        options = optimoptions('linprog','Algorithm','dual-simplex');
        [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,options);
        x
        abs(fval)
        
    case 4
        % --- Problema robusto LHS ---
        syms z1 z2 p11 p12 p21 p22
        n_variables = n_variables + 6;

        J = 2;
        epsilon = 0.1;
        method = 0;
        flag = 0;
        [Gamma, B] = g_def(J,epsilon,method,flag);


        ub = Inf(n_variables,1);
        lb = zeros(n_variables,1);
        A = [10 20 Gamma 0 1 1 0 0;
             epsilon*10 0 -1 0 -1 0 0 0;
             0  epsilon*20 -1 0 0 -1 0 0;
             6 8 0 Gamma 0 0 1 1;
             epsilon*6 0 0 -1 0 0 -1 0;
             0 epsilon*8 0 -1 0 0 0 -1];
        b = [140; 0; 0; 72; 0; 0];

        Aeq = [];
        beq = [];
        f = [-8 -12 0 0 0 0 0 0];

        f*[x1; x2; z1; z2; p11; p12; p21; p22]
        A*[x1; x2; z1; z2; p11; p12; p21; p22]

        %Calcoliamo il problema
        options = optimoptions('linprog','Algorithm','dual-simplex');
        [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,options);
        x
        abs(fval)
        
    case 5
        % --- Problema robusto LHS opposto 
        % (max con vincoli upper-equal) ---
        syms z1 z2 p11 p12 p21 p22
        n_variables = n_variables + 6;

        J = 2;
        epsilon = 0.1;
        method = 0;
        flag = 0;
        [Gamma, B] = g_def(J,epsilon,method,flag);


        ub = Inf(n_variables,1);
        lb = zeros(n_variables,1);
        A = [-10 -20 -Gamma 0 -1 -1 0 0;
             -epsilon*10 0 1 0 1 0 0 0;
             0  -epsilon*20 1 0 0 1 0 0;
             -6 -8 0 -Gamma 0 0 -1 -1;
             -epsilon*6 0 0 1 0 0 1 0;
             0 -epsilon*8 0 1 0 0 0 1];
        b = [-140; 0; 0; -72; 0; 0];

        Aeq = [];
        beq = [];
        f = [8 12 0 0 0 0 0 0];

        f*[x1; x2; z1; z2; p11; p12; p21; p22]
        A*[x1; x2; z1; z2; p11; p12; p21; p22]
        b

        %Calcoliamo il problema
        options = optimoptions('linprog','Algorithm','dual-simplex');
        [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,options);
        x
        abs(fval)
        
    case 6
        % --- Problema robusto LHS + RHS ---
        syms z1 z2 p11 p12 p21 p22 p01 p02
        n_variables = n_variables + 8;

        J = 3;
        epsilon = .1;
        epsilon_0 = .1;
        method = 0;
        flag = 0;
        [Gamma, B] = g_def(J,epsilon,method,flag);


        ub = Inf(n_variables,1);
        lb = zeros(n_variables,1);
        A = [10         20          Gamma   0       1   1   0   0  1    0;
             epsilon*10 0           -1      0       -1  0   0   0  0    0;
             0          epsilon*20  -1      0       0   -1  0   0  0    0;
             0          0           -1      0       0   0   0   0  -1   0;
             6          8           0       Gamma   0   0   1   1  0    1;
             epsilon*6  0           0       -1      0   0   -1  0  0    0;
             0          epsilon*8   0       -1      0   0   0   -1 0    0;
             0          0           0       -1      0   0   0   0  0   -1];
        b = [140; 0; 0; -epsilon*140; 72; 0; 0; -epsilon*72];

        Aeq = [];
        beq = [];
        f = [-8 -12 0 0 0 0 0 0 0 0];

        f*[x1; x2; z1; z2; p11; p12; p21; p22; p01; p02]
        A*[x1; x2; z1; z2; p11; p12; p21; p22; p01; p02]
        b

        %Calcoliamo il problema
        options = optimoptions('linprog','Algorithm','dual-simplex');
        [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,options);
        x
        abs(fval)
        
    case 7
        % --- Problema robusto LHS con vincoli opposti E con vincoli 
        % equality splittati in inequality constraints------
        syms z1 z2 p11 p12 p21 p22 % in più? p31 p32 p41 p42 (z3 z4)
        n_variables = n_variables + 6;

        J = 2;
        epsilon = .1;
        epsilon_0 = .1;
        method = 0;
        flag = 0;
        [Gamma, B] = g_def(J,epsilon,method,flag);


        ub = Inf(n_variables,1);
        lb = zeros(n_variables,1);
        A = [10          20         Gamma     0      1      1    0   0;
             epsilon*10  0           -1       0     -1      0    0   0;
             0           epsilon*20  -1       0      0     -1    0   0;
            -10         -20        -Gamma     0     -1     -1    0   0;
            -epsilon*10 0             1       0      1      0    0   0;
             0          -epsilon*20   1       0      0      1    0   0;
             6           8            0     Gamma    0      0    1   1;
             epsilon*6   0            0      -1      0      0   -1   0;
             0           epsilon*8    0      -1      0      0    0  -1;
             -6          -8           0    -Gamma    0      0   -1  -1;
            -epsilon*6   0            0       1      0      0    1   0;
             0          -epsilon*8    0       1      0      0    0   1];
        b = [140; 0; 0; -140; 0; 0; 72; 0; 0; -72; 0; 0];

        Aeq = [];
        beq = [];
        f = [-8 -12 0 0 0 0 0 0];

        f*[x1; x2; z1; z2; p11; p12; p21; p22]
        A*[x1; x2; z1; z2; p11; p12; p21; p22]
        b

        %Calcoliamo il problema
        options = optimoptions('linprog','Algorithm','dual-simplex');
        [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,options);
        x
        abs(fval)
        
    case 8
        % --- Problema robusto LHS con vincoli opposti E con vincoli 
        % equality splittati in inequality constraints ma 
        % con variabili slack ---
        syms z1 z2 p11 p12 p21 p22 s1 s2
        n_variables = n_variables + 8;

        J = 2;
        epsilon = .1;
        epsilon_0 = .1;
        method = 0;
        flag = 0;
        [Gamma, B] = g_def(J,epsilon,method,flag);


        ub = Inf(n_variables,1);
        lb = zeros(n_variables,1);
        % Le variabili ausiliare s1 ed s2 devono allargare l'intervallo e
        % non restringerlo, quindi vanno con il segno cambiato
        A = [10          20         Gamma     0      1      1    0   0   1   0;
             epsilon*10  0           -1       0     -1      0    0   0   0   0;
             0           epsilon*20  -1       0      0     -1    0   0   0   0;
            -10         -20        -Gamma     0     -1     -1    0   0   1   0;
            -epsilon*10 0             1       0      1      0    0   0   0   0;
             0          -epsilon*20   1       0      0      1    0   0   0   0;
             6           8            0     Gamma    0      0    1   1   0   1;
             epsilon*6   0            0      -1      0      0   -1   0   0   0;
             0           epsilon*8    0      -1      0      0    0  -1   0   0;
             -6          -8           0    -Gamma    0      0   -1  -1   0   1;
            -epsilon*6   0            0       1      0      0    1   0   0   0;
             0          -epsilon*8    0       1      0      0    0   1   0   0];
        b = [140; 0; 0; -140; 0; 0; 72; 0; 0; -72; 0; 0];

        Aeq = [];
        beq = [];
        f = [-8 -12 0 0 0 0 0 0 0 0];

        f*[x1; x2; z1; z2; p11; p12; p21; p22; s1; s2]
        A*[x1; x2; z1; z2; p11; p12; p21; p22; s1; s2]
        b

        %Calcoliamo il problema
        options = optimoptions('linprog','Algorithm','dual-simplex');
        [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,options);
        x
        abs(fval)
end


%% Test on function handle

clc, clear, close all

a = [1:5];
b = [1:3];

f = @(x) myobj(a,b,x);
% f2 = @(x) b*x;

nvars = length(a)+length(b);
lb = zeros(nvars,1);
ub = 100*ones(nvars,1);
x = gamultiobj(f,nvars,[],[],[],[],lb,ub)

function f = myobj(a,b,x)
    f(1) = a*x(1:5);
    f(2) = b*x(3:5);
end