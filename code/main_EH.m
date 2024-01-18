% Clear the command window, clear all variables, and close all figures
clc, clear, close all

% Define the number of time intervals
time = 24;

% --- Opening files and reading data from an Excel file --- %
% Read load data from the specified Excel sheet and range
Load = xlsread('Dati Hub.xlsx','Foglio1','B1:B8761');
% Read thermal data
Thermal = xlsread('Dati Hub.xlsx','Foglio1','C1:C8761');
% Read gas data
Gas = xlsread('Dati Hub.xlsx','Foglio1','D1:D8761');
% Read cooling data
Cool = xlsread('Dati Hub.xlsx','Foglio1','E1:E8761');
% Read irradiation data and scale it
Irr = xlsread('Dati Hub.xlsx','Foglio1','F1:F8761');
Irr = 0.001 * Irr;

% --- Reshaping variables --- %
% Reshape data into matrices with 'time' rows, facilitating further processing
L = reshape(Load, [time, size(Load, 1) / time]);
T = reshape(Thermal, [time, size(Thermal, 1) / time]);
G = reshape(Gas, [time, size(Gas, 1) / time]);
C = reshape(Cool, [time, size(Cool, 1) / time]);
I = reshape(Irr, [time, size(Irr, 1) / time]);

% --- kMeans clustering to identify a typical day for each energy type --- %
% Apply kMeans with one centroid (typical day) for each energy type
[Cl_k, C_l] = kmeans(L', 1);
[Ct_k, C_t] = kmeans(T', 1);
[Cg_k, C_g] = kmeans(G', 1);
[Cc_k, C_c] = kmeans(C', 1);
[irr_k, irr] = kmeans(I', 1);

% Define the number of time frames (hours in a year)
T = time;
% List of devices/energy vectors in the configuration
comps = {'hp'; 'ess'; 'tes'; 'boiler'; 'chp'; 'pv'; 'stc'};

% Robust optimization flag
ro = true;
% Auxiliary variable for future developments (leave as true)
symmetric = true;
% Flag for printing output equations in the terminal
PLOT = false;

% Number of typical days (set to 1 here)
Ns = 1;
% Probability of a typical day occurring
ws = 365 / 365;

% Economic parameters
r = 0.08; % Discount rate
n = [20; 12; 20; 15; 30; 30; 20]; % Useful life of components in years
N = 20; % Time horizon in years
a = (r * (1 + r).^n) ./ ((1 + r).^n - 1); % Actualization factor
R = (1 - ((N - 1) ./ n)); % Residual rate

% --- Optimization problem setup --- %
% If robust optimization
if ro
    [f, A, b, Aeq, beq, lb, ub, n_var, Cost, Emiss, variables] = robust_hub(T, Ns, ws, PLOT, comps, a, R, C_l, C_t, C_g, C_c, irr, symmetric);
else
    % If deterministic optimization
    [f, A, b, Aeq, beq, lb, ub, n_var, Cost, Emiss, variables] = deterministic_hub(T, Ns, ws, PLOT, comps, a, R, C_l, C_t, C_g, C_c, irr, obj, f1f2);
end

% --- Solve the linear programming problem --- %
% Objective: minimize costs
[x, fval, exitflag, output, lambda] = linprog(f, A, b, Aeq, beq, lb, ub);

% --- Write output file --- %
% Convert numeric array 'x' to a cell array
x_cell = num2cell(x);
% Create a table from 'variables' and 'x_cell'
output_table = table(variables, x_cell);
% Write the table to an Excel file
writetable(output_table, 'output.xlsx', 'WriteVariableNames', false);
