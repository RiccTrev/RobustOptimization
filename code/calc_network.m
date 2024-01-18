function [baseMVA, Base_kV, Branch, Bus, Slack, Gen, Vpu, Ybus_sp,...
    Yf_sp, Yt_sp, Ybus, V, Vm, Y, Ibus, Ibus_m,...
    Ibranch, Im] = calc_network(net, VM, VA, BASE_KV)

% Extract useful data from the network
baseMVA = net.baseMVA;
Base_kV = net.bus(:,BASE_KV);
Branch = net.branch;
Bus = net.bus;
Slack = find(Bus(:,2) == 3);
Gen = net.gen;

% Calculate voltages (p.u.)
% Vpu = VM * exp(1j * (pi/180) * VA);
Vpu = net.bus(:, VM) .* exp(sqrt(-1) * pi/180 * Bus(:, VA));
Vm = abs(Vpu);

% Find voltages on brancehs and admittance matrix
[~, ~, ~, ~, ~, Ysc, Vdrop, ~] = get_losses(baseMVA, net.bus, net.branch);
%%%% NOTE: La Ysc does not consider capacitive susceptance values but only of conductance and inductive susceptance

% Calculate the currents in the branches
Ibranch = Ysc .* Vdrop;
Im_b = abs(Ibranch) * 1e3;

[Ybus_sp, Yf_sp, Yt_sp] = makeYbus(baseMVA, Bus, Branch);
Ybus = full(Ybus_sp);
Yf = full(Yf_sp);
Yt = full(Yt_sp);
V = Vpu .* Base_kV; % Values in kV
Ibranch_Yf = Yf * V;
Ibranch_Yt = Yt * V;


Y = abs(Ybus);
Ibus = Ybus * V; % Values in kV
Ibus_m = abs(Ibus) * 1e3; % ATTENTION: values in A

% Calculate currents
for j=1:size(Branch,1)
    idx = Branch(j,1:2);
    %%% CAUTION: voltage values are in kV,
    %%% therefore current values will be in kA
    Ibranch(j,1) = Ybus(idx(1),idx(2)) * (V(idx(1))-V(idx(2)));
end

Im = abs(Ibranch) * 1e3; % WARNING: data transformed into A


end