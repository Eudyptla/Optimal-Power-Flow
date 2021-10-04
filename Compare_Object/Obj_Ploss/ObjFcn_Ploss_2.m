function [f,df] = ObjFcn_Ploss_2(x,bus,branch)
% f:   Objective function (scalar) as the first output
% df:  Gradient (vector) as the second output

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus Yf Yt baseMVA  gencost idx_G idx_L 

% ====================================================================
% Initialization
% ====================================================================
define_constants;           % MATPOWER use only

% Problem dimensions
nb = numAll(1);            % number of buses
nl = numAll(2);            % number of branches
ng = numAll(3);            % number of generators
nvars = 2*nb + 2*ng;           % number of variables
nxtra = nvars - 2*nb;          % number of variables excluding Va and Vm

% Index of Va, Vm, Pg and Qg
iVa = 1:nb;                % index of Va
iVm = nb+1:2*nb;           % index of Vm
iPg = 2*nb+1:2*nb+ng;      % index of Pg
iQg = 2*nb+ng+1:2*nb+2*ng; % index of Qg
% ====================================================================
% Update optimization variables x = [Va;Vm;Pg;Qg] during iteration
% ====================================================================
Va = x(iVa);               % in rad
Vm = x(iVm);               % in p.u.
V = Vm.*exp(1j*Va);
% ====================================================================
% Evaluate objective function
% ====================================================================
% construct complex bus voltage vector
bus(:, VA) = Va * 180 / pi;        % in degree
bus(:, VM) = Vm;                   % in p.u.
[loss, fchg, tchg, dloss_dV] = get_losses(baseMVA,bus,branch);
clear fchg tchg 
real_Ploss = real(loss);
f = sum(real_Ploss);

% ====================================================================
% Evaluate gradient vector
% Gradient is the vector of first derivatives of the objective function
% ====================================================================
df = sparse(nvars,1);
dPloss_dVa = (dloss_dV.a')*(ones(nl,1));
dPloss_dVm = (dloss_dV.m')*(ones(nl,1));
df(iVa,:) = real(dPloss_dVa);
df(iVm,:) = real(dPloss_dVm);

