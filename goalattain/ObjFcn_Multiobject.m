
function [f,df] = ObjFcn_Multiobject(x,bus)

% f:   Objective function (scalar) as the first output
% df:  Gradient (vector) as the second output

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA branch idx_G idx_L gencost 

% ====================================================================
% Initialization
% ====================================================================
define_constants;           % MATPOWER use only

% Problem dimensions
nb = numAll(1);            % number of buses
nl = numAll(2);            % number of branches
ng = numAll(3);            % number of generators
nvars = numAll(4);         % number of variables
nxtra = numAll(5);         % number of variables excluding Va and Vm

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
f = sparse(3,1);
%f1 ED
Pg = x(iPg) * baseMVA;     % in MW

% Coefficients of n-th order polynomial cost
coeff2 = gencost(:,5);     % ($/MW^2)
coeff1 = gencost(:,6);     % ($/MW)
coeff0 = gencost(:,7);     % ($)

% Polynomial cost of Pg 
f(1,1) = (Pg') * diag(coeff2) * Pg + (coeff1') * Pg + sum(coeff0); % ($)

%f2 Ploss
% construct complex bus voltage vector
bus(:, VA) = Va * 180 / pi;        % in degree
bus(:, VM) = Vm;                   % in p.u.
[loss, fchg, tchg, dloss_dV] = get_losses(baseMVA,bus,branch);
clear fchg tchg 
real_Ploss = real(loss);
f(2,1) = sum(real_Ploss);

%f3 Lindex
f(3,1) =x(nvars);

% ====================================================================
% Evaluate gradient vector
% Gradient is the vector of first derivatives of the objective function
% ====================================================================

df = sparse(nvars,3);
%df1 ED
df(iPg,1) = 2 * coeff2 .* Pg * baseMVA + coeff1 * baseMVA;

%df2 Ploss
dPloss_dVa = (dloss_dV.a')*(ones(nl,1));
dPloss_dVm = (dloss_dV.m')*(ones(nl,1));
df(iVa,2) = real(dPloss_dVa);
df(iVm,2) = real(dPloss_dVm);

%df3 Lindex
df(:,3) = sparse(nvars,1,1,nvars,1);
