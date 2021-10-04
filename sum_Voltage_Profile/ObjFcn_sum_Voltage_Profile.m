function [f,df] = ObjFcn_sum_Voltage_Profile(x)

% f:   Objective function (scalar) as the first output
% df:  Gradient (vector) as the second output
% d2f: Hessian (matrix) as the third output

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA branch idx_G idx_L

% ====================================================================
% Initialization
% ====================================================================
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

Vm = x(iVm);                       % in p.u.
% Update voltage phasors
Voltage_Dev=Vm(idx_L)-ones(nb-ng,1);
Voltage_Profile = abs(Voltage_Dev);
f=sum(Voltage_Profile);

% ====================================================================
% Evaluate gradient vector
% Gradient is the vector of first derivatives of the objective function
% ====================================================================
df = sparse(nvars,1);

negative=find(Voltage_Dev<0);
L= ones(nb-ng,1);
L(negative,1)=-1*L(negative,1);
df(idx_L+nb,1)=L;



