function [f,df] = ObjFcn_ED_3(x)

% f:   Objective function (scalar) as the first output
% df:  Gradient (vector) as the second output
% d2f: Hessian (matrix) as the third output

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus  Yf Yt baseMVA gencost idx_G idx_L

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
% Evaluate objective function
% ====================================================================
Pg = x(iPg) * baseMVA;     % in MW

% Coefficients of n-th order polynomial cost
coeff2 = gencost(:,5);     % ($/MW^2)
coeff1 = gencost(:,6);     % ($/MW)
coeff0 = gencost(:,7);     % ($)

% Polynomial cost of Pg 
f = (Pg') * diag(coeff2) * Pg + (coeff1') * Pg + sum(coeff0); % ($)

% ====================================================================
% Evaluate gradient vector
% Gradient is the vector of first derivatives of the objective function
% ====================================================================
 df = sparse(nvars,1);
 df(iPg,1) = 2 * coeff2 .* Pg * baseMVA + coeff1 * baseMVA;
 

