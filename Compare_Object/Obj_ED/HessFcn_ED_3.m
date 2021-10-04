function Lxx = HessFcn_ED_3(x,lambda,cost_mult)

% Hessian of objective
% Hessian of nonlinear inequality and nonlinear equality constraint 

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus Yf Yt baseMVA gencost idx_G idx_L

% ====================================================================
% Initialization
% ====================================================================
define_constants;          % MATPOWER use only

% Problem dimensions
nb = numAll(1);            % number of buses
nl = numAll(2);            % number of branches
ng = numAll(3);            % number of dispatchable injections
nvars = 2*nb + 2*ng;           % number of variables
nxtra = nvars - 2*nb;          % number of variables excluding Va and Vm

% Index of Va, Vm, Pg and Qg
iVa = 1:nb;                % index of Va
iVm = nb+1:2*nb;           % index of Vm
iPg = 2*nb+1:2*nb+ng;      % index of Pg
iQg = 2*nb+ng+1:2*nb+2*ng; % index of Qg

% ====================================================================
% Hessian of objective
% ====================================================================
% Polynomial cost of P 
coeff2 = gencost(:,5);    % ($/MW^2)

% Evaluate d2f
d2f = sparse(nvars,nvars);
d2f(iPg,iPg) = 2 * diag(coeff2) * baseMVA * baseMVA;

d2f = d2f * cost_mult;

% ====================================================================
% Hessian of nonlinear equality constraint
% ====================================================================
Va = x(iVa);                       % in rad
Vm = x(iVm);                       % in p.u.

% Update voltage phasors
V = Vm .* exp(1j * Va);

nlam = length(lambda.eqnonlin) / 2;
lamP = lambda.eqnonlin(1:nlam);
lamQ = lambda.eqnonlin((1:nlam)+nlam);
[Gpaa, Gpav, Gpva, Gpvv] = d2Sbus_dV2(Ybus, V, lamP);
[Gqaa, Gqav, Gqva, Gqvv] = d2Sbus_dV2(Ybus, V, lamQ);

% Evaluate Hessian of power balance constraints
d2G = [
    real([Gpaa Gpav; Gpva Gpvv]) + imag([Gqaa Gqav; Gqva Gqvv]) sparse(2*nb, nxtra);
    sparse(nxtra, 2*nb + nxtra)
];


Lxx = d2f + d2G ;
