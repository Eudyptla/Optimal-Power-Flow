function Lxx =HessFcn_QReserve(x,lambda,cost_mult)

% Hessian of objective
% Hessian of nonlinear inequality and equality constraint 

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA branch 
% ====================================================================
% Initialization
% ====================================================================
define_constants;          % MATPOWER use only

% Problem dimensions
nb = numAll(1);            % number of buses
nl = numAll(2);            % number of branches
ng = numAll(3);            % number of dispatchable injections
nvars = numAll(4);         % number of variables
nxtra = numAll(5);         % number of variables excluding Va and Vm


iVa = 1:nb;                      % index of Va
iVm = nb+1:2*nb;                 % index of Vm
iPg = 2*nb+1:2*nb+ng;            % index of Pg
iQg = 2*nb+ng+1:2*nb+2*ng;       % index of Qg


% ====================================================================
%  Set  constrained lines
% ====================================================================

% Data at operating point
Va = x(iVa);                       % in rad
Vm = x(iVm);                       % in p.u.
Pg = x(iPg);                       % in p.u.
Qg = x(iQg);                       % in p.u.

gen(:, PG) = Pg * baseMVA;         % in MW
gen(:, QG) = Qg * baseMVA;         % in MVAr

V = Vm .* exp(1j * Va);            % voltage phasors at operating point

% ====================================================================
% Hessian of objective
% ====================================================================

%d2f=[];

% ====================================================================
% Hessian of nonlinear equality constraint
% ====================================================================
nlam = length(lambda.eqnonlin)/2;
lamP = lambda.eqnonlin(1:nlam);
lamQ = lambda.eqnonlin((1:nlam)+nlam);
[Gpaa, Gpav, Gpva, Gpvv] = d2Sbus_dV2(Ybus, V, lamP);
[Gqaa, Gqav, Gqva, Gqvv] = d2Sbus_dV2(Ybus, V, lamQ);

% Evaluate Hessian of power balance constraints
d2G = [
    real([Gpaa Gpav; Gpva Gpvv]) + imag([Gqaa Gqav; Gqva Gqvv]) sparse(2*nb, nxtra);
    sparse(nxtra, 2*nb + nxtra)
];

% Evaluate Hessian of equality constraints
%d2H = [];

% Hessian is the matrix of second derivatives of the Lagrangian
Lxx = d2G ;
