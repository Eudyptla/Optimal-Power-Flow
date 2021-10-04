function Lxx =HessFcn_Qcost(x,lambda,cost_mult,gencost)

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
Pg = x(iPg) * baseMVA;     % in MW
Qg = x(iQg) * baseMVA;     % in MVAR
S2 = Pg.^2+ Qg.^2;
PQ = Pg.* Qg;
Pg2 = Pg.^2;
Qg2 = Qg.^2;
Pg3 = Pg.^3;
Qg3 = Qg.^3;
Pg4 = Pg.^4;
Qg4 = Qg.^4;
Qg5 = Qg.^5;
Qg6 = Qg.^6;
S22 = S2.^2;
S23 = S2.^3;
S24 = S2.^4;
% Coefficients of n-th order polynomial cost
coeff2 = gencost(:,5);     % ($/MW^2)
coeff1 = gencost(:,6);     % ($/MW)


V = Vm .* exp(1j * Va);            % voltage phasors at operating point

% ====================================================================
% Hessian of objective
% ====================================================================

d2f=sparse(nvars,nvars);
d2f_dPg2 = (2*coeff2.*(Pg4./S22-10*Pg2.*Qg4./S23 + 12*Pg4.*Qg4./S24)...
          +2*coeff1.*(4*Pg3.*Qg2./S23-3*Pg.*Qg2./S22))*baseMVA*baseMVA;

d2f_dPgdQg = (8*coeff2.*(Pg.*Qg3./S22-Pg.*Qg5./S23-2*Pg3.*Qg3./S23+3*Pg3.*Qg5./S24)...
             +2*coeff1.*(Qg./S2 - 5*Pg.*Qg2./S22+4*Pg.*Qg4./S23))*baseMVA*baseMVA;
         
d2f_dQg2 = (12*coeff2.*(Pg2.*Qg2./S22 - 3*Pg2.*Qg4./S23 + 2*Pg2.*Qg6./S24)...
           +2*coeff1.*(Pg./S2 +5*Pg.*Qg2./S22 +4*Pg.*Qg4./S23))*baseMVA*baseMVA;  

d2f(iPg,iPg) = diag(d2f_dPg2);
d2f(iPg,iQg) = diag(d2f_dPgdQg);
d2f(iQg,iPg) = d2f(iPg,iQg);
d2f(iQg,iQg) = diag(d2f_dQg2);
d2f= d2f*cost_mult;

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
Lxx = d2f+d2G ;
