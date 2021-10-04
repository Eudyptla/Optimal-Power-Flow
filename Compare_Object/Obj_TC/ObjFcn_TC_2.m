function [f,df] = ObjFcn_TC_2(x)

% f:   Objective function (scalar) as the first output
% df:  Gradient (vector) as the second output
% d2f: Hessian (matrix) as the third output

% ====================================================================
% Declare variables as global
% ====================================================================
global  numAll Ybus Yf Yt baseMVA  gencost idx_G idx_L

% ====================================================================
% Initialization
% ====================================================================
define_constants;            % MATPOWER use only
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

Pg = x(iPg) * baseMVA;     % in MW
Qg = x(iQg) * baseMVA;     % in MVAR

S2 = Pg.^2+ Qg.^2;
PQ = Pg.* Qg;

% Coefficients of n-th order polynomial cost
coeff2 = gencost(:,5);     % ($/MW^2)
coeff1 = gencost(:,6);     % ($/MW)
coeff0 = gencost(:,7);     % ($)

f=sum( coeff2.*((PQ./S2).^2).*(Qg.^2)+coeff1.*(PQ./S2).*Qg+coeff0...
    + coeff2.*(Pg.^2)+coeff1.*Pg);

% ====================================================================
% Evaluate gradient vector
% Gradient is the vector of first derivatives of the objective function
% ====================================================================

df=sparse(nvars,1);
Pg2 = Pg.^2;
Qg2 = Qg.^2;
Pg3 = Pg.^3;
Qg3 = Qg.^3;
Pg4 = Pg.^4;
Qg4 = Qg.^4;
Qg5 = Qg.^5;
S22 = S2.^2;
S23 = S2.^3;


df(iPg,1)=(2*coeff2.*(Pg.*Qg4./S22 - 2*(Pg3.*Qg4./S23))...
         +coeff1.*(Qg2./S2 - 2*(Pg2.*Qg2./S22)))*baseMVA...
         +2*coeff2.*Pg*baseMVA +coeff1*baseMVA;

df(iQg,1)=(4*coeff2.*(Pg2.*Qg3./S22-Pg2.*Qg5./S23)...
          +2*coeff1.*(Pg.*Qg./S2-Pg.*Qg3./S22))*baseMVA;

