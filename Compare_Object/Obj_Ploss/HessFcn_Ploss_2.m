function Lxx = HessFcn_Ploss_2(x,lambda,cost_mult,bus,branch)

% Hessian of objective
% Hessian of nonlinear inequality and nonlinear equality constraint 

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA Yf Yt  gencost idx_G idx_L

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

Va = x(iVa);                       % in rad
Vm = x(iVm);                       % in p.u.

% Update voltage phasors
V = Vm .* exp(1j * Va);
% ====================================================================
% Hessian of objective
% ====================================================================
% Polynomial cost of P 


% Evaluate d2f

% construct complex bus voltage vector
% create map of external bus numbers to bus indices
% create map of external bus numbers to bus indices
i2e = bus(:, BUS_I);
e2i = sparse(max(i2e), 1);
e2i(i2e) = (1:size(bus, 1))';
% parameters
Cf = sparse(1:nl, e2i(branch(:, F_BUS)), branch(:, BR_STATUS), nl, nb);
Ct = sparse(1:nl, e2i(branch(:, T_BUS)), branch(:, BR_STATUS), nl, nb);

tap = ones(nl, 1);                              % default tap ratio = 1 for lines
xfmr = find(branch(:, TAP));                    % indices of transformers
tap(xfmr) = branch(xfmr, TAP);                  % include transformer tap ratios
tapVm = tap ; 
tapVa = exp(1j*pi/180 * branch(:, SHIFT));% add phase shifters
CfVm = spdiags(1 ./ tapVm, 0, nl, nl) * Cf;
CfVa = spdiags(1 ./ tapVa, 0, nl, nl) * Cf;

Ysc = 1 ./ (branch(:, BR_R) - 1j * branch(:, BR_X)); %branch attidance
d2S_Va = (CfVa*exp(1j*x(iVa))).*(Ct*exp(-1j*x(iVa)));

d2S_dVaf2 = Ysc.*(CfVm*x(iVm)).*(Ct*x(iVm)).*(d2S_Va+conj(d2S_Va));

d2S_dVafVat = -d2S_dVaf2;

d2S_dVatVaf = d2S_dVafVat;

d2S_dVafVmf = 1j*Ysc.*((CfVm*ones(nb,1)).*(Ct*x(iVm)).*(conj(d2S_Va)-d2S_Va));
d2S_dVmfVaf = d2S_dVafVmf;

d2S_dVafVmt = 1j*Ysc.*(CfVm*x(iVm)).*(conj(d2S_Va)-d2S_Va);
 
d2S_dVat2 = d2S_dVaf2;

d2S_dVatVmf = -d2S_dVmfVaf;

d2S_dVatVmt = -d2S_dVafVmt;

d2S_dVmf2 = 2*Ysc.*(CfVm*ones(nb,1));
d2S_dVmt2 = 2*Ysc.*(Ct*ones(nb,1));

d2S_dVmfVmt = Ysc.*(-d2S_Va-conj(d2S_Va));
d2S_dVmtVmf = d2S_dVmfVmt;
d2f = sparse(nvars,nvars);
d2f(iVa,iVa) = Cf'*diag(d2S_dVaf2)*Cf+Cf'*diag(d2S_dVafVat)*Ct+Ct'*diag(d2S_dVatVaf)*Cf+Ct'*diag(d2S_dVat2)*Ct;
d2f(iVa,iVm) = Cf'*diag(d2S_dVafVmf)*Cf+Cf'*diag(d2S_dVafVmt)*Ct+Ct'*diag(d2S_dVatVmt)*Ct+Ct'*diag(d2S_dVatVmf)*Cf;
d2f(iVm,iVa) = d2f(iVa,iVm)';
d2f(iVm,iVm) = Cf'*diag(d2S_dVmf2)*Cf+Ct'*diag(d2S_dVmt2)*Ct+Cf'*diag(d2S_dVmfVmt)*Ct+Ct'*diag(d2S_dVmtVmf)*Cf;
d2f = baseMVA*d2f;
d2f= real(d2f);

d2f = d2f * cost_mult;

% ====================================================================
% Hessian of nonlinear equality constraint
% ====================================================================


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

% ====================================================================
% Hessian of objective, nonlinear inequality, and nonlinear equality constraints
% Hessian is the matrix of second derivatives of the Lagrangian
% ====================================================================
Lxx = d2f + d2G  ;