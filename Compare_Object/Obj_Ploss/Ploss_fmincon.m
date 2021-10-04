function Obj_Ploss = Ploss_fmincon(mpc)

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus Yf Yt baseMVA  gencost idx_G idx_L

% ====================================================================
% Load data from MATPOWER format
% ====================================================================
define_constants;            % MATPOWER use only
nb = numAll(1);            % number of buses
nl = numAll(2);            % number of branches
ng = numAll(3);            % number of generators
nvars = 2*nb + 2*ng;           % number of variables
nxtra = nvars - 2*nb;          % number of variables excluding Va and Vm

% ====================================================================
% Initialization
% ====================================================================

% Convert external to internal indexing
mpc = ext2int(mpc);        % MATPOWER Fun: ext2int

% Distribute inputs to outputs
[ bus, gen, branch] =...
             deal(mpc.bus, mpc.gen, mpc.branch);

% ====================================================================
% Initial points in p.u. values
% ====================================================================
Va = bus(:, VA) * (pi / 180);   % in rad
Vm = bus(:, VM);                % in p.u.
Pg = gen(:, PG) / baseMVA;      % in p.u.
Qg = gen(:, QG) / baseMVA;      % in p.u.

x0 = [Va;Vm;Pg;Qg];             % column vector of initial points

% ====================================================================
% Scalar objective function
% ====================================================================
f_fcn = @(x) ObjFcn_Ploss_2(x,bus,branch);

% ====================================================================
% Bound constraints (Eq.6.4)
% ====================================================================
refs = find(bus(:, BUS_TYPE) == REF);
Vau = Inf(nb, 1);                  % upper bounds of Va 
Val = -Vau;                        % lower bounds of Va

% Voltage angle bounds for swing bus (Eq.6.10)
Vau(refs) = Va(refs);              % upper bounds of Va_swing
Val(refs) = Va(refs);              % lower bounds of Va_swing

% Voltage magnitude bounds (Eq.6.11)
Vmu = bus(:, VMAX);                % upper bounds of Vm
Vml = bus(:, VMIN);                % lower bounds of Vm

% Gen real power bounds (Eq.6.12)
PGmax = gen(:, PMAX) / baseMVA;    % upper bounds of Pg
PGmin = gen(:, PMIN) / baseMVA;    % lower bounds of Pg

% Gen reactive power bounds (Eq.6.13)
QGmax = gen(:, QMAX) / baseMVA;    % upper bounds of Qg
QGmin = gen(:, QMIN) / baseMVA;    % lower bounds of Qg

% Variable limits (Eq.6.4)
xmax = [Vau;Vmu;PGmax;QGmax];      % upper bounds of variables
xmin = [Val;Vml;PGmin;QGmin];      % lower bounds of variables

% ====================================================================
% Linear constraints:
% Linear inequality constraints (A*x <= b)
% Linear equality constraints (Aeq*x = beq)
% ====================================================================
A = [];
b = [];

Aeq = [];
beq = [];

% ====================================================================
%  Nonlinear constraints: (Eq.6.2, Eq.6.3)
%  Nonlinear equality constraints:   power balance equations (Eq.4.2, Eq.4.3)
%  Nonlinear inequality constraints: branch flow limits      (Eq.6.7, Eq.6.8)
% ====================================================================
gh_fcn = @(x) NlinFcn_Ploss_2(x,bus,gen);

% ====================================================================
% Hessian matrix
% ====================================================================  
cost_mult = 1;   % For interior-point algorithm
hess_fcn = @(x,lambda) HessFcn_Ploss_2(x,lambda,cost_mult,bus,branch);

% ====================================================================
% Set nondefault options
% ====================================================================
  
  options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
                         'GradObj','on',...      % Gradient of fun
                       'GradConstr','on',...   % Gradients of nonlcon
                       'Hessian','user-supplied','HessFcn',hess_fcn,... % Hessian of fun and nonlcon;
                       'TolCon',5e-6,'TolX',1e-4,'TolFun',1e-4);
% ====================================================================
% Find minimum of constrained nonlinear multivariable function
% ====================================================================
[x,fval,exitflag,Output] = fmincon(f_fcn, x0, A, b, Aeq, beq,...
                                          xmin, xmax, gh_fcn, options);

% ====================================================================
% Summarize the results
% ====================================================================
success = (exitflag > 0);

% Index of Va, Vm, Pg and Qg
iVa = 1:nb;                        % index of Va
iVm = nb+1:2*nb;                   % index of Vm
iPg = 2*nb+1:2*nb+ng;              % index of Pg
iQg = 2*nb+ng+1:2*nb+2*ng;         % index of Qg

% Optimization Variables
Va = x(iVa);                       % in rad
Vm = x(iVm);                       % in p.u.
Pg = x(iPg);                       % in p.u.
Qg = x(iQg);                       % in p.u.

% Update some bus data
bus(:, VA) = Va * 180 / pi;        % in degree
bus(:, VM) = Vm;                   % in p.u.

% Update some generator data
gen(:, PG) = Pg * baseMVA;         % in MW
gen(:, QG) = Qg * baseMVA;         % in MVAr
gen(:, VG) = Vm(gen(:, GEN_BUS));  % in p.u.

% Compute branch flows
V = Vm .* exp(1j * Va);            % Voltage phasors

Sf = V(branch(:, F_BUS)) .* conj(Yf * V);  % cplx pwr at "from" bus, p.u. (Eq. 3.15)
St = V(branch(:, T_BUS)) .* conj(Yt * V);  % cplx pwr at "to" bus, p.u.   (Eq. 3.16)

% Update some branch data
branch(:, PF) = real(Sf) * baseMVA;        % P injected at "from" bus end (MW)
branch(:, QF) = imag(Sf) * baseMVA;        % Q injected at "from" bus end (MVAr)
branch(:, PT) = real(St) * baseMVA;        % P injected at "to" bus end (MW)
branch(:, QT) = imag(St) * baseMVA;        % Q injected at "to" bus end (MVAr)

% Update the results (bus data, generator data, and branch data)
results = mpc;
[results.bus, results.branch, results.gen, ...
    results.x, results.fval,results.success] = ...  % New items
        deal(bus, branch, gen, x, fval, success);

% Convert internal to external bus numbering
results = int2ext(results);                % MATPOWER Fun: int2ext

% ====================================================================
% Calculate Original Pcost
% ====================================================================
coeff2 = gencost(:,5);     % ($/MW^2)
coeff1 = gencost(:,6);     % ($/MW)
coeff0 = gencost(:,7);     % ($)
Pg = results.gen(:,2);
Qg = results.gen(:,3);
Pcost= sum( coeff2.*(Pg.^2)+coeff1.*Pg+coeff0);
% ====================================================================
% Calculate Original Ploss
% ====================================================================
Ploss = sum(real(get_losses(results)));% MATPOWER Fun: get_losses
% ====================================================================
% Calculate Original Voltage Profile
% ====================================================================
Vm = results.bus(:,VM);               % in p.u.
Voltage_Deviations=Vm(idx_L)-ones(nb-ng,1);
Voltage_Profile= sum(abs(Voltage_Deviations));
% ====================================================================
% Calculate Original L-index 
% ====================================================================
Va = results.bus(:,VA)*(pi/180);               % in rad
V = Vm .* exp(1j * Va);
Vgen = V(idx_G);
Vload = V(idx_L);
Y_LL = Ybus(idx_L,idx_L);
Y_LG = Ybus(idx_L,idx_G);
F_LG = (-1)*( Y_LL\Y_LG );
Lj =ones(nb-ng,1)-(F_LG*Vgen)./Vload;
Lindex = abs(Lj);
Lindex_max = max(Lindex);
% ====================================================================
% Calculate Original Q Reserve 
% ====================================================================
QReserve= sum(results.gen(:,3));
% ====================================================================
% Calculate Original Qcost
% ====================================================================
PQ=Pg.*Qg;
S2=Pg.^2+Qg.^2;

Qcost= sum( coeff2.*((PQ./S2).^2).*(Qg.^2)+coeff1.*(PQ./S2).*Qg+coeff0);


 Obj_Ploss = [Pcost,Ploss,Voltage_Profile,Lindex_max,QReserve...
                  Qcost];
