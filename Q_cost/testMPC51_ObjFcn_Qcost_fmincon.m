clear all; close all;

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA branch 

% ====================================================================
% Input the name of the test case
% ====================================================================
% testcase ='case6ww'; lam =1.0; 
%                 Pcost       Lmax      Qcost 
%                 ______    ________    ______
%     Original    3189.5    0.094362    1507.9
%     Optimal     3174.5    0.096612    1330.8

%testcase ='case6ww'; lam =1.4; 
%                 Pcost      Lmax      Qcost 
%                 ______    _______    ______
%     Original    4364.6    0.16252    2292.9
%     Optimal     4281.2    0.15215    2127.5

%testcase = 'case9'; lam =1.0;
%                 Pcost      Lmax      Qcost 
%                 ______    _______    ______
%     Original    5438.3     0.1673    1128.8
%     Optimal     5299.1    0.13848    1085.6

%testcase = 'case9'; lam =2.0;
%                 Pcost     Lmax      Qcost 
%                 _____    _______    ______
%     Original    17394    0.45017    2542.9
%     Optimal     16955    0.33633    1662.7

% testcase = 'case14';lam =1.0; 
%初值界外 只有兩台發電機工作，三台不工作。

%                 Pcost       Lmax      Qcost 
%                 ______    ________    ______
%     Original      8166    0.076707    308.27
%     Optimal     8169.9    0.079917    182.09

%testcase = 'case14';lam =1.5;
%                 Pcost     Lmax      Qcost 
%                 _____    _______    ______
%     Original    14909    0.12112    1334.8
%     Optimal     13754    0.12874    388.84

%  testcase='case30';lam=1.0; 
%                 Pcost       Lmax      Qcost 
%                 ______    ________    ______
%     Original    584.27    0.055286    106.43
%     Optimal     600.82    0.052478    15.648

%testcase='case30';lam=1.7; 
%                 Pcost       Lmax      Qcost 
%                 ______    ________    ______
%     Original    1217.5     0.10292    247.32
%     Optimal     1184.2    0.095071    206.79

%testcase = 'case39';lam =1.0;
%                 Pcost     Lmax      Qcost 
%                 _____    _______    ______
%     Original    45077    0.20098    308.84
%     Optimal     41891     0.2024    237.59

%testcase = 'case39';lam =1.1;
%                 Pcost     Lmax      Qcost 
%                 _____    _______    ______
%     Original    55251    0.22767    522.85
%     Optimal     51607    0.22321    339.54

%testcase='case57'; lam =1.0;
%                 Pcost     Lmax      Qcost 
%                 _____    _______    ______
%     Original    51348     0.3099    1875.7
%     Optimal     42215    0.30109    871.25

%testcase='case57'; lam =1.4;
%                 Pcost     Lmax      Qcost 
%                 _____    _______    ______
%     Original    99998    0.58586      3817
%     Optimal     70796    0.48166    3068.6

%testcase='case118'; lam=1.0;
%                   Pcost         Lmax      Qcost
%                 __________    ________    _____
%     Original    1.3121e+05    0.069311     6229
%     Optimal      1.301e+05    0.069696    27.62

%testcase='case118'; lam=1.5; 
%                   Pcost        Lmax      Qcost 
%                 __________    _______    ______
%     Original    2.5493e+05    0.11353     13615
%     Optimal     2.1797e+05    0.10274    3558.1

%testcase='case300'; lam=1.0;
%                   Pcost        Lmax      Qcost
%                 __________    _______    _____
%     Original    7.2475e+05     0.4202    56513
%     Optimal      7.228e+05    0.40731    35833

%testcase='case1354pegase';lam =1.0;

% testcase ='case2383wp'; lam=1.0 ;
%diverge
%testcase ='case3012wp';lam=1.0; %多台發電機在同一條bus上

% ====================================================================
% Load data from MATPOWER format
% ====================================================================
define_constants;            % MATPOWER use only
mpc = loadcase( testcase );  % MATPOWER Fun: loadcase
mpc = ext2int(mpc);
mpopt=mpoption('pf.enforce_q_lims',1,'out.all',0);

% ====================================================================
% Global variables
% ====================================================================
nb = size(mpc.bus, 1);         % number of buses
nl = size(mpc.branch, 1);      % number of branches
ng = size(mpc.gen, 1);         % number of generators
nvars = 2*nb + 2*ng;       % number of variables
nxtra = nvars - 2*nb;      % number of variables excluding Va and Vm
idx_G = find((mpc.bus(:,BUS_TYPE) == 3) | (mpc.bus(:,BUS_TYPE) == 2));    % index of slack &PV bus
idx_L = find(mpc.bus(:,BUS_TYPE) == 1 );                              % index of PQ bus
numAll = [nb;nl;ng;nvars;nxtra;];
[Ybus,Yf,Yt] = makeYbus(mpc); % MATPOWER Fun: makeYbus
% ====================================================================
% Power flow solution from MATPOWER
% ====================================================================
mpc.gen(:,[PG QG])=mpc.gen(:,[PG QG])*lam;
mpc.bus(:,[PD QD])=mpc.bus(:,[PD QD])*lam;
for i=1:ng
    if mpc.gen(i,PG)>mpc.gen(i,PMAX)
        mpc.gen(i,PG)=mpc.gen(i,PMAX);
    end
end
results_0 = runpf( mpc,mpopt);
% Build admittance matrix

% ====================================================================
% Initialization
% ====================================================================
t0 = clock;                % start timer
% Distribute inputs to outputs
mpc = ext2int(results_0);
[baseMVA, bus, gen, branch,gencost] =...
             deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch,mpc.gencost);
% Convert external to internal indexing
% ====================================================================
% Initial points in p.u. values
% ====================================================================
Va = bus(:, VA) * (pi / 180);   % in rad
Vm = bus(:, VM);                % in p.u.
Pg = gen(:, PG) / baseMVA;      % in p.u.
Qg = gen(:, QG) / baseMVA;      % in p.u.
x0 = [Va;Vm;Pg;Qg];             % column vector of initial points

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
Vml = 0.5*bus(:, VMIN);                % lower bounds of Vm

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
% Scalar objective function
% ====================================================================
f_fcn = @(x) ObjFcn_Qcost(x,gencost,gen,bus);

% ====================================================================
%  Nonlinear constraints: (Eq.6.2, Eq.6.3)
%  Nonlinear equality constraints:   power balance equations (Eq.4.2, Eq.4.3)
% ====================================================================
gh_fcn = @(x) NlinFcn_Qcost(x,bus,gen);

% ====================================================================
% Hessian matrix
% ====================================================================  
% cost_mult = 1;   % For interior-point algorithm
% hess_fcn = @(x,lambda) HessFcn_Qcost(x,lambda,cost_mult,gencost);

% ====================================================================
% Set nondefault options
% ====================================================================
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
%                        'GradObj','on',...      % Gradient of fun
%                        'GradConstr','on',...   % Gradients of nonlcon
%                        'DerivativeCheck','on',... %Check gradient dont use on case>100 
%                        'TolCon',5e-6,'TolX',1e-4,'TolFun',1e-4,...
%                        'PlotFcn',@optimplotfval); 
%  options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
%                        'GradObj','on',...      % Gradient of fun
%                        'GradConstr','on',...   % Gradients of nonlcon
%                        'FinDiffType','central',...
%                        'Hessian','user-supplied',...% Hessian of fun and nonlcon;
%                        'HessFcn',hess_fcn,... % user set;
%                        'TolCon',5e-6,'TolX',1e-4,'TolFun',1e-4,...
%                        'PlotFcn',@optimplotfval);
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
                       'GradObj','on',...      % Gradient of fun
                       'GradConstr','on',...   % Gradients of nonlcon
                       'FinDiffType','central',... %interior point only
                       'TolCon',5e-6,'TolX',1e-4,'TolFun',1e-4,...
                       'PlotFcn',@optimplotfval);
%  options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
%                        'FinDiffType','central',... %interior point only
%                        'TolCon',5e-6,'TolX',1e-4,'TolFun',1e-4,...
%                        'PlotFcn',@optimplotfval);
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
% Print out power flow results in MATPOWER output format
% ====================================================================
et = etime(clock, t0);  % Compute elapsed time
results.et = et;
%printpf(results);       % MATPOWER Fun: printpf
% ====================================================================
% Calculate Original Qcost
% ====================================================================
Pg_0 = results_0.gen(:,2);
Qg_0 = results_0.gen(:,3);
PQ_0=Pg_0.*Qg_0;
S2_0=Pg_0.^2+Qg_0.^2;
coeff2 = gencost(:,5);     % ($/MW^2)
coeff1 = gencost(:,6);     % ($/MW)
coeff0 = gencost(:,7);     % ($)

Qcost_0= sum( coeff2.*((PQ_0./S2_0).^2).*(Qg_0.^2)+coeff1.*(PQ_0./S2_0).*Qg_0+coeff0);
% ====================================================================
% Calculate Original L-index 
% ====================================================================
Va_0 = results_0.bus(:,VA)*(pi/180);               % in rad
Vm_0 = results_0.bus(:,VM);               % in p.u.
V_0 = Vm_0 .* exp(1j * Va_0);
Vgen_0 = V_0(idx_G);
Vload_0 = V_0(idx_L);
Y_LL = Ybus(idx_L,idx_L);
Y_LG = Ybus(idx_L,idx_G);
F_LG = (-1)*( Y_LL\Y_LG );
Lj_0 =ones(nb-ng,1)-(F_LG*Vgen_0)./Vload_0;
Lindex_0 = abs(Lj_0);
Lindex_0_max = max(Lindex_0);
% ====================================================================
% Calculate Original Pcost
% ====================================================================
Pcost_0= sum( coeff2.*(Pg_0.^2)+coeff1.*Pg_0+coeff0);
% ====================================================================
% Calculate OPF Qcost
% ====================================================================
PQ = gen(:, PG).*gen(:, QG);
S2= gen(:, PG).^2 + gen(:, QG).^2;
Qcost_N= sum(coeff2.*((PQ./S2).^2).*(gen(:, QG).^2)+coeff1.*(PQ./S2).*gen(:, QG)+coeff0);
% ====================================================================
% Calculate OPF Pcost
% ====================================================================
Pcost_N= sum( coeff2.*((Pg * baseMVA).^2)+coeff1.*Pg * baseMVA+coeff0);

% ====================================================================
% Calculate OPF L-index 
% ====================================================================

Vgen = V(idx_G);
Vload = V(idx_L);
Lj =ones(nb-ng,1)-(F_LG*Vgen)./Vload;
Lindex_O = abs(Lj);
Lindex_max = max(Lindex_O);

% ====================================================================
% Print out the results of Voltage Profile
% ====================================================================
Pcost = [Pcost_0;Pcost_N];
Qcost = [Qcost_0;Qcost_N];
Lmax = [Lindex_0_max;Lindex_max];

rowName = {'Original';'Optimal'};
% Tabulate the results 
table(Pcost,Lmax,Qcost,...
    'rowNames',rowName)

