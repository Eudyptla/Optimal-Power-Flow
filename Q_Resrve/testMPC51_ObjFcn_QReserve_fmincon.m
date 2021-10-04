clear all; close all;

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA branch 

% ====================================================================
% Input the name of the test case
% ====================================================================
%testcase ='case6ww'; lam =1.0; 
%QR_0=179.94 QR_O=176.09 
% L_0=0.094362 L_O=0.094572
%%發電機輸出電壓不可調度
%testcase ='case6ww'; lam =1.4; 
%QR_0=diverge QR_O=290.49
%L_0=diverge L_O=0.14998

%testcase = 'case9'; lam =1.0;
%QR_0=34.88 QR_O=-13.386
%L_0=1.673 L_0=0.13337 

%testcase = 'case9'; lam =2.0;
%QR_0=371.47 QR_O=241.5
%L_0=0.45017 L_0=0.31866

%testcase = 'case14';lam =1.0; 
%初值界外 只有兩台發電機工作，三台不工作。 
%QR_0=81.862 QR_O=34.726
%L_0=0.076707 L_0=0.078362

%testcase = 'case14';lam =1.5;
%QR_0=diverge QR_O=84.377
%L_0=0.12112 L_0=0.12354

%testcase='case30';lam=1.0; 
%QR_0= 100.41 QR_O=95.573
%L_0=0.055286 L_0=0.049071

%testcase='case30';lam=1.7; 
%QR_0= 198.29 QR_O=188.41
%L_0=0.10292 L_0=0.089446

%testcase = 'case39';lam =1.0;
%QR_0=1247.7 QR_O=1095.8
%L_0=0.20098 L_0=0.1935

%testcase = 'case39';lam =1.1;
%QR_0=1709.9 QR_O=1479.8
%L_0=0.22767 L_0=0.21494

%testcase='case57'; lam =1.0;
%QR_0=321.08 QR_O=248.6
%L_0=0.3099 L_0=0.30362

% testcase='case57'; lam =1.4;
%QR_0= 648.48 QR_O=489.32
%L_0=0.58586 L_0=0.48067

% testcase='case118'; lam=1.0;
%QR_0=793.92 QR_O=-68.541
%L_0=0.069311 L_0=0.061469

%testcase='case118'; lam=1.5; 
%QR_0=2757.1 QR_O=792.89
%L_0=0.11353 L_0=0.096314

%testcase='case300'; lam=1.0;
%QR_0=7998.3 QR_O=5627.2
%L_0=0.4202 L_0=0.3932

%testcase='case1354pegase';lam =1.0;
%QR_0=19858 QR_O=9317.9
%L_0=0.40928 L_0=0.38736

testcase ='case2383wp'; lam=1.0 ;
%diverge
%testcase ='case3012wp';lam=1.0; %多台發電機在同一條bus上

% ====================================================================
% Load data from MATPOWER format
% ====================================================================
define_constants;            % MATPOWER use only
mpc = loadcase( testcase );  % MATPOWER Fun: loadcase
mpc = ext2int(mpc);
mpopt=mpoption('pf.enforce_q_lims',1);

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
[baseMVA, bus, gen, branch] =...
             deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);
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
f_fcn = @(x) ObjFcn_QReserve(x,QGmax);

% ====================================================================
%  Nonlinear constraints: (Eq.6.2, Eq.6.3)
%  Nonlinear equality constraints:   power balance equations (Eq.4.2, Eq.4.3)
% ====================================================================
gh_fcn = @(x) NlinFcn_QReserve(x,bus,gen);

% ====================================================================
% Hessian matrix
% ====================================================================  
cost_mult = 1;   % For interior-point algorithm
hess_fcn = @(x,lambda) HessFcn_QReserve(x,lambda,cost_mult);

% ====================================================================
% Set nondefault options
% ====================================================================
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
%                        'GradObj','on',...      % Gradient of fun
%                        'GradConstr','on',...   % Gradients of nonlcon
%                        'DerivativeCheck','on',... %Check gradient dont use on case>100 
%                        'TolCon',5e-6,'TolX',1e-4,'TolFun',1e-4,...
%                        'PlotFcn',@optimplotfval); 
 options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
                       'GradObj','on',...      % Gradient of fun
                       'GradConstr','on',...   % Gradients of nonlcon
                       'FinDiffType','central',...
                       'Hessian','user-supplied',...% Hessian of fun and nonlcon;
                       'HessFcn',hess_fcn,... % user set;
                       'TolCon',5e-6,'TolX',1e-4,'TolFun',1e-4,...
                       'PlotFcn',@optimplotfval);
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
%                        'GradObj','on',...      % Gradient of fun
%                        'GradConstr','on',...   % Gradients of nonlcon
%                        'FinDiffType','central',... %interior point only
%                        'TolCon',5e-6,'TolX',1e-4,'TolFun',1e-4,...
%                        'PlotFcn',@optimplotfval);
%  options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
%                        'FinDiffType','central',... %interior point only
%                        'TolCon',5e-6,'TolX',1e-4,'TolFun',1e-4,...
%                        'PlotFcn',@optimplotfval);
% ====================================================================
% Find minimum of constrained nonlinear multivariable function
% ====================================================================
[x,fval,exitflag,Output] = fmincon(f_fcn, x0, A, b, Aeq, beq,...
                                          xmin, xmax, gh_fcn, options);
% [x,fval,exitflag,Output,lambda] = fmincon(f_fcn, x0, A, b, Aeq, beq,...
%                                           xmin, xmax, gh_fcn, options);

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
printpf(results);       % MATPOWER Fun: printpf
% ====================================================================
% Calculate Original QReserve
% ====================================================================
QReserve_0= sum(results_0.gen(:,3));
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
% Calculate OPF L-index 
% ====================================================================

Vgen = V(idx_G);
Vload = V(idx_L);
Lj =ones(nb-ng,1)-(F_LG*Vgen)./Vload;
Lindex = abs(Lj);
Lindex_max = max(Lindex);

% ====================================================================
% Print out the results of Voltage Profile
% ====================================================================


% Tabulate the results of max L-index
table(QReserve_0,fval*baseMVA,...
    'VariableNames',{'Q_Reserve_Original' 'Q_Reserve_Optimal'}) 

table(Lindex_0_max,Lindex_max,...
    'VariableNames',{'Lindex_Original' 'Lindex_Optimal'}) 