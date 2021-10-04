clear all; close all;

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA branch idx_G idx_L

% ====================================================================
% Input the name of the test case
% ====================================================================
lam =1;

%testcase = 'case5Lab';   %Error original=0.0066257 L2sum = 0.0066524
%testcase = 'case9Lab';   %L2sum=0.0066061
%testcase = 'case26Lab';   %L2sum=0.085586

%testcase ='case6ww';  %L2sum = 0.020042 pf=0.020055 
%testcase = 'case9';  %L2sum =0.049348
%testcase = 'case14'; %Error pf sol not feasible answer 0.020405

%testcase = 'case39'; %L2sum = 0.58787
%testcase='case57';    %L2sum = 0.91849

%testcase='case300';   %L2sum = 4.0615 52sec

%testcase='case1354pegase'; %L2sum = 10.922 2335sec diverge

%testcase ='case2383wp';  %L2sum = 2.0859 6508sec diverge
%testcase ='case3012wp'; 

%case30 max_lam: 2.9859

%testcase='case30';lam=1.0; %L2ori=0.03302    L2opt=0.025954   
%testcase='case30';lam=1.3; %L2ori=0.058722   L2opt=0.045342
%testcase='case30';lam=1.6; %L2ori=0.094293   L2opt=0.07094
%testcase='case30';lam=1.9; %L2ori=0.14598    L2opt=diverge 
%testcase='case30';lam=2.0; %L2ori=0.1674    L2opt=diverge
%testcase='case30';lam=2.1; %L2ori=0.19153    L2opt=diverge
testcase='case30';lam=2.2; %L2ori=0.53642(diverge)    L2opt=diverge
%testcase='case30';lam=2.3; %L2ori=1.5911(diverge)    L2opt=diverge


%case118  max_lam: 1.4581 

%testcase='case118'; lam=1.0;  %L2ori=0.041209     L2opt=0.034756  
%testcase='case118'; lam=1.1;  %L2ori=0.050614     L2opt=0.041717   
%testcase='case118'; lam=1.2;  %L2ori=0.061059     L2opt=0.049762  
%testcase='case118'; lam=1.3;  %L2ori=0.072759     L2opt=0.058454
%testcase='case118'; lam=1.4;  %L2ori=0.086282       L2opt=0.068305  
%testcase='case118'; lam=1.5;  %L2ori=0.10337     L2opt=0.079183 
%testcase='case118'; lam=1.6;  %L2ori=0.12526     L2opt=0.091158
%testcase='case118'; lam=1.7;  %L2ori=0.17735(diverge)     L2opt=0.10387
%testcase='case118'; lam=1.8;  %L2ori=72764(diverge)     L2opt=diverge


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
results_0 = runpf(mpc,mpopt);

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
Lindex2_0_sum = sum(Lj_0.*conj(Lj_0));
% ====================================================================
% Initialization
% ====================================================================
t0 = clock;                % start timer
% Convert external to internal indexing
        % MATPOWER Fun: ext2int
mpc = ext2int(results_0);
% Distribute inputs to outputs
[baseMVA, bus, gen, branch] =...
             deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);
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
f_fcn = @(x) ObjFcn_sumL2(x);

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
Vml =0.5*bus(:, VMIN);                % lower bounds of Vm

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
% ====================================================================
gh_fcn = @(x) NlinFcn_sumL2(x,bus,gen);

% ====================================================================
% Hessian matrix
% ====================================================================  
cost_mult = 1;   % For interior-point algorithm
hess_fcn = @(x,lambda) HessFcn_sumL2(x,lambda,cost_mult);

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
% Calculate OPF L-index 
% ====================================================================

Vgen = V(idx_G);
Vload = V(idx_L);
Lj =ones(nb-ng,1)-(F_LG*Vgen)./Vload;
Lindex = abs(Lj);
Lindex2_sum =sum(Lj.*conj(Lj)) ;

% ====================================================================
% Print out the results of L-index and Max L-index
% ====================================================================

% Tabulate the results of each bus L-index
table(Lindex_0,Lindex,...
    'VariableNames',{'Lindex_Original' 'Lindex_Optimal'})   

% Tabulate the results of max L-index
table(Lindex2_0_sum,Lindex2_sum,...
    'VariableNames',{'Lindex2_sum_Original' 'Lindex2_sum_Optimal'}) 