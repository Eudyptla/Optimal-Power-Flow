clear all; close all;

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA branch idx_G idx_L gencost 

% ====================================================================
% Input the name of the test case
% ====================================================================
% testcase ='case6ww'; lam =1.0; %%發電機輸出電壓不可調度
% goal = [3144,0.09447]; weight =[0.01,0.1].*goal;

% Original    Optimal      Goal  
%     ________    ________    _______
% 
%       3189.5      3143.9       3144
%     0.094362    0.094431    0.09447

% testcase ='case6ww'; lam =1.4; %%發電機輸出電壓不可調度
% goal = [4137.9,0.15061]; weight =[0.01,0.1].*goal;
% Original    Optimal     Goal  
%     ________    _______    _______
% 
%      4364.6      4262.9     4137.9
%     0.16252     0.15221    0.15061
%Original diverge
 
% testcase = 'case9'; lam =1.0;
% goal = [5296.7,0.13327]; weight =[0.01,0.1].*goal;
%  Original    Optimal     Goal  
%     ________    _______    _______
% 
%     5438.3       5299.1     5296.7
%     0.1673      0.13388    0.13327

% testcase = 'case9'; lam =2.4;
% goal = [23896,0.43224]; weight =[0.01,0.1].*goal;
% Original    Optimal     Goal  
%     ________    _______    _______
% 
%       24707      23901       23896
%     0.65899     0.4331     0.43224


% testcase='case30';lam=1.0; 
% goal=[576.89,0.049124];weight =[0.01,0.1].*goal;
% 
% Original    Optimal      Goal  
%     ________    _______    ________
% 
%       593.45    576.85       576.89
%     0.055286     0.049     0.049124

% testcase='case30';lam=1.7; 
% goal=[1178.1,0.088856];weight =[0.01,0.1].*goal;

% Original    Optimal       Goal  
%     ________    ________    ________
% 
%      1200.9       1178.2      1178.1
%     0.10292     0.089077    0.088856

% testcase = 'case39';lam =1.0; 
% goal = [41864,0.19331] ; weight = [0.01,0.1].*goal;
%  
% Original    Optimal     Goal  
%     ________    _______    _______
% 
%       45077       41867      41864
%     0.20098     0.19343    0.19331

% testcase = 'case39';lam =1.15; 
% goal = [56643,0.22681] ; weight = [0.01,0.1].*goal;
% 
% Original    Optimal     Goal  
%     ________    _______    _______
% 
%       61451       57501      56643
%     0.28403     0.25351    0.22681


% testcase='case57'; lam =1.0; 
% goal = [41738,0.2991];weight = [0.01,0.1].*goal;

%  Original    Optimal     Goal 
%     ________    _______    ______
% 
%      51348        41743     41738
%     0.3099      0.29945    0.2991


% testcase='case300'; lam=1.0;  
% goal=[719730,0.38506];weight = [0.01,0.1].*goal;

%  Original      Optimal         Goal   
%     __________    __________    __________
% 
%     7.2475e+05    7.1992e+05    7.1973e+05
%         0.4202        0.3861       0.38506

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
nvars = 2*nb + 2*ng+1;       % number of variables
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
% ====================================================================
% Calculate Original ED
% ====================================================================
totalcost =sum(totcost(results_0.gencost, results_0.gen(:,PG)));
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
% Original
% ====================================================================
 Original_value=[totalcost,Lindex_0_max]; 
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
T0 = Lindex_0_max;                       % Max L_index guess         

x0 = [Va;Vm;Pg;Qg;T0];             % column vector of initial points

% ====================================================================
% Scalar objective function
% ====================================================================
f_fcn = @(x) ObjFcn_ED_Lindex(x);

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

%Lindex bound
Tmax = 1;
Tmin = 0;
% Variable limits (Eq.6.4)
xmax = [Vau;Vmu;PGmax;QGmax;Tmax];      % upper bounds of variables
xmin = [Val;Vml;PGmin;QGmin;Tmin];      % lower bounds of variables

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
gh_fcn = @(x) NlinFcn_ED_Lindex(x,bus,gen);

% ====================================================================
% Set nondefault options
% ====================================================================

options = optimoptions('fgoalattain','Display','iter',...
                       'GradObj','on',...      % Gradient of fun
                       'GradConstr','on',...   % Gradients of nonlcon
                       'TolCon',5e-6,'TolX',1e-4,'TolFun',1e-4,...
                       'PlotFcn',@optimplotstepsize);
                 

% ====================================================================
% Find minimum of constrained nonlinear multivariable function
% ====================================================================
[x,fval,attainfactor,exitflag,Output,lambda] = fgoalattain(f_fcn,x0,goal,weight, A, b, Aeq, beq,...
                                          xmin, xmax, gh_fcn, options);

% ====================================================================
% Summarize the results
% ====================================================================
success = (exitflag > 0);

%Index of Va, Vm, Pg and Qg
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
% Print out the results of L-index and Max L-index
% ====================================================================

table(Original_value',fval,goal',...
    'VariableNames',{'Original' 'Optimal' 'Goal'})   

attainfactor

