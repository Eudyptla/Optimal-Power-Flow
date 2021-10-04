clear all; close all;

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA branch idx_L idx_G

% ====================================================================
% Input the name of the test case
% ====================================================================
%testcase ='case6ww'; lam =1.0; %VD_0=0.27021 VD=0.028034 
%%發電機輸出電壓不可調度
testcase ='case6ww'; lam =1.4; %VD_0=0.29968 VD=0.18635

%testcase = 'case9'; lam =1.0; %VD_0=0.10145 VD=0.052376 
%testcase = 'case9'; lam =2.0; %VD_0=0.60183 VD=0.1624

%testcase = 'case14';lam =1.0; %VD_0 =0.4078 VD=0.034321
%初值界外 只有兩台發電機工作，三台不工作。 
%testcase = 'case14';lam =1.5; %VD_0= 0.24182  VD=0.064979

testcase='case30';lam=1.0; %VD_0=0.5417   VD=0.14464  
%testcase='case30';lam=1.7; %VD_0=1.1664    VD=0.26176

%testcase = 'case39';lam =1.0;%VD_0= 0.82815   VD= 0.21675 
%testcase = 'case39';lam =1.1;%VD_0= 0.7493   VD=0.2864
%奇怪的結果
%testcase='case57'; lam =1.0;%VD_0=1.2336   VD=1.0986
%testcase='case57'; lam =1.4;%VD_0=3.9074   VD=3.2498

%testcase='case118'; lam=1.0;%VD_0=1.4221   VD=0.4159
%testcase='case118'; lam=1.5;%VD_0=2.2946   VD=0.69452 

%testcase='case300'; lam=1.0; %VD_0=5.5411   VD=3.7106

%testcase='case1354pegase';lam =1.0; %VD_0=43.886 VD=11.139
%testcase ='case2383wp'; lam=1.0 ;%VD_0=diverge VD=461.87
%PF無解，但是OPF有解
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
% Scalar objective function
% ====================================================================
f_fcn = @(x) ObjFcn_sum_Voltage_Profile(x);

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
%  Nonlinear constraints: (Eq.6.2, Eq.6.3)
%  Nonlinear equality constraints:   power balance equations (Eq.4.2, Eq.4.3)
% ====================================================================
gh_fcn = @(x) NlinFcn_sum_Voltage_Profile(x,bus,gen);

% ====================================================================
% Hessian matrix
% ====================================================================  
cost_mult = 1;   % For interior-point algorithm
hess_fcn = @(x,lambda) HessFcn_VoltageProfile(x,lambda,cost_mult);

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
% Calculate Original sum Voltage Profile 
% ====================================================================
Vm_0 = results_0.bus(:,VM);               % in p.u.
Voltage_Profile_0=Vm_0(idx_L)-ones(nb-ng,1);
sum_Voltage_Profile_0= sum(abs(Voltage_Profile_0));

% ====================================================================
% Print out the results of Voltage Profile
% ====================================================================


% Tabulate the results of max L-index
table(sum_Voltage_Profile_0,fval,...
    'VariableNames',{'Voltage_Profile_Original' 'Voltage_Profile_Optimal'}) 