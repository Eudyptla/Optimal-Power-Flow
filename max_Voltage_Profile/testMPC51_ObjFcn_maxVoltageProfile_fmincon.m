clear all; close all;

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA branch idx_L idx_G

% ====================================================================
% Input the name of the test case
% ====================================================================
%testcase ='case6ww'; lam =1.0; %VD_0=0.014555 VD=0.014481
%%發電機輸出電壓不可調度
%testcase ='case6ww'; lam =1.3; %VD_0=0.076254 VD=0.046394

%testcase = 'case9'; lam =1.0; %VD_0=0.042379 VD=0.013418 
testcase = 'case9'; lam =2.0; %VD_0=0.1809 VD=0.040018

%testcase = 'case14';lam =1.0; %VD_0 =0.061945 VD=0.011309 
%初值界外 只有兩台發電機工作，三台不工作。 
%testcase = 'case14';lam =1.5; %VD_0= 0.046117  VD=0.017895

%testcase='case30';lam=1.0; %VD_0=0.039376   VD=0.011207  
%testcase='case30';lam=1.7; %VD_0=0.0844    VD=0.024445

%testcase = 'case39';lam =1.0;%VD_0=0.057896   VD= 0.028279 
%testcase = 'case39';lam =1.15;%VD_0=0.057896   VD= diverge 

%testcase='case57'; lam =1.0;%VD_0=0.064068   VD=0.053708 
%testcase='case57'; lam =1.4;%VD_0= 0.21756   VD=0.12674

%testcase='case118'; lam=1.0;%VD_0=0.054018   VD=Diverge 
%testcase='case118'; lam=1.5;%VD_0=0.10375   VD=Diverge 

%testcase='case300'; lam=1.0; %VD_0=0.079355   VD=0.054318


%testcase='case1354pegase';lam=1.0; %VD_0=0.094517 VD=0.063512
%testcase ='case2383wp';  lam=1.0;VD_0=diverge VD=diverge
%give up
%testcase ='case3012wp'; 

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
% Calculate Original sum Voltage Profile 
% ====================================================================
Vm_0 = results_0.bus(:,VM);               % in p.u.
Voltage_Profile_0=Vm_0(idx_L)-ones(nb-ng,1);
max_Voltage_Profile_0= max(abs(Voltage_Profile_0));

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
T0 = max_Voltage_Profile_0;
x0 = [Va;Vm;Pg;Qg;T0];             % column vector of initial points

% ====================================================================
% Scalar objective function
% ====================================================================
f_fcn = @(x) ObjFcn_max_Voltage_Profile(x);

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
Tmax=99;
Tmin=0;

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
gh_fcn = @(x) NlinFcn_max_Voltage_Profile(x,bus,gen);

% ====================================================================
% Hessian matrix
% ====================================================================  
cost_mult = 1;   % For interior-point algorithm
hess_fcn = @(x,lambda) HessFcn_maxVoltage(x,lambda,cost_mult);

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
% Print out the results of Voltage Profile
% ====================================================================


% Tabulate the results of max L-index
table(max_Voltage_Profile_0,fval,...
    'VariableNames',{'Voltage_Profile_Original' 'Voltage_Profile_Optimal'}) 