clear all; close all;

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA branch

% ====================================================================
% Input the name of the test case
% ====================================================================
%testcase = 'case5Lab';   %Ploss=3.05
%testcase = 'case9Lab';   %Ploss=3.2575
%testcase = 'case26Lab';  %Ploss=10.794 

%testcase ='case6ww';lam=1.0;  %Ploss=6.67777
% testcase ='case6ww';lam=1.4;  %Ploss=13.923

%testcase = 'case9';lam=1.0;   %Ploss=2.3158
testcase = 'case9';lam=2.0;   %Ploss=12.558

%testcase='case30'; lam=1.0;   %Ploss=1.6036
%testcase='case30'; lam=1.7;   %Ploss=6.4855 

%testcase = 'case39';lam=1.0;  %Ploss=29.821
%testcase = 'case39';lam=1.15;  %Ploss=42.571
%testcase='case57';    %Ploss=11.306
%testcase='case118';   %Ploss=9.2405
%testcase='case300';   %Ploss=210.61

%testcase='case1354pegase'; %2023¬í Ploss =1000.8 48.52sec

%testcase ='case2383wp';  %Ploss=432.83 118.72sec
%testcase ='case3012wp';  %Ploss=472.54 241.49sec

% ====================================================================
% Load data from MATPOWER format
% ====================================================================
define_constants;            % MATPOWER use only
mpc = loadcase( testcase );  % MATPOWER Fun: loadcase

% ====================================================================
% Original power flow solution obtained from MATPOWER
% ====================================================================
% Set and retrieve a MATPOWER options struct
mpopt = mpoption('verbose', 2, 'out.all', 0,'pf.enforce_q_lims',1);   % MATPOWER Fun: mpoption

% ====================================================================
% Global variables
% ====================================================================
mpc = ext2int(mpc);
nb = size(mpc.bus, 1);         % number of buses
nl = size(mpc.branch, 1);      % number of branches
ng = size(mpc.gen, 1);         % number of generators
nvars = 2*nb + 2*ng;       % number of variables
nxtra = nvars - 2*nb;      % number of variables excluding Va and Vm

numAll = [nb;nl;ng;nvars;nxtra];

% Build admittance matrix
[Ybus,Yf,Yt] = makeYbus(mpc.baseMVA,mpc.bus,mpc.branch); % MATPOWER Fun: makeYbus
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
result0 = runpf(mpc, mpopt);                    % MATPOWER Fun: runpf
% Ploss cost (Original)
Plosses_0 = sum(real(get_losses(result0)));% MATPOWER Fun: get_losses

% ====================================================================
% Initialization
% ====================================================================
t0 = clock;                % start timer

% Convert external to internal indexing
mpc = ext2int(mpc);        % MATPOWER Fun: ext2int

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
f_fcn = @(x) ObjFcn_Ploss(x,bus);

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

% xmax = [Vau',Vmu',PGmax',QGmax'];    % upper bounds of Variables
% xmin = [Val',Vml',PGmin',QGmin'];    % lower bounds of Variables

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
% Find branches with flow limits (Eq.6.7, Eq.6.8)
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
nl2 = length(il);           % number of constrained lines
gh_fcn = @(x) NlinFcn_Ploss(x,bus,gen);

% ====================================================================
% Hessian matrix
% ====================================================================  
cost_mult = 1;   % For interior-point algorithm
hess_fcn = @(x,lambda) HessFcn_Ploss(x,lambda,cost_mult,bus);

% ====================================================================
% Set nondefault options
% ====================================================================
  
  options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
                       'AlwaysHonorConstraints','none',...
                       'GradObj','on',...      % Gradient of fun
                       'GradConstr','on',...   % Gradients of nonlcon
                       'Hessian','user-supplied','HessFcn',hess_fcn,... % Hessian of fun and nonlcon;
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
% Print out the results of Pg and Cost
% ====================================================================

Plosses_1 = fval;               % Pg cost (Economic)

% Tabulate the results of total real power loss
table(Plosses_0,Plosses_1,...
    'VariableNames',{'Plosses_Original' 'Plosses_Optimal'}) % in $/h