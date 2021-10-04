clear all; close all;

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA

% ====================================================================
% Input the name of the test case
% ====================================================================
%testcase = 'case5Lab';   %cost=1,603.7
%testcase = 'case9Lab';   %cost=3,188 
%testcase = 'case26Lab';  %cost=15,423

%testcase ='case6ww';  %cost=3,144
testcase = 'case9';   %cost=5,296.7
%testcase='case30';    %cost=576.89
%testcase = 'case39';  %cost=41,864
%testcase='case57';    %cost=41,738
%testcase='case300';   %cost=719,730

%testcase ='case2383wp';  %cost=1,868,500
%testcase ='case3012wp';  %cost=2,591,768

% ====================================================================
% Load data from MATPOWER format
% ====================================================================
define_constants;            % MATPOWER use only
mpc = loadcase( testcase );  % MATPOWER Fun: loadcase

% ====================================================================
% Original power flow solution obtained from MATPOWER
% ====================================================================
% Set and retrieve a MATPOWER options struct
mpopt = mpoption('verbose', 2, 'out.all', 0);   % MATPOWER Fun: mpoption
result0 = runpf(mpc, mpopt);                    % MATPOWER Fun: runpf

genPout0_ALL = result0.gen(:,PG);    % Pg (Original)
genPout0_SUM = sum( genPout0_ALL );  % Pd + Ploss (Original)

% Pg cost (Original)
genCost0_ALL = sum(totcost(mpc.gencost,genPout0_ALL));  % MATPOWER Fun: totcost

% ====================================================================
% Initialization
% ====================================================================
t0 = clock;                % start timer

% Convert external to internal indexing
mpc = ext2int(result0);        % MATPOWER Fun: ext2int

% Distribute inputs to outputs
[baseMVA,bus,gen,branch,gencost] =...
deal(mpc.baseMVA,mpc.bus,mpc.gen,mpc.branch,mpc.gencost);

% ====================================================================
% Set global variables
% ====================================================================
% No. of bus, branch, gen, variables
nb = size(bus, 1);         % number of buses
nl = size(branch, 1);      % number of branches
ng = size(gen, 1);         % number of dispatchable injections
nvars = 2*nb + 2*ng;       % number of variables
nxtra = nvars - 2*nb;

numAll = [nb;nl;ng;nvars;nxtra];

% Build admittance matrices
[Ybus,Yf,Yt] = makeYbus(baseMVA,bus,branch); % MATPOWER Fun: makeYbus

% ====================================================================
% Set up initial variables in p.u. values
% ====================================================================
Va = bus(:, VA) * (pi / 180);  % in rad
Vm = bus(:, VM);               % in p.u.
Pg = gen(:, PG) / baseMVA;     % in p.u.
Qg = gen(:, QG) / baseMVA;     % in p.u.

% x0 = [Va;Vm;Pg;Qg];          % initial point

% ====================================================================
% Scalar objective function
% fun should accept a row vector of length nvars and return a scalar value
% ====================================================================
fun = @(x) ObjFcn_ED2_PSO(x,gencost);

% ====================================================================
% Bound constraints (Eq.6.4)
% ====================================================================

% Voltage angle bounds for swing bus (Eq.6.10)
refs = find(bus(:, BUS_TYPE) == REF);
Vau = Inf(nb, 1);                % upper bounds of Va 
Val = -Vau;                      % lower bounds of Va 

Vau(refs) = Va(refs);            % upper bounds of Va_swing
Val(refs) = Va(refs);            % lower bounds of Va_swing

% Voltage magnitude bounds (Eq.6.11)
Vmu = bus(:,VMAX);               % upper bounds of Vm
Vml = bus(:,VMIN);               % lower bounds of Vm

% Gen real power bounds (Eq.6.12)
PGmin = gen(:, PMIN) / baseMVA;  % lower bounds of Pg
PGmax = gen(:, PMAX) / baseMVA;  % upper bounds of Pg

% Gen reactive power bounds (Eq.6.13)
QGmin = gen(:, QMIN) / baseMVA;  % lower bounds of Qg
QGmax = gen(:, QMAX) / baseMVA;  % upper bounds of Qg

%Variable limits (Eq.6.4)
UB = [Vau;Vmu;PGmax;QGmax];      % upper bounds of Variables
LB = [Val;Vml;PGmin;QGmin];      % lower bounds of Variables

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
gh_fcn = @(x) NlinFcn_ED2_PSO(x,bus,gen);

% ====================================================================
% Set Nondefault Options
% ====================================================================
% options = gaoptimset('Display','iter',...
%                      'PopulationSize',2000,...
%                      'Vectorized','on',...
%                      'PlotFcns',@gaplotbestf);

options = psooptimset('ConstrBoundary','penalize',...
'UseParallel','ture' ,'Display','iter','TolCon',1e-04,...
'TolFun',1e-04,'Vectorized','on','PopulationType','doubleVector',...
'HybridFcn',@fmincon);
options.PlotFcns = {@psoplotbestf,@psoplotswarmsurf} ;
                  
% ====================================================================
% Find minimum of function using genetic algorithm
% ====================================================================
[x,fval,exitflag,output] = pso(fun,nvars,A,b,Aeq,beq,LB,UB,...
                              gh_fcn,options);

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
Va = x(iVa)';                      % in rad
Vm = x(iVm)';                      % in p.u.
Pg = x(iPg)';                      % in p.u.
Qg = x(iQg)';                      % in p.u.

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
genPout1_ALL = gen(:, PG);         % Pg (Economic)
genCost1_ALL = fval;               % Pg cost (Economic)

% Tabulate the results of each generation
table(genPout0_ALL,genPout1_ALL,...
    'VariableNames',{'Pg_Original' 'Pg_Economic'})           % in MW

% Tabulate the results of total generation cost
table(genCost0_ALL,genCost1_ALL,...
    'VariableNames',{'Pg_Cost_Original' 'Pg_Cost_Economic'}) % in $/h