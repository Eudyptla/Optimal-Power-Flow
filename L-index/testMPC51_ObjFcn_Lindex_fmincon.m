clear all; close all;

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA branch idx_G idx_L

% ====================================================================
% Input the name of the test case
% ====================================================================

%testcase = 'case5Lab';lam =1.0; %Liori= 0.075504   Liopt= 0.067357  
%testcase = 'case5Lab';lam =1.1; %Liori= diverge   Liopt= 0.067357 
%testcase = 'case5Lab';lam =1.4; %Liori= diverge   Liopt= 0.1045 
%testcase = 'case5Lab';lam =2.6; %Liori= diverge   Liopt= diverge

%testcase = 'case9Lab'; lam =1.0; %Liori= 0.054293   Liopt=  0.04699  
%testcase = 'case9Lab'; lam =1.3; %Liori= 0.072316   Liopt=  0.063147 
%testcase = 'case9Lab'; lam =1.4; %Liori= diverge   Liopt=  0.078189 
%testcase = 'case9Lab'; lam =1.5; %Liori= diverge   Liopt=diverge

%testcase = 'case26Lab';lam =1.0; %Liori= diverge  Liopt= 0.10943
%testcase = 'case26Lab';lam =1.1; %Liori= diverge  Liopt= diverge

%testcase ='case6ww'; lam =1.0; %Liori= 0.094362  Liopt= 0.09447 % Error
%%發電機輸出電壓不可調度
%testcase ='case6ww'; lam =1.3; %Liori= 0.13941  Liopt= 0.13083 
%testcase ='case6ww'; lam =1.4; %Liori= diverge  Liopt= 0.15061 
%testcase ='case6ww'; lam =1.5; %Liori= diverge  Liopt= diverge 

%testcase = 'case9'; lam =1.0; %Liori= 0.167  Liopt= 0.133
%testcase = 'case9'; lam =1.3; %Liori= 0.23316  Liopt= 0.18047
%testcase = 'case9'; lam =1.6; %Liori= 0.31301  Liopt= 0.23392
%testcase = 'case9'; lam =1.9; %Liori= 0.41368  Liopt= 0.29562
%testcase = 'case9'; lam =2.2; %Liori= 0.53741  Liopt= 0.37013
%testcase = 'case9'; lam =2.0; %Liori= 0.4202   Liopt= 0.31869
%testcase = 'case9'; lam =2.6; %Liori= diverge   Liopt= diverge
%testcase = 'case9'; lam =2.5; %Liori= 0.4202   Liopt= 0.31869

%testcase = 'case14';lam =1.0; %Liori= 0.076707  Liopt= 0.077219 
%初值界外 只有兩台發電機工作，三台不工作。
%testcase = 'case14';lam =1.3; %Liori= diverge  Liopt= 0.10311
%testcase = 'case14';lam =1.6; %Liori= diverge  Liopt= 0.1315
%testcase = 'case14';lam =1.9; %Liori= diverge  Liopt= 0.16837
%testcase = 'case14';lam =2.0; %Liori= diverge  Liopt= diverge
%testcase='case30';lam=1.0; %Liori=0.055286   Liopt=0.049124  
%testcase='case30';lam=1.3; %Liori=0.074035   Liopt=0.064563
testcase='case30';lam=1.75; %Liori=0.094943   Liopt=0.084815
%testcase='case30';lam=2.; %Liori=0.10292    Liopt=0.088856 
%原因不明
%testcase='case30';lam=2.2; %Liori=diverge    Liopt=diverge

%testcase = 'case39';lam =1.0; %Liori= 0.20098  Liopt= 0.19331
%testcase = 'case39';lam =1.1; %Liori= 0.22767  Liopt= 0.2141
%testcase = 'case39';lam =1.15; %Liori= diverge  Liopt= diverge

%testcase='case57'; lam =1.0; %Liori= 0.3099  Liopt= 0.2991
%testcase='case57'; lam =1.1; %Liori= 0.35736  Liopt= 0.3391
%testcase='case57'; lam =1.2; %Liori= 0.41373  Liopt= 0.38187
%testcase='case57'; lam =1.3; %Liori= 0.482585  Liopt= 0.42893
%testcase='case57'; lam =1.4; %Liori= 0.58586  Liopt= 0.48064
%testcase='case57'; lam =1.5; %Liori= diverge  Liopt= 0.54646
%testcase='case57'; lam =1.5; %Liori= diverge  Liopt= diverge

%testcase='case118'; lam=1.0;  %Liori=0.069389     Liopt=0.061445  
%testcase='case118'; lam=1.3;  %Liori=0.09327     Liopt=0.08195
%testcase='case118'; lam=1.6;  %Liori=0.12667     Liopt=0.10426  
%testcase='case118'; lam=1.7;  %Liori=diverge     Liopt=0.11194 
%testcase='case118'; lam=1.8;  %Liori=diverge     Liopt=diverge 

%testcase='case300'; lam=1.0;  %Liori=0.4202     Liopt=0.38506 
%testcase='case300'; lam=1.1;  %Liori=diverge     Liopt=diverge

%testcase='case89pegase'; lam=1.0;  %Liori=diverge     Liopt=diverge
%give up
%testcase='case1354pegase'; lam=1.0;  %Liori=diverge     Liopt=diverge
%testcase ='case2383wp';  
%testcase ='case3012wp'; 
%testcase='case24_ieee_rts'; lam=1.0;  %Liori=diverge     Liopt=diverge
% ====================================================================
% Load data from MATPOWER format
% ====================================================================
define_constants;            % MATPOWER use only
mpc = loadcase( testcase );  % MATPOWER Fun: loadcase
mpc = ext2int(mpc);
mpopt=mpoption;

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
T0 = Lindex_0_max;                       % Max L_index guess         

x0 = [Va;Vm;Pg;Qg;T0];             % column vector of initial points

% ====================================================================
% Scalar objective function
% ====================================================================
f_fcn = @(x) ObjFcn_Lindex(x);

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
QGmax =lam*gen(:, QMAX) / baseMVA;    % upper bounds of Qg
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
gh_fcn = @(x) NlinFcn_Lindex(x,bus,gen);

% ====================================================================
% Hessian matrix
% ====================================================================  
% cost_mult = 1;   % For interior-point algorithm
% hess_fcn = @(x,lambda) HessFcn_Ploss(x,lambda,cost_mult,Yf(il,:), Yt(il,:),bus);

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
                       'TolCon',5e-6,'TolX',1e-4,'TolFun',1e-4,...
                       'AlwaysHonorConstraints','none',...
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
Lindex_max = max(Lindex);

% ====================================================================
% Print out the results of L-index and Max L-index
% ====================================================================

% Tabulate the results of each bus L-index
table(Lindex_0,Lindex,...
    'VariableNames',{'Lindex_Original' 'Lindex_Optimal'})   

% Tabulate the results of max L-index
table(Lindex_0_max,Lindex_max,...
    'VariableNames',{'Lindex_MAX_Original' 'Lindex_MAX_Optimal'}) 