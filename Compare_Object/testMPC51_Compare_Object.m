clear all; close all;
% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus Yf Yt baseMVA  gencost idx_G idx_L

% ====================================================================
% Input the name of the test case
% ====================================================================
% testcase ='case6ww'; lam =1.0; 

%testcase ='case6ww'; lam =1.4; 

%testcase = 'case9'; lam =1.0;

%testcase = 'case9'; lam =2;

%testcase = 'case14';lam =1.0; 

%testcase = 'case14';lam =1.5;

%testcase='case30';lam=1.0; 

%testcase='case30';lam=1.7; 

%testcase = 'case39';lam =1.0;

%testcase = 'case39';lam =1.1;

%testcase='case57'; lam =1.0;

%testcase='case57'; lam =1.1;

%testcase='case118'; lam=1.0;
%ERROR 
%testcase='case118'; lam=1.5; 
%ERROR 

testcase='case300'; lam=1.0;

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
%mpopt=mpoption('pf.enforce_q_lims',1,'out.all',0);
mpopt=mpoption;
% ====================================================================
% Global variables
% ====================================================================
nb = size(mpc.bus, 1);         % number of buses
nl = size(mpc.branch, 1);      % number of branches
ng = size(mpc.gen, 1);         % number of generators
nvars = 2*nb + 2*ng;           % number of variables
nxtra = nvars - 2*nb;          % number of variables excluding Va and Vm
idx_G = find((mpc.bus(:,BUS_TYPE) == 3) | (mpc.bus(:,BUS_TYPE) == 2));    % index of slack &PV bus
idx_L = find(mpc.bus(:,BUS_TYPE) == 1 );                              % index of PQ bus
numAll = [nb;nl;ng];
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
[baseMVA, bus, gen, branch,gencost] =...
             deal(results_0.baseMVA, results_0.bus, results_0.gen, results_0.branch,results_0.gencost);


% ====================================================================
% Calculate Original Pcost
% ====================================================================
coeff2 = gencost(:,5);     % ($/MW^2)
coeff1 = gencost(:,6);     % ($/MW)
coeff0 = gencost(:,7);     % ($)
Pg_0 = results_0.gen(:,2);
Qg_0 = results_0.gen(:,3);
Pcost_0= sum( coeff2.*(Pg_0.^2)+coeff1.*Pg_0+coeff0);
% ====================================================================
% Calculate Original Ploss
% ====================================================================
Ploss_0 = sum(real(get_losses(results_0)));% MATPOWER Fun: get_losses
% ====================================================================
% Calculate Original Voltage Profile
% ====================================================================
Vm_0 = results_0.bus(:,VM);               % in p.u.
Voltage_Deviations_0=Vm_0(idx_L)-ones(nb-ng,1);
Voltage_Profile_0= sum(abs(Voltage_Deviations_0));
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
% Calculate Original Q Reserve 
% ====================================================================
QReserve_0= sum(results_0.gen(:,3));
% ====================================================================
% Calculate Original Qcost
% ====================================================================

PQ_0=Pg_0.*Qg_0;
S2_0=Pg_0.^2+Qg_0.^2;
Qcost_0= sum( coeff2.*((PQ_0./S2_0).^2).*(Qg_0.^2)+coeff1.*(PQ_0./S2_0).*Qg_0+coeff0);


Original_value = [Pcost_0,Ploss_0,Voltage_Profile_0,Lindex_0_max,QReserve_0...
                  Qcost_0];
% ====================================================================
% Calculate Other Object
% ====================================================================
Obj_ED = ED_fmincon(results_0);
Obj_Ploss = Ploss_fmincon(results_0);
Obj_VP = VP_fmincon(results_0);
Obj_Lindex = Lindex_fmincon(results_0,Lindex_0_max);
Obj_QR = QR_fmincon(results_0);
Ogj_TC = TC_fmincon(results_0);
% ====================================================================
% Print out comparison of soltions
% ====================================================================

A=[Original_value;Obj_ED;Obj_Ploss;Obj_VP;Obj_Lindex;Obj_QR;Ogj_TC];

% Tabulate the results 
T = array2table(A,...
    'VariableNames',{'PgCost','Ploss','VoltageProfile','Lindex','QReserve','QgCost'},...
    'RowNames',{'Original','ED','Ploss','VoltageProfile','Lindex','QReserve','TotalCost'})
