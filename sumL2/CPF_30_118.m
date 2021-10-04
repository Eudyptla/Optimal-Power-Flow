clear all; close all;
% ====================================================================
% Input the name of the test case
% ====================================================================
%testcase = 'case9';  %max_lam: 0.9876
%testcase='case30';   %max_lam: 2.9859
%testcase = 'case39';  %max_lam: 0.7571
testcase='case57';    %max_lam: 0.5945

%testcase='case118';   %max_lam: 1.4581

% ====================================================================
% Input the name of the test case
% ====================================================================
R =2.5;
define_constants;
mpopt= mpoption('out.all',0,'verbose',2,'cpf.stop_at','FULL','cpf.step',0.2);
mpopt= mpoption(mpopt,'cpf.plot.level',2);
mpcb=loadcase(testcase);
mpct=mpcb;
mpct.gen(:,[PG QG])=mpcb.gen(:,[PG QG])*R;
mpct.bus(:,[PD QD])=mpcb.bus(:,[PD QD])*R;
results = runcpf(mpcb,mpct,mpopt);
results.cpf
