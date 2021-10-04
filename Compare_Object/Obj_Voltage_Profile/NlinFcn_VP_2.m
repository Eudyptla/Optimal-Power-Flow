function [h,g,dh,dg] =  NlinFcn_VP_2(x,bus,gen)
h=[];
dh=[];
% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus Yf Yt baseMVA  gencost idx_G idx_L

% ====================================================================
% Initialization
% ====================================================================
define_constants;          % MATPOWER use only

% Problem dimensions
nb = numAll(1);            % number of buses
nl = numAll(2);            % number of branches
ng = numAll(3);            % number of dispatchable injections
nvars = 2*nb + 2*ng;           % number of variables
nxtra = nvars - 2*nb;          % number of variables excluding Va and Vm

% Index of Va, Vm, Pg and Qg
iVa = 1:nb;                % index of Va
iVm = nb+1:2*nb;           % index of Vm
iPg = 2*nb+1:2*nb+ng;      % index of Pg
iQg = 2*nb+ng+1:2*nb+2*ng; % index of Qg

% ====================================================================
% Update optimization variables x = [Va;Vm;Pg;Qg] during iteration
% ====================================================================
Va = x(iVa);                       % in rad
Vm = x(iVm);                       % in p.u.
Pg = x(iPg);                       % in p.u.
Qg = x(iQg);                       % in p.u.

% Update some generator data
gen(:, PG) = Pg * baseMVA;         % in MW
gen(:, QG) = Qg * baseMVA;         % in MVAr

% Update voltage phasors
V = Vm .* exp(1j * Va);

% ====================================================================
% Evaluate equality constraint function (Eq.6.2)
% ====================================================================

% Build the vector of complex bus power injections (CgSg-Sd, Eq.3.17)
Sbus = makeSbus(baseMVA, bus, gen); % MATPOWER Fun: makeSbus

% Evaluate power balance equations
mis = V .* conj(Ybus * V) - Sbus;   % (Eq.3.17)

g = [ real(mis);       % active power mismatch for all buses   (Eq.4.2)
      imag(mis)];      % reactive power mismatch for all buses (Eq.4.3)
% ====================================================================
% Evaluate partials of equality constraints (Gradient of g)
% ==================================================================== 

% Computes partial derivatives of power injection w.r.t. voltage
[dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);           % w.r.t. V
neg_Cg = sparse(gen(:, GEN_BUS), 1:ng, -1, nb, ng);   % Pbus w.r.t. Pg

% construct Jacobian of equality (power flow) constraints and transpose it
dg = sparse(2*nb, nvars);
dg(:, [iVa iVm iPg iQg]) = [
   real([dSbus_dVa dSbus_dVm]) neg_Cg sparse(nb, ng);  % P mismatch w.r.t Va, Vm, Pg, Qg
   imag([dSbus_dVa dSbus_dVm]) sparse(nb, ng) neg_Cg;  % Q mismatch w.r.t Va, Vm, Pg, Qg
   ];

dg = dg';
