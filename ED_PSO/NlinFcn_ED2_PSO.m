function [h,g] = NlinFcn_ED2_PSO(x,bus,gen)
h=[];

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA 

% ====================================================================
% Initialization
% ====================================================================
define_constants;           % MATPOWER use only

% Problem dimensions
nb = numAll(1);            % number of buses
nl = numAll(2);            % number of branches
ng = numAll(3);            % number of dispatchable injections
nvars = numAll(4);         % number of variables
nxtra = numAll(5);         

% Index of Va, Vm, Pg and Qg
iVa = 1:nb;                % index of Va
iVm = nb+1:2*nb;           % index of Vm
iPg = 2*nb+1:2*nb+ng;      % index of Pg
iQg = 2*nb+ng+1:2*nb+2*ng; % index of Qg

% ====================================================================
% Update x = [Va;Vm;Pg;Qg] during iteration
% ====================================================================
% Optimization Variables
Va = x(iVa)';                     % in rad
Vm = x(iVm)';                     % in p.u.
Pg = x(iPg)';                     % in p.u.
Qg = x(iQg)';                     % in p.u.

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

