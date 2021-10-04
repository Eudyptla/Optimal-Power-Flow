 function [h,g,dh,dg] = NlinFcn_Multiobject(x,bus,gen)

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA branch idx_G idx_L

% ====================================================================
% Initialization
% ====================================================================
define_constants;          % MATPOWER use only

% Problem dimensions
nb = numAll(1);            % number of buses
nl = numAll(2);            % number of branches
ng = numAll(3);            % number of dispatchable injections
nvars = numAll(4);         % number of variables
nxtra = numAll(5);         % number of variables excluding Va and Vm

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
% Evaluate inequality constraints (L_index<=L_index_max) 
% ==================================================================== 
Y_LL = Ybus(idx_L,idx_L);
Y_LG = Ybus(idx_L,idx_G);
F_LG = (-1)*( Y_LL\Y_LG );
Vgen = V(idx_G);
Vload = V(idx_L);
Lj =ones(nb-ng,1)-(F_LG*Vgen)./Vload;
Lindex = abs(Lj);
h = Lindex-x(nvars)*ones(nb-ng,1);


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
% ====================================================================
% Evaluate partials of equality constraints (Gradient of h)
% ==================================================================== 
dh = sparse(nvars,nb-ng);
dh_den = Lindex;
Vga = exp(1j*x(idx_G));
Vlm = x(idx_L+nb);
Z =diag(Vload)\F_LG*diag(Vgen);
%Z = Fgl*Vgm*e^(1j*Vga)/Vlm*e^(1j*Vla);
L = F_LG*Vgen./Vload;
Lc = conj(L);
%sum(Fgl*Vgen)/Vload
dZ_dVga = diag(dh_den)\(0.5*1j*Z); %(7)
dZc_dVga = conj(dZ_dVga); %(8)
dZZc_dVga = diag(Lc)*dZ_dVga+diag(L)*dZc_dVga;
dh(idx_G,:) = (-dZ_dVga-dZc_dVga+dZZc_dVga).';

dZ_dVgm =0.5*( diag(dh_den)\(diag(Vload)\F_LG*diag(Vga))); %(9)
dZc_dVgm = conj(dZ_dVgm); %(10)
dZZc_dVgm = diag(Lc)*dZ_dVgm+diag(L)*dZc_dVgm;
dh(idx_G+nb,:) = (-dZ_dVgm-dZc_dVgm+dZZc_dVgm).';

dZ_dVla =0.5*( diag(dh_den)\((-1j)*Z)); %(11)
dZc_dVla = conj(dZ_dVla); %(12)
dh(idx_L,:) = diag((-dZ_dVla-dZc_dVla)*ones(ng,1));

dZ_dVlm =0.5*(diag(dh_den)\(-diag(Vlm)\Z)); %(13)
dZc_dVlm = conj(dZ_dVlm); %(14)
dZZc_dVlm = diag(Lc)*dZ_dVlm+diag(L)*dZc_dVlm;
dh(idx_L+nb,:) = diag((-dZ_dVlm-dZc_dVlm+dZZc_dVlm)*ones(ng,1));
dh(nvars,:)= -ones(1,nb-ng);

