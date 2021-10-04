function [f,df] = ObjFcn_sumL2(x)

% f:   Objective function (scalar) as the first output
% df:  Gradient (vector) as the second output
% d2f: Hessian (matrix) as the third output

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA branch idx_G idx_L

% ====================================================================
% Initialization
% ====================================================================
define_constants;           % MATPOWER use only

% Problem dimensions
nb = numAll(1);            % number of buses
nl = numAll(2);            % number of branches
ng = numAll(3);            % number of generators
nvars = numAll(4);         % number of variables
nxtra = numAll(5);         % number of variables excluding Va and Vm


% Index of Va, Vm, Pg and Qg
iVa = 1:nb;                % index of Va
iVm = nb+1:2*nb;           % index of Vm


% ====================================================================
% Update optimization variables x = [Va;Vm;Pg;Qg] during iteration
% ====================================================================
Va = x(iVa);               % in rad
Vm = x(iVm);               % in p.u.
V = Vm .* exp(1j * Va);

% ====================================================================
% Evaluate objective function
% ====================================================================

Vgen = V(idx_G);
Vload = V(idx_L);
Y_LL = Ybus(idx_L,idx_L);
Y_LG = Ybus(idx_L,idx_G);
F_LG = (-1)*( Y_LL\Y_LG );
Lj =ones(nb-ng,1)-(F_LG*Vgen)./Vload;
f=sum(Lj.*conj(Lj));

% ====================================================================
% Evaluate gradient vector
% Gradient is the vector of first derivatives of the objective function
% ====================================================================

Vga = exp(1j*x(idx_G));
Vlm = x(idx_L+nb);
Z =diag(Vload)\F_LG*diag(Vgen);
Zc = conj(Z);
%Z = Fgl*Vgm*e^(1j*Vga)/Vlm*e^(1j*Vla);
L = F_LG*Vgen./Vload;
Lc = conj(L);
%sum(Fgl*Vgen)/Vload

df = sparse(nvars,1);

dZ_dVga = 1j*Z; %(7)
dZc_dVga = conj(dZ_dVga); %(8)
dZZc_dVga = diag(Lc)*dZ_dVga+diag(L)*dZc_dVga;
df(idx_G,:) = (-dZ_dVga-dZc_dVga+dZZc_dVga).'*ones(nb-ng,1);

dZ_dVgm = diag(Vload)\F_LG*diag(Vga); %(9)
dZc_dVgm = conj(dZ_dVgm); %(10)
dZZc_dVgm = diag(Lc)*dZ_dVgm+diag(L)*dZc_dVgm;
df(idx_G+nb,:) = (-dZ_dVgm-dZc_dVgm+dZZc_dVgm).'*ones(nb-ng,1);

dZ_dVla = -1j*Z; %(11)
dZc_dVla = conj(dZ_dVla); %(12)
df(idx_L,:) = (-dZ_dVla-dZc_dVla)*ones(ng,1);

dZ_dVlm =-diag(Vlm)\Z; %(13)
dZc_dVlm = conj(dZ_dVlm); %(14)
dZZc_dVlm = diag(Lc)*dZ_dVlm+diag(L)*dZc_dVlm;
df(idx_L+nb,:) = (-dZ_dVlm-dZc_dVlm+dZZc_dVlm)*ones(ng,1);



 