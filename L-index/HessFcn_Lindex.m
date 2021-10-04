function Lxx = HessFcn_Lindex(x,lambda,cost_mult)

% Hessian of objective
% Hessian of nonlinear inequality and nonlinear equality constraint 

% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus baseMVA  branch  idx_G idx_L

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

Va = x(iVa);                       % in rad
Vm = x(iVm);                       % in p.u.

% Update voltage phasors
V = Vm .* exp(1j * Va);
% ====================================================================
% Hessian of objective
% ====================================================================

Va = x(iVa);               % in rad
Vm = x(iVm);               % in p.u.
V = Vm .* exp(1j * Va);

Vgen = V(idx_G);
Vload = V(idx_L);
Y_LL = Ybus(idx_L,idx_L);
Y_LG = Ybus(idx_L,idx_G);
F_LG = (-1)*( Y_LL\Y_LG );
% Evaluate d2f

% construct complex bus voltage vector
% create map of external bus numbers to bus indices
% create map of external bus numbers to bus indices


Vga = exp(1j*x(idx_G));
Vgm = x(idx_G+nb);
Vla = exp(1j*x(idx_L));
Vlm = x(idx_L+nb);
Z =diag(Vload)\F_LG*diag(Vgen);
%Z = Fgl*Vgm*e^(1j*Vga)/Vlm*e^(1j*Vla);
ZZc_Vga = Z.*conj(Z*(ones(ng,ng)-eye(ng)));
L = F_LG*Vgen./Vload;
%sum(Fgl*Vgen)/Vload
dZ_dVga = -1j*Z+1j*conj(Z);
%dZ/dVga = 1j*Z 
%dZZc/dVga = 1j*ZZc
dZ_Vgm = diag(Vload)\F_LG*diag(Vga);
%dZ/dVgm = Fgl*e^(1j*Vga)/Vlm*e^(1j*Vla)
dZZc_Vgm = dZ_Vgm.*conj(Z*(ones(ng,ng)-eye(ng)));
%dZZc/sVgm = sum((Vgm2*e^(Vga1-Vga2))/Vlm^2)
dZ_dVla = 1j*L-1j*conj(L);
%dZ/dVla = 1j*L
dZ_Vlm =L./Vlm;
dZ_dVlm = dZ_Vlm+conj(dZ_Vlm);
%dZ/dVlm=L/Vlm;
dZZc_dVlm = -2*L.*conj(L)./Vlm;

d2f=sparse(nvars,nvars);
d2Z_Vga2 = Z+conj(Z);
d2Z_dVga2 =diag((d2Z_Vga2.')*ones(nb-ng,1));
d2ZZc_Vga2 = -ZZc_Vga-conj(ZZc_Vga);
d2ZZc_dVga2 = diag((d2ZZc_Vga2.')*ones(nb-ng,1));
d2ZZcd_Vg1aVg2a = (Z.')*conj(Z).*(ones(ng,ng)-eye(ng));
d2ZZcd_dVg1aVg2a= +d2ZZcd_Vg1aVg2a+conj(d2ZZcd_Vg1aVg2a);
d2f(idx_G,idx_G) = d2Z_dVga2+d2ZZc_dVga2+d2ZZcd_dVg1aVg2a;

d2Z_dVgaVla = -Z-conj(Z);


d2Z_dVla2 = L+conj(L);


d2Z_VgaVgm = (-1j*dZ_Vgm+1j*conj(dZ_Vgm)).';
d2Z_dVgadVgm = diag(d2Z_VgaVgm*ones(nb-ng,1));
d2ZZc_VgaVgm = (1j*dZZc_Vgm-1j*conj(dZZc_Vgm)).';
d2ZZc_dVgadVgm = diag(d2ZZc_VgaVgm*ones(nb-ng,1));
d2ZZc_Vga1Vgm2 =(dZ_Vgm'*Z).*(ones(ng,ng)-eye(ng));
d2ZZc_dVga1dVgm2 = 1j*d2ZZc_Vga1Vgm2-1j*conj(d2ZZc_Vga1Vgm2);




d2Z_dVgaVlm = -(diag(Vlm)\dZ_dVga);
d2ZZc_VgaVlm = diag(Vlm)\ZZc_Vga;
d2ZZc_dVgadVlm = -2j*d2ZZc_VgaVlm+2j*conj(d2ZZc_VgaVlm);


d2Z_dVlaVgm = 1j*dZ_Vgm-1j*conj(dZ_Vgm);
d2f(idx_L,idx_G+nb)=d2Z_dVlaVgm;
d2Z_dVlaVlm = diag(-dZ_dVla./Vlm);
d2f(idx_L,idx_L+nb) = d2Z_dVlaVlm;



d2ZZc_Vgm2 =  2*diag((dZ_Vgm.*conj(dZ_Vgm))'*ones(nb-ng,1));
d2ZZc_Vgm1Vgm2 = dZ_Vgm'*dZ_Vgm*(ones(ng,ng)-eye(ng));
d2ZZc_dVgm1dVgm2 = d2ZZc_Vgm1Vgm2+conj(d2ZZc_Vgm1Vgm2);


d2Z_dVlm2 = -2*dZ_dVlm./Vlm;
d2ZZc_dVlm2 = -3*dZZc_dVlm./Vlm;


d2Z_VgmVlm = (-diag(Vlm)\dZ_Vgm).';
d2Z_dVgmdVlm = -d2Z_VgmVlm-conj(d2Z_VgmVlm);

d2ZZc_VgmVlm = -2*(diag(Vlm)\dZZc_Vgm).';

d2ZZc_dVgmdVlm = d2ZZc_VgmVlm+conj(d2ZZc_VgmVlm)-2*((diag(Vlm)\(dZ_Vgm.*conj(Z)+conj(dZ_Vgm).*Z)).');

d2f = d2f * cost_mult;

% ====================================================================
% Hessian of nonlinear equality constraint
% ====================================================================


nlam = length(lambda.eqnonlin) / 2;
lamP = lambda.eqnonlin(1:nlam);
lamQ = lambda.eqnonlin((1:nlam)+nlam);
[Gpaa, Gpav, Gpva, Gpvv] = d2Sbus_dV2(Ybus, V, lamP);
[Gqaa, Gqav, Gqva, Gqvv] = d2Sbus_dV2(Ybus, V, lamQ);

% Evaluate Hessian of power balance constraints
d2G = [
    real([Gpaa Gpav; Gpva Gpvv]) + imag([Gqaa Gqav; Gqva Gqvv]) sparse(2*nb, nxtra);
    sparse(nxtra, 2*nb + nxtra)
];

% ====================================================================
% Hessian of objective, nonlinear inequality, and nonlinear equality constraints
% Hessian is the matrix of second derivatives of the Lagrangian
% ====================================================================
Lxx = d2f + d2G  ;
