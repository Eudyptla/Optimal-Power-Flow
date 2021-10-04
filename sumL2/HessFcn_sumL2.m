function Lxx = HessFcn_sumL2(x,lambda,cost_mult)

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

% Update voltage V 

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



% ====================================================================
% Evaluate gradient vector
% Gradient is the vector of first derivatives of the objective function
% ====================================================================

Vga = exp(1j*x(idx_G));
Vlm = x(idx_L+nb);
Z =diag(Vload)\F_LG*diag(Vgen);
Zc = conj(Z);
Zf = diag(Vload)\F_LG*diag(Vga);
Zfc= conj(Zf);
%Z = Fgl*Vgm*e^(1j*Vga)/Vlm*e^(1j*Vla);
L = F_LG*Vgen./Vload;
Lc = conj(L);
%sum(Fgl*Vgen)/Vload
% ====================================================================
% Evaluate Hessian matrix
% Hessian is the matrix of second derivatives of the objective function
% ====================================================================
d2f=sparse(nvars,nvars);

d2Z_dVga2 = -Z; %(15)
d2Zc_dVga2 = -Zc; %(16)
d2ZZc_dVgaa=diag((Z.*(Zc*(ones(ng,ng)-eye(ng)))).'*ones(nb-ng,1));
d2ZZc_dVgab2 =(Z'*Z+conj(Z'*Z)).*(ones(ng,ng)-eye(ng))-d2ZZc_dVgaa-conj(d2ZZc_dVgaa);
d2f(idx_G,idx_G) = diag((-d2Z_dVga2-d2Zc_dVga2).'*ones(nb-ng,1))+d2ZZc_dVgab2;

d2Z_dVgadVgm = 1j*diag(Vload)\F_LG*diag(Vga); %(17)
d2Zc_dVgadVgm = conj(d2Z_dVgadVgm); %(18)
d2ZZC_dVgaadVgm =diag((1j*Zf.*(Zc*(ones(ng,ng)-eye(ng)))).'*ones(nb-ng,1));
d2ZZc_dVgadVgm = (-1j*(Zc')*Zf+conj(-1j*(Zc')*Zf)).*(ones(ng,ng)-eye(ng))+conj(d2ZZC_dVgaadVgm)+d2ZZC_dVgaadVgm;
d2f(idx_G,idx_G+nb) =diag((-d2Z_dVgadVgm-d2Zc_dVgadVgm).'*ones(nb-ng,1))+d2ZZc_dVgadVgm;

d2Z_dVgadVla = Z; %(19)
d2Zc_dVgadVla = Zc; %(20)
d2f(idx_G,idx_L) = (-d2Z_dVgadVla-d2Zc_dVgadVla).' ;

d2Z_dVgadVlm = -1j*(diag(Vlm)\Z); %(21)
d2Zc_dVgadVlm = conj(d2Z_dVgadVlm);%(22)
d2ZZc_dVgadVlm = -2j*(diag(Vlm)\Z.*(Zc*(ones(ng,ng)-eye(ng))));
d2ZZc_dVgadVlm = d2ZZc_dVgadVlm+conj(d2ZZc_dVgadVlm);
d2f(idx_G,idx_L+nb) =(-d2Z_dVgadVlm-d2Zc_dVgadVlm+d2ZZc_dVgadVlm).';

%d2Z_dVgm2 = d2Zc_dVgm2 = 0 % (23) %(24)
d2ZZc_dVgm2 = diag((2*Zf.*conj(Zf))'*ones(nb-ng,1))+(Zf'*Zf+conj(Zf'*Zf)).*(ones(ng,ng)-eye(ng));
d2f(idx_G+nb,idx_G+nb) =d2ZZc_dVgm2 ;%ok

d2Z_dVgmdVla =-1j*(diag(Vload)\F_LG*diag(Vga)); %(25)
d2Zc_dVgmdVla = conj(d2Z_dVgmdVla);%(26)
d2f(idx_G+nb,idx_L) = (-d2Z_dVgmdVla-d2Zc_dVgmdVla).';

d2Z_dVgmdVlm = -diag(Vlm)\diag(Vload)\F_LG*diag(Vga); %(27)
d2Zc_dVgmdVlm =conj(d2Z_dVgmdVlm); %(28)
d2ZZC_dVgmVlm = -2*(diag(Vlm)\Zf.*(Zc*(ones(ng,ng)-eye(ng))));
d2ZZC_dVgmVlm = d2ZZC_dVgmVlm+conj(d2ZZC_dVgmVlm);
d2f(idx_G+ng,idx_L+nb) = (-d2Z_dVgmdVlm -d2Zc_dVgmdVlm+d2ZZC_dVgmVlm).';

d2Z_dVla2 = -L; %(29)
d2Zc_dVla2 = -Lc;%(30)
d2f(idx_L,idx_L) = diag(-d2Z_dVla2-d2Zc_dVla2);

d2Z_dVlaVlm = 1j*L./Vlm; %(31)
d2Zc_dVlaVlm = -1j*Lc./Vlm; %(32)
d2f(idx_L,idx_L+nb) = diag(-d2Z_dVlaVlm-d2Zc_dVlaVlm);

d2Z_dVlm2 = 2*L./Vlm./Vlm; %(33)
d2Zc_dVlm2 = 2*Lc./Vlm./Vlm; %(34)
d2ZZc_dVlm2 = 6*(L.*Lc./Vlm./Vlm);
d2f(idx_L+nb,idx_L+nb) = diag(-d2Z_dVlm2-d2Zc_dVlm2+d2ZZc_dVlm2);

d2f(1:nb+nb,1:nb)=d2f(1:nb,1:nb+nb).';

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
