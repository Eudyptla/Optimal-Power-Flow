function [f,df] = ObjFcn_Lindex_2(x)

% f:   Objective function (scalar) as the first output
% df:  Gradient (vector) as the second output


% ====================================================================
% Declare variables as global
% ====================================================================
global numAll Ybus Yf Yt baseMVA  gencost idx_G idx_L
% ====================================================================
% Initialization
% ====================================================================
% Problem dimensions
nb = numAll(1);            % number of buses
nl = numAll(2);            % number of branches
ng = numAll(3);            % number of generators
nvars = 2*nb + 2*ng+1;           % number of variables
nxtra = nvars - 2*nb;          % number of variables excluding Va and Vm

f=x(nvars);
% ====================================================================
% Evaluate gradient vector
% Gradient is the vector of first derivatives of the objective function
% ====================================================================
 df = sparse(nvars,1,1,nvars,1);


