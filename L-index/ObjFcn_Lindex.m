function [f,df,d2f] = ObjFcn_Lindex(x)

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
% Problem dimensions
nb = numAll(1);            % number of buses
nl = numAll(2);            % number of branches
ng = numAll(3);            % number of generators
nvars = numAll(4);         % number of variables
nxtra = numAll(5);         % number of variables excluding Va and Vm

f=x(nvars);
% ====================================================================
% Evaluate gradient vector
% Gradient is the vector of first derivatives of the objective function
% ====================================================================
 df = sparse(nvars,1,1,nvars,1);
% ====================================================================
% Evaluate Hessian matrix
% Hessian is the matrix of second derivatives of the objective function
% ====================================================================
d2f = sparse(nvars,nvars);


