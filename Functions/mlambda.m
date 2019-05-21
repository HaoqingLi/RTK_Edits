function [X,r] = mlambda(W,a,p)
%
% [X,r] = mlambda(W,a,p) produces p optimal solutions to the integer 
%         least squares problem min_{x}(x-a)'*W^{-1}*(x-a) and the 
%         corresponding values of the objective function 
%
% Input arguments:
%    W ---- n by n symmetric positive definite matrix. In GNSS, it is  
%           the covariane matrix of the real least squares estimator  
%           of the integer ambiguity vector.
%    a ---- n-dimensional real vector. In GNSS, it is the real least 
%           squares estimator of the ambiguity vector.
%    p ---- number of required optimal solutions with default value 1.
%
% Output arguments:
%    X ---- n by p matrix. Its j-th column is the j-th optimal solution,
%           i.e., its objective function is the j-th smallest.
%    r ---- p-dimensional real vector for the p values of the objective 
%           function at the p optimal solutions, i.e.,
%           r(j)=(X(:,j)-a)'*W^{-1}*(X(:,j)-a)
%

% ----------------------------------------------------------------------
% Copyright (c) 2011. Xiao-Wen Chang and Xiaohu Xie
% Version 1.0, September 2011.
% ----------------------------------------------------------------------

% Check input arguments
if nargin < 2 % input error
    error('Not enough input arguments!')
end

if nargin < 3
    p = 1;
end

if p <= 0 % input error
    error('Third input argument must be an integer bigger than 0!')
end

%%%%%%%%%%Daniel Medina modification
% if isempty(W) || isempty(a)
%     X = [[], []];
%     r = [1,1];
%     return;
% end
%%%%%%%%%%

[m,n] = size(W);
if m ~= n || m ~= size(a,1) || size(a,2) ~= 1  % input error
    error('Input arguments have a matrix dimension error!')
end

% Transform the ILS problem by reduction
[L,d,Z,a] = reduction(W,a);

% Search to find the p optimal solutions of the reduced ILS problem
[Optis,r] = search(L,d,a,p);

% Find the p optimal solutions of the original ILS problem by transformation
X = Z'*Optis;
