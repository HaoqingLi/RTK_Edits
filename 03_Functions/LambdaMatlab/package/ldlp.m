function [L,d,P,a] = ldlp(W,a)
%
% [L,d,P,a] = ldlp(W,a) computes the L'DL factorization of W with minimum
%             symmetric piovting: W = P'L'DLP and computes a:=Pa.
%
% Input arguments:
%    W ---- n by n symmetric positive definite matrix. 
%    a ---- n-dimensional vector.
%
% Output arguments:
%    L ---- n by n unit lower triangular matrix.
%    d ---- vector formed by the diagonal of the n by n diagonal matrix D.
%    P ---- n by n permutation matrix.
%    a ---- updated a.

% ----------------------------------------------------------------------
% Copyright (c) 2011. Xiao-Wen Chang and Xiaohu Xie
% Version 1.0, September 2011.
% ----------------------------------------------------------------------

n = size (W,1);
% n
perm = 1:n;  % record the permutation information
% perm
for k = n:-1:1
    [mindiag, q] = min(diag(W(1:k,1:k)));
%	if mindiag <= 0
%		error ('Input matrix W is not positive definite!');
%    end
    perm([k q]) = perm([q k]);
    W([k q],:) = W([q k],:);
    W(:,[k q]) = W(:, [q k]);
    W(k,1:k-1) = W(k,1:k-1)/W(k,k);
    W(1:k-1,1:k-1) = W(1:k-1,1:k-1) - W(k,1:k-1)'*W(k,k)*W(k,1:k-1);   
    a([k q]) = a([q k]);
end

d = diag(W)';
L = tril(W,-1) + eye(n);
P = zeros(n);
for j = 1:n
    P(j,perm(j))=1;
end

