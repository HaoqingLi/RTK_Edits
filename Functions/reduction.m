function [L,d,Z,a] = reduction(W,a)
%
% [L,d,Z,a] = reduction(W,a) computes the factorization W=Z'*L'*D*L*Z
%             and a:=Z^{-T}a.
%
% Input arguments:
%    W ---- n by n symmetric positove definite matrix
%    a ---- n-dimensional real vector.
%
% Output arguments:
%    L ---- n by n unit lower triangular matrix.
%    d ---- vector formed by the diagonal of the n by n diagonal matrix D.
%    Z ---- n by n unimodular matrix.
%    a ---- updated a.

% ----------------------------------------------------------------------
% Copyright (c) 2011. Xiao-Wen Chang and Xiaohu Xie
% Version 1.0, September 2011.
% ----------------------------------------------------------------------

n = size(W,2);

% Compute the L'DL factorization of W with minimum symmetric pivoting
[L,d,Z,a] = ldlp(W,a);

k = n - 1;

while k > 0
    
    kp1=k+1;
    
    tmp = L(kp1,k);
    mu = round(tmp);
    if mu ~= 0
        tmp = tmp - mu*L(kp1,kp1);
    end
    
    delta = d(k) + tmp^2*d(kp1);
    
    if delta < d(kp1)
        if mu ~=0  % Perform size reductions
            L(kp1:n,k) = L(kp1:n,k) - mu*L(kp1:n,kp1);
            Z(kp1,:) = Z(kp1,:) + mu*Z(k,:);
            a(k) = a(k) - mu*a(kp1);
            for i = k+2:n
                mu = round(L(i,k));
                if mu ~= 0 
                    L(i:n,k) = L(i:n,k) - mu*L(i:n,i);
                    Z(i,:) = Z(i,:) + mu*Z(k,:);
                    a(k) = a(k) - mu*a(i);
                end
            end
        end

        % Perform a symmetric permutation 
        eta = d(k)/delta;
        lambda = d(kp1)*L(kp1,k)/delta;
        d(k) = eta*d(kp1);
        d(kp1) = delta;
        L(k:kp1,1:k-1) = [-L(kp1,k),1; eta,lambda]*L(k:kp1,1:k-1);
        L(kp1,k) = lambda;
        L(k+2:n,[k kp1]) = L(k+2:n,[kp1 k]);
        Z([k kp1],:) = Z([kp1 k],:);
        a([k kp1]) = a([kp1 k]);

        if k < n-1
            k = k+1;
        end
    else
        k = k-1;
    end
end
