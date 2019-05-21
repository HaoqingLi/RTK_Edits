function [Optis,fun] = search(L,d,a,p)
%
% [Optis,fun] = search(L,d,a,p) produces the p optimal solutions to the ILS
%             problem  min_{z}f(z)=(z-a)'*(L'DL)^{-1}*(z-a) by search.
%
% Input arguments:
%    L ---- n by n unit lower triangular matrix.
%    d ---- n-dimensional real positive vector, d = diag(D).
%    a ---- n-dimensional real vector.
%    p ---- number of required optimal solutions with default value 1
%
% Output arguments:
%    Optis -- n by p matrix whose j-th column is the j-th optimal solution.
%    fun ---- p-dimensional real vector for the p-values of the objective 
%             function at the p optimal solutions, referred to as residuals 
%             or distances, i.e., fun(j)=f(Optis(:,j)).

% ----------------------------------------------------------------------
% Copyright (c) 2011. Xiao-Wen Chang and Xiaohu Xie
% Version 1.0, September 2011.
% ----------------------------------------------------------------------

n = length(d);

% Initialization 
z = zeros(n,1);      % current candidate solution
zb = zeros(n,1);     % used to compute z
step = zeros(n,1);   % step(k) is used to enumerate z(k) at level k
dist = zeros(n,1);   % partial distance of a candidate solution
S = zeros(n);        % used for computing zb
Optis = zeros(n,p);  % to store the p candidate solutions  
fun = zeros(1,p);      % residuals of the p candidate solutions
imax = p;

% Determine which level to move to after z(1) is chosen at level 1.
if p == 1            
    k0 = 2;          
else
    k0 = 1;           
end

maxDist = inf;       % initial ellipsoid bound 
count = 0;           % initial number of candidate solutions
endSearch = false;   % termination indicator  

zb(n) = a(n);
z(n) = round(zb(n));  
y = zb(n) - z(n);  
step(n) = sign(y); 

k=n;

while ~endSearch

    newDist = dist(k) + y*y/d(k);
    if newDist < maxDist
        if k ~= 1      % move to level k-1
            k = k - 1;
            dist(k) = newDist;
            S(k,1:k) = S(k+1,1:k) + (z(k+1)-zb(k+1))*L(k+1,1:k);
            zb(k) = a(k) + S(k,k);
            z(k) = round(zb(k));
            y = zb(k) - z(k);
            step(k) = sign(y);
        else      % store the found point and try next valid integer
            if count < p-1     % store the first p-1 initial points
                count  =  count + 1;
                Optis(:,count) = z;
                fun(count) = newDist;  % store f(z)
            else
                Optis(:,imax) = z;
                fun(imax) = newDist;
                [maxDist,imax] = max(fun);
            end
            k = k0;
            z(k) = z(k) + step(k);   % next valid integer      
            y = zb(k) - z(k);
            step(k) = -step(k) - sign(step(k));
        end
    else   % exit or move to level k+1
        if k == n            
            endSearch = true;
        else
            k = k + 1;      % move to level k+1
            z(k) = z(k) + step(k);    % next integer
            y = zb(k) - z(k);
            step(k) = -step(k) - sign(step(k));
        end
    end
end

% Sort the solutions by their corresponding residuals
[fun,p] = sort(fun);
Optis = Optis(:,p);


