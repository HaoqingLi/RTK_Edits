function bigM = augmentCovarianceMatrix (M, index)
%**************************************************************************
%
% Date: 09.01.2017
% DLR Neustrelitz
% Author: Daniel Arias Medina
%
% This function is used to augment a covariance matrix (square matrix) "M"
% by adding ones in the diagonal values indicated in the input "index",
% while filling with zeros the cross-related values.
% An example for this:
% 
%  M = [1  2  3   , index = 2 ---> bigM = [1  0  2  3
%       4  5  6                            0  1  0  0
%       7  8  9]                           4  0  5  6
%                                          7  0  8  9] 
%
% ver 0.1 - Basic Implementation
%
%
%**************************************************************************

if det(M) == 0 
    error('Your matrix is not diagonal!')
end

if length(index) < 1 
%     disp('There is nothing to augment!')
    bigM = M;
   return
end


for z = 1:length(index)
    n       = length(M);
    bigM    = eye(n+1);
    c = index(z);
    for i=1:n+1
        for j=1:n+1
            if i<c && j<c
                bigM(i,j) = M(i,j);
            elseif i<c && j>c
                bigM(i,j) = M(i,j-1);
            elseif i>c && j<c
                bigM(i,j) = M(i-1,j);
            elseif i>c && j>c
                bigM(i,j) = M(i-1,j-1);
            end
        end
    end
    M = bigM;
end



