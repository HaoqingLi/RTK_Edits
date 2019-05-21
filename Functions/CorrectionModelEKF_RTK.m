function [ h ] = CorrectionModelEKF_RTK( state, satPos, satRefPos, ambiguitiesN, waveLength, D, baseline )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

n = size(satPos,1);                                                         % number of satellites observed
h = zeros(2*n,1);                                                           % initialization of the expected/observed measurements

roverPos    = state(1:3);          % Position of the rover, obtained by suming the position of the base to the baseline between base and rover
for i=1:n
    d_iR(i)   = norm(roverPos' - satPos(i,:));                               % Distance from the receiver to the ith satellite
    d_cR(i)   = norm(roverPos' - satRefPos(i,:));                            % Distance from the reference satellite to the receiver
    
    
    h(i)    = waveLength(i)*D(i,:)*ambiguitiesN;
    h(i+n)  = 0;
end

end

