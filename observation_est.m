function [ h ] = observation_est(state, satPos, satRefPos, ambiguitiesN, waveLength, D, basePosition )

n = size(satPos,1);                                                         % number of satellites observed
h = zeros(2*n,1);                                                           % initialization of the expected/observed measurements
d_iR=zeros(2*n,1);  
d_cR=zeros(2*n,1); 
d_iB=zeros(2*n,1);  
d_cB=zeros(2*n,1); 
roverPos    = state(1:3);          % Position of the rover, obtained by suming the position of the base to the baseline between base and rover
for i=1:n
    d_iR(i)   = norm(roverPos' - satPos(i,:));                               % Distance from the receiver to the ith satellite
    d_cR(i)   = norm(roverPos' - satRefPos(i,:));                            % Distance from the reference satellite to the receiver
    d_iR(i+n) = norm(roverPos' - satPos(i,:));                               % Distance from the receiver to the ith satellite
    d_cR(i+n) = norm(roverPos' - satRefPos(i,:));                            % Distance from the reference satellite to the receiver
    
    d_iB(i)   = norm(basePosition - satPos(i,:));                               % Distance from the Basement to the ith satellite
    d_cB(i)   = norm(basePosition - satRefPos(i,:));                            % Distance from the reference satellite to the Basement
    d_iB(i+n) = norm(basePosition - satPos(i,:));                               % Distance from the Basement to the ith satellite
    d_cB(i+n) = norm(basePosition - satRefPos(i,:));                            % Distance from the reference satellite to the Basement
    
    
    h(i)    = waveLength(i)*D(i,:)*ambiguitiesN;
    h(i+n)  = 0;
end
    D_est=d_iR-d_cR-(d_iB-d_cB);
    h=h+D_est;
end