function [ Jh, D ] = BuildJacobianEKF_RTK( state, sizeState, satPos, satRefPos, satRefIndex, wavelengthVector, typeObs, leverArmIMU2GNSS )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if nargin==6,
   leverArmIMU2GNSS = [0,0,0]'; 
end
n = size(satPos,1);
Jh = zeros(2*n, length(state));


SDAmb               = state(7:end);
roverPos            = state(1:3)';                                        % Position of the rover, obtained by suming the position of the base to the baseline between base and rover
for i = 1:n  
    d_i(i)             = norm(roverPos - satPos(i,:));                      % normalised LOS vector pointing from the rover to the ith satellite
    d_Ref(i)           = norm(roverPos - satRefPos(i,:));                   % normalised LOS vector pointing from the rover to the reference satellite
    Jh(i,1:3)       = [ (satPos(i,:) - roverPos)/d_i(i) - (satRefPos(i,:) - roverPos)/d_Ref(i) ];
    Jh(i+n,1:3)     = Jh(i,1:3);
    
    d_Ref           = (roverPos - satRefPos(i,:))/norm(roverPos - satRefPos(i,:));  % normalised LOS vector pointing from the rover to the reference satellite
    d_i             = (roverPos - satPos(i,:))/norm(roverPos - satPos(i,:));        % normalised LOS vector pointing from the rover to the ith satellite
    Jh(i,1:3)       = d_i - d_Ref;
    Jh(i+n,1:3)     = d_i - d_Ref;
    
    
    
%     d_cR   = sqrt( (roverPos(1) - satRefPos(i,1))^2 + (roverPos(2) - satRefPos(i,2))^2 + (roverPos(3) - satRefPos(i,3))^2 );  % normalised LOS vector pointing from the rover to the reference satellite
%     d_cB   = sqrt( (basePos(1) - satRefPos(i,1))^2 + (basePos(2) - satRefPos(i,2))^2 + (basePos(3) - satRefPos(i,3))^2 );  % normalised LOS vector pointing from the rover to the reference satellite
%     d_iR   = sqrt( (roverPos(1) - satPos(i,1))^2 + (roverPos(2) - satPos(i,2))^2 + (roverPos(3) - satPos(i,3))^2 );  % normalised LOS vector pointing from the rover to the reference satellite
%     d_iB   = sqrt( (basePos(1) - satPos(i,1))^2 + (basePos(2) - satPos(i,2))^2 + (basePos(3) - satPos(i,3))^2 );  % normalised LOS vector pointing from the rover to the reference satellite
%     z(i)    = d_iR - d_cR - (d_iB - d_cB) + wavelengthVector(i)*(SDAmb(i) - SDAmb(satRefIndex(typeObs(i))));
%     z(i+n)  = d_iR - d_cR - (d_iB - d_cB);
end

% z([satRefIndex, satRefIndex+n]) = [];

% Creating the matrix D which relates the satellites to the reference satellite
D = eye(n);
Daux = diag(wavelengthVector);
D(satRefIndex,:)    = [];
Daux(satRefIndex,:) = [];
for i=1:max(typeObs)
    if i==1
        Daux([1:1:length(find(typeObs==1))-1], satRefIndex(i)) = -1*wavelengthVector(1);
        D([1:1:length(find(typeObs==1))-1], satRefIndex(i)) = -1;
    else
        tmp = find(typeObs==i);
        Daux([length(find(typeObs<i))-(i-1)+1:length(find(typeObs<i))-(i-1)+1+length(find(typeObs==i))-i],satRefIndex(i)) = -1*wavelengthVector(tmp(1));
        D([length(find(typeObs<i))-(i-1)+1:length(find(typeObs<i))-(i-1)+1+length(find(typeObs==i))-i],satRefIndex(i)) = -1;
    end
end

% Inserting the D matrix into the Jacobian for the observation model
Jh([satRefIndex, n+satRefIndex],:) = [];
Jh(1:n-max(typeObs),7:end)     = Daux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jh = -Jh;
% D = -D;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

