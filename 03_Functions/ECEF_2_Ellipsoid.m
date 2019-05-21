function [ latitude, longitude, height ] = ECEF_2_Ellipsoid( Reference_position )
%**************************************************************************
%
% Date: 10.08.2015
% DLR 
% Author: Daniel Arias Medina
%
% A function to calculate the ellipsoid coordinates of a point from ECEF
% coordinates:
%       (x,y,z) --> (latitude,longitude,height)
%  
% The procedure is explained on the page 186 of the ESA Book: "GNSS Data
% Processing, Vol I: Fundamentals and Algorithms" 
%
% Longitude can be obtained directly, while an iterative method is required
% to estimate the latitude and the height.
%
%
%**************************************************************************

longitude = atan(Reference_position(2)/Reference_position(1));

iterationCounter = 1;
max_iteration = 20;
latitude_change = 1;
latitude_accuracy = 1e-06;
eccentricity = 0.08181919;                                                  % First eccentricity
R_EA = 6378137;                                                             % Semi-major axis

p = sqrt(Reference_position(1)^2+Reference_position(2)^2);

% First value of the latitude for the iteration procedure
latitude = atan( ( Reference_position(3)/p ) / (1-eccentricity^2));

while max_iteration > iterationCounter && latitude_change > latitude_accuracy
    
    latitude_aux = latitude;
    N_e = R_EA/(1-eccentricity^2*sin(latitude)^2);                          % Prime vertical radius of curvature
    height = p/cos(latitude) - N_e;
    latitude = atan( (Reference_position(3)/p) / ( 1-eccentricity^2*( N_e/(N_e+height) ) )   );
    latitude_change = abs(latitude - latitude_aux);
    iterationCounter = iterationCounter + 1;
    
end
