function R_out = ECEF_2_NED_rotation( Reference_position )
%**************************************************************************
%
% A function to rotate the ECEF frame to the local NED frame 
%
%
% Formula is taken from Coordinate Systems and Transformations Eq. 2.24
%
% 23/November 2017. v.1.2: Changed xyz2ell function for the homemade
% function ECEF_2_Ellipsoid, as it requires quite less time and we are not
% concern about the accuracy for this application
%**************************************************************************

[latitude, longitude, ~] = ECEF_2_Ellipsoid(Reference_position);
% [latitude, longitude, ~] = xyz2ell(Reference_position(1),Reference_position(2),Reference_position(3));

R_out = zeros(3,3);

R_out(1,1) = -sin(latitude) * cos(longitude); 
R_out(1,2) = -sin(latitude) * sin(longitude);
R_out(1,3) = cos(latitude);

R_out(2,1) = -sin(longitude);
R_out(2,2) = cos(longitude);
R_out(2,3) = 0;

R_out(3,1) = -cos(latitude) * cos(longitude); 
R_out(3,2) = -cos(latitude) * sin(longitude);
R_out(3,3) = -sin(latitude);
