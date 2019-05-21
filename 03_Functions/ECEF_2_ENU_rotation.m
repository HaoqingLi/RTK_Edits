function Rotation_Matrix_ECEF2ENU = ECEF_2_ENU_rotation( Reference_position )
%**************************************************************************
%
% Date: 10.08.2015
% DLR 
% Author: Daniel Arias Medina
%
% This functions returns the rotation matrix to accomplish the rotation
% from ECEF to ENU coordinates. It requires the latitude and longitude of
% the reference position. 
% 
% As latitude and longitude values belong to the ellipsoid coordinate
% system, we need to get the values of latitude and longitude using the
% function "ECEF_2_Ellipsoid". 
%  
% The procedure is explained on the page 188 of the ESA Book: "GNSS Data
% Processing, Vol I: Fundamentals and Algorithms" 
%
%
%**************************************************************************

[latitude, longitude, ~] = ECEF_2_Ellipsoid(Reference_position);
% [latitude, longitude, ~] = xyz2ell(Reference_position(1),Reference_position(2),Reference_position(3));


Rotation_Matrix_ECEF2ENU = [           -sin(longitude),                    cos(longitude),                   0         ;...
                             -cos(longitude) * sin(latitude),     -sin(longitude) * sin(latitude),  	cos(latitude) ;...
                              cos(longitude) * cos(latitude),      sin(longitude) * cos(latitude),      sin(latitude) ];
