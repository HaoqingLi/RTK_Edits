function [ Jh, D ] = jacobianINSRTK_EKF( state, sizeState, satPos, satRefPos, satRefIndex, wavelengthVector, typeObs, leverArmIMU2GNSS )
%**************************************************************************
%
% Date: 09.01.2017
% DLR Neustrelitz
% Author: Daniel Arias Medina
%
% Function to build the Jacobian of the measurement model for the EKF,
% for the case of having an INS/GNSS (RTK) system, in which the single
% differences of the phase cycles are tracked in the state.
% 
% The measurements being used are:
% - Code and Phase double differences for RTK-based positioning
% - Doppler? 
%
% Additional considerations:
% - The state consists of: x = [ x_INS | x_GNSS ] = [ pos, vel, quat, biasAcc, biasGyr | SD_{1:m}, \dot{clk}? ]
% The observation model is as follows:
% DD_{i,r} = p^i - p^{GNSS} / | p^i - p^{GNSS} |  + SD_i - SD_r
% The issue is when the position of the GNSS antenna does not correspond
% with the tracked position in the EKF. Normally what it is tracked is the
% location of the INS. Therefore, we need to estimate the position of the
% GNSS from the INS one:
% p^{GNSS} = p^{INS} + R_eb * l^{b}_{INS-GNSS}  =  p^{INS} + quat \otimes l^{b}_{INS-GNSS} \otimes quat^*
% Where l^{b}_{INS-GNSS} is the baseline from the INS to the GNSS position,
% measured in the body frame. The rotation applied can be expressed using
% the rotation matrix or the quaternion rotation. 
% 
% Another remark: the quaternion is expressed as [q_w; q_u], being q_w the
% real part and q_u the imaginary part.
% \partial quat \otimes l^{b}_{INS-GNSS}  \ otimes quat^* / \partial quat --->
% d(q x c x q)/ dq = 2 [q_w c + q_u x c | q_u^T c I_3 + q_u c^T - c q_u^T + q_w[c]_x ]
% For more details about the above formula and notation, please refer to
% the "Weekly Presentation" slides for MSS from Daniel Medina or
% "Quaternion Dynamics for Error State KF" from Joan Sola, version 2017.
%
% Similarly, for the estimation of the velocity using the Doppler
% measurements, we have an observation model such as: 
% c {\Delta f}^i / f = (v^i - v)*(p^i-p)|p^i - p| + \dot{clk}
% Where we also suffer from the fact that the position and the velocity
% refer to the GNSS position and not the INS, as they are tracked in the
% state. Here, the influence of the quaternion over the position is very
% very soft (the unit vector to the satellite will not change considering
% while considering the position of the GNSS antenna or the INS), but the
% velocity is, indeed, important:
% v^{GNSS} = v{INS} + R{q}( \omega^b x l^{b}_{INS-GNSS} ) = v{INS} + \quat \otimes( \omega^b x l^{b}_{INS-GNSS} ) \otimes quat^* 
% also, we will face this relatively complicated partial derivative with
% respect to the quaternion, seen a bit before.
% ver 0.1 - Basic Implementation
%
%
%**************************************************************************


n = size(satPos,1);
Jh = zeros(2*n, length(state));
positionINS = state(1:3);
quat = state(7:10);
q_w = quat(1);
q_u = quat(2:4);

positionGNSS = positionINS + quat_rotate( quat, leverArmIMU2GNSS );

SDAmb               = state(sizeState+1:end);
roverPos            = positionGNSS';                                        % Position of the rover, obtained by suming the position of the base to the baseline between base and rover

for i = 1:n  
    d_i(i)             = norm(roverPos - satPos(i,:));                      % normalised LOS vector pointing from the rover to the ith satellite
    d_Ref(i)           = norm(roverPos - satRefPos(i,:));                   % normalised LOS vector pointing from the rover to the reference satellite
    Jh(i,1:3)       = [ (satPos(i,:) - roverPos)/d_i(i) - (satRefPos(i,:) - roverPos)/d_Ref(i) ];
    Jh(i+n,1:3)     = Jh(i,1:3);
    
    d_Ref           = (roverPos - satRefPos(i,:))/norm(roverPos - satRefPos(i,:));  % normalised LOS vector pointing from the rover to the reference satellite
    d_i             = (roverPos - satPos(i,:))/norm(roverPos - satPos(i,:));        % normalised LOS vector pointing from the rover to the ith satellite
    Jh(i,1:3)       = d_i - d_Ref;
    Jh(i+n,1:3)     = d_i - d_Ref;
    
    F_q_bIMU2GNSS = 2 * [ q_w*leverArmIMU2GNSS + cross(q_u, leverArmIMU2GNSS), q_u'*leverArmIMU2GNSS*eye(3) + q_u*leverArmIMU2GNSS' - leverArmIMU2GNSS*q_u' - q_w*skewMatrix(leverArmIMU2GNSS)  ];
    Jh(i,7:10)      =  -(d_i-d_Ref) * F_q_bIMU2GNSS;
    Jh(i+n,7:10)    = -(d_i-d_Ref) * F_q_bIMU2GNSS;
    
end

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
Jh(1:n-max(typeObs),sizeState+1:end)     = Daux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jh = -Jh;
% D = -D;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

