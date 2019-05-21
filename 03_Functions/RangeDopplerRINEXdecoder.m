function [ pseudorange, phase, doppler, SV_pos, SV_vel, iono_corr, trop_corr, SV_CLK, SV_SNR, SV_elev, Number_of_satellites, N_aux, satLOS, GNSS_type ] = RangeDopplerRINEXdecoder( i, Eph, C1C, L1C, D1C, snr, time_GPS, nSatTot, pos_receiv, phiR, lamR, hR, iono, elevationMask, receiver_CLKoffset, frequency )
% -------------------------------------------------------------------------
%
% Date: 14.12.2015
% DLR 
% Author: Daniel Arias Medina
%
% A function to save some space and extra lines of code in the main loop
% of the scripts, where it is called the pertinent dedicated SPP/filter 
% classes for GNSS.
%  
% To understand the procedure, it is suggested looking at the scripts using
% RINEX, such as "02_Scripts/RUN_Koblenz_2016_07_26.m"
%
% -------------------------------------------------------------------------

if nargin == 12
    elevationMask = 0;
    receiver_CLKoffset = 0;
    frequency = goGNSS.FL1;
elseif nargin == 13
    receiver_CLKoffset = 0;
    frequency = goGNSS.FL1;
elseif nargin == 14
    frequency = goGNSS.FL1;
end

if isnan(receiver_CLKoffset)
    receiver_CLKoffset = 0;
end

sat0 = find(C1C(:,i) ~= 0);
Eph_t = rt_find_eph (Eph, time_GPS(i), nSatTot);    % CHANGES INSIDE???
[XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(time_GPS(i), C1C(sat0,i), sat0, Eph_t, [], [], zeros(nSatTot,1), zeros(nSatTot,1), receiver_CLKoffset);
[az, SV_elev, dist] = topocent(pos_receiv, XS);     % obtain the elevation and azimuth of the satellites
trop_corr = tropo_error_correction(SV_elev,hR*ones(length(SV_elev),1));     % Troposphere corrections
iono_corr  = iono_error_correction(phiR, lamR, az, SV_elev, time_GPS(i), iono, []); % Ionosphere corrections
iono_corr = iono_corr * (goGNSS.FL1/frequency)^2;

% Second iteration to refine the transmission time based on the knowledge of the ionospheric and tropospheric corrections
[XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(time_GPS(i), C1C(sat0,i), sat0, Eph_t, [], [], trop_corr, iono_corr, receiver_CLKoffset);
[az, SV_elev, dist] = topocent(pos_receiv, XS);     % obtain the elevation and azimuth of the satellites
trop_corr = tropo_error_correction(SV_elev,hR*ones(length(SV_elev),1));     % Troposphere corrections
iono_corr  = iono_error_correction(phiR, lamR, az, SV_elev, time_GPS(i), iono, []); % Ionosphere corrections
iono_corr = iono_corr * (goGNSS.FL1/frequency)^2;

lowElevation_satellites = find(abs(SV_elev)<elevationMask);
no_eph(lowElevation_satellites) = 1;

index = find(no_eph == 0);
sat          =    sat0(index);
c1cA         =    C1C(sat0,i);
pseudorange  =    c1cA(index); pseudorange = pseudorange';
doppler      =    D1C(sat,i);
l1cA         =    L1C(sat0,i);
phase        =    l1cA(index); phase = phase';
snrA         =    snr(sat0,i);
SV_SNR          =    snrA(index); SV_SNR = SV_SNR';
SV_elev           =    SV_elev(index);
SV_elev           =    SV_elev*pi/180; SV_elev = SV_elev';
az           =    az(index);
dist         =    dist(index);
XS           =    XS(index,:);
XS_tx        =    XS_tx(index,:);
VS_tx        =    VS_tx(index,:);
SV_pos       =    XS_tx; %SV_pos = SV_pos';
SV_vel       =    VS_tx;
dtS          =    dtS(index);
SV_CLK       =    dtS; SV_CLK = SV_CLK';
sys          =    sys(index);
nsat         =    size(pseudorange,1);
trop_corr    =    trop_corr(index); trop_corr = trop_corr';
iono_corr     =    iono_corr(index); iono_corr = iono_corr';
n            =    length(index);
N_aux = index;
Number_of_satellites = n;
satLOS = sat;

GNSS_type = [length(find(sys==1)); length(find(sys==2)); length(find(sys==3))];



end

