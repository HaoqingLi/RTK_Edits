% -------------------------------------------------------------------------
%                             DLR Neustrelitz
%
%                          Daniel Arias Medina
%
%   GENERAL RTK Solver                                 
%
%  Script designed to ONLY solve the positioning problem of a 
% rover vehicle. It does NOT include any comparison with a reference
%
% TO DO LIST:
%  - addition of GLONASS not yet done
%  - Graphical Interface
%  - Contineous Mode
%  - Better cycle slip detector
%  - Selection of the positioning method: SPP, KF, RTK, ...
%  - Automatic selection of the base position or manual input?
%  - Selection of the time to use & frequency of algorithm trigger
%  - Working with delayed corrections (age of corrections and so on...)
% 
%
% Date: 15.03.2018
% -------------------------------------------------------------------------
clearvars -except C1C_R L1C_R C2C_R L2C_R D1C_R D2C_R S1C_R S2C_R time_GPS_R timeSOW_R week_R date_R pos_receiv_R interval_R codeC1_R C1C_B L1C_B C2C_B L2C_B D1C_B D2C_B S1C_B S2C_B time_GPS_B timeSOW_B week_B date_B pos_receiv_B interval_B codeC1_B Eph iono
clear all;
% close all;
% clc;
dbstop if warning
dbstop if error
addpath 03_Functions
addpath EKF_Classes
addpath Functions
% Global Variables
FreqL1 = 1575.42e6;
FreqL2 = 1227.60e6;
WavelengthL1 = 0.19029367;
WavelengthL2 = 0.24421021;
speedOfLight = 299792458;
dt = 1;

%% Configuration Parameters
roverObservationFile = 'BugBackbord_0845-101360.17O'; 
navigationFile = 'brdc1360.17n';
% Base Reference Position
baseObservationFile  = 'LAND_SAPOS1_0845-101360.17O';
basePosition = [4041839.1018   537121.6018  4888452.5105]; 
% % Base Reference Position
% baseObservationFile  = 'LAND_SAPOS3_0845-101360.17O';
% basePosition = [4047394.9519   514396.2386  4886511.5423]; 


%%%%%%%%%%% This should be (public) options!
useVariationalFilter = 1;
NumberEpochs2Use = 50;%4500;
GPS_flag = 1; GLO_flag = 0; GAL_flag = 0;
elevationMask = 15;
%%%%%%%%%%% These should NOT be (public) options!
useCycleSlipDetection = 0; useCycleSlipDetectionFreqComb = 1;
useInstantaneousMode = 0;
SlipThres   = 0.005;
Ratio2FixAmb= 3;
useFixHold = 0;
useCorrectionDebug = 1;
observationModel = 'elevation-based' ; % 'elevation-based' 'CN0-based'
%%%%%%%%%%%


%%%%%%% Variance models
a_elev = 0.01; % Values from Eling "Development of an RTK-GPS system for precise ... "
b_elev = 0.01; % Values from Eling
f = 100^2;
elevationVarianceModel = @(a,b,El) 2 * (a^2+(b./sin(El)).^2);
a_CN0 = 10; % Values from Kuusniemi
b_CN0 = 150^2; % Values from Kuusniemi
a_CN0 = 0.3; % Values from ION 2019
b_CN0 = 144^2; % Values from ION 2019
cn0VarianceModel = @(a,b,CN0) a + b * 10.^( -CN0/10 );
ScalingPhase =1; 
ScalingCode = f; 
%%%%%%% 



fprintf('ElevationMask: %i\n',elevationMask)
%%%%%%%%% HERE IT SHOULD INDICATE WHICH DATA ARE BEING USED!





%% Loading GNSS Data
[constellations] = initConstellation(GPS_flag,GLO_flag,GAL_flag,0,0,0);
nSatTot = constellations.nEnabledSat;
n_sys = GPS_flag + GLO_flag + GAL_flag; % Number of constellations used
if ~exist('Eph')
    [C1C_R, L1C_R, C2C_R, L2C_R, D1C_R, D2C_R, S1C_R, S2C_R, time_GPS_R, timeSOW_R, week_R, date_R, pos_receiv_R, interval_R, ~, ~, codeC1_R] = load_RINEX_obs(roverObservationFile, constellations);
    [C1C_B, L1C_B, C2C_B, L2C_B, D1C_B, D2C_B, S1C_B, S2C_B, time_GPS_B, timeSOW_B, week_B, date_B, pos_receiv_B, interval_B, ~, ~, codeC1_B] = load_RINEX_obs(baseObservationFile, constellations);
    [Eph, iono] = load_RINEX_nav(navigationFile, constellations, 0);
end
nEpochs = length(date_R);
[phiR, lamR, hR] = cart2geod(pos_receiv_R(1), pos_receiv_R(2), pos_receiv_R(3));
phiR = phiR * 180 / pi;
lamR = lamR * 180 / pi;
dtR = zeros(nEpochs,1);         % receiver clock error
min_nsat_LS = 3 + n_sys;

%% Loading Outlier Data
load('outlier_simulation1.mat')
outlength_C1=4;
outlength_C2=5;
outlength_L1=4;
outlength_L2=3;
outlier_ind=[1:3 5:10];
%% Constructors for the classes

%%%% MOST OF THESE FUNCTIONS AND OPTIONS SHOULD NOT BE VISIBLE AND SHOULD
%%%% BE CONTAINED INSIDE THE CLASS ITSELF // IN THE FUNCTIONS FROM THE
%%%% CLASS

% EKF Class
sizeState   = 6;
initState   = zeros(6,1);                        
P_          = diag([10^2*ones(3,1);5^2*ones(3,1)]);
% Dynamical model for the prediction step
Dynamic_f   = @DynamicModelEKF_consVel;                                     % Constant speed dynamical model for a state compound by position, velocity and the phase ambiguities.
% Jacobian for the prediction step
Jf          = @(dt,state) [eye(3), dt*eye(3), zeros([3,length(state)-6]);...
                           zeros(3),  eye(3), zeros([3,length(state)-6]);...
                           zeros([length(state)-6,6]), eye(length(state)-6)];                
% System Noise
Q_velocity  = (0.1);
Q_amb       = (0.01);
Q_position  = 1000; 
Q_          = blkdiag( Q_position * eye(3) );    %blkdiag( Q_velocity * eye(3) );            
% Jacobian for the system noise
Jq          = @(dt,state) [dt^2/2*eye(3),  zeros([3,length(state)-6]);     dt*eye(3),  zeros([3,length(state)-6]);     zeros([length(state)-6,3]),  dt*eye(length(state)-6)];
% Jacobian and process model for the correction step
Jh          = @BuildJacobianEKF_RTK;
Observ_h    = @CorrectionModelEKF_RTK;                                       

RTK = RTK_Variational(...
        '-basePosition', basePosition,...
        '-sizeState', sizeState,...
        '-dynModel',Dynamic_f,...
        '-measModel',Observ_h,...
        '-FJ',Jf,...
        '-CJ',Jq,...
        '-HJ',Jh,...
        '-Q',Q_, ...
        '-state', initState,...  
        '-P',P_);
    
% WLS Classes 
LS_Class =      WLS_SPP([0,0,0]', 2^2, 0.02^2, 30, 10^-6, 2, [0.60006, 50.6392, 33.8385], 2 );
LS_Class.debugMode=1;
Pos_ECEF_R      = [];
Vel_ECEF_R      = [];
Pos_ECEF_EEKF   = [];
Vel_ECEF_EEKF   = [];    
    
%% Main Loop
CLKoffset_R     = 0;
CLKoffset_B     = 0;
PhaseAmbiguities{1}          = [];
satRefPRNL1     = -2;
satRefPRNL2     = -2;
satRef_ind      = 1;
satPRN          = [];

NumberEpochs2Use = min(NumberEpochs2Use,nEpochs);
Pos_ECEF_EKF = NaN(NumberEpochs2Use, 3);
Pos_ELL_EKF = NaN(NumberEpochs2Use, 3);
Pos_ENU_EKF = NaN(NumberEpochs2Use, 3);
Pos_ECEF_EKFfix = NaN(NumberEpochs2Use, 3);
Pos_ELL_EKFfix = NaN(NumberEpochs2Use, 3);
Pos_ENU_EKFfix = NaN(NumberEpochs2Use, 3);
PhaseAmbiguities= NaN(NumberEpochs2Use,nSatTot*2);
PhaseAmbiguitiesFix= NaN(NumberEpochs2Use,nSatTot*2);
Cov_Amb_EKF = NaN(1,nSatTot*2);
Accuracy_Hor_EKF = NaN(NumberEpochs2Use,1);
numberObservations = [];

dif=zeros(36,NumberEpochs2Use);

Pos_ECEF_EKF    = NaN(NumberEpochs2Use, 3);
Vel_ECEF_EKF    = NaN(NumberEpochs2Use, 3);
Pos_ECEF_EKFfix = NaN(NumberEpochs2Use, 3);
Vel_ECEF_EKFfix = NaN(NumberEpochs2Use, 3);
EulerAngles_EKF = NaN(NumberEpochs2Use, 3);
State_PredEKF   = [];
Cov_Pos_EKF = [];
Cov_Vel_EKF = [];
Cova_PredEKF = [];

for iAux = 1:nSatTot
    freqCombi_memory{iAux} = [];
    phaseCodeL1Combi_memory{iAux} = [];
    phaseCodeL2Combi_memory{iAux} = [];
    cycleSlipCounterL1 = 0;
    cycleSlipCounterL2 = 0;
end

for iGNSS=1:NumberEpochs2Use
    
    % Extracting the data from the RINEX data
    [ C1C_R_aux, L1C_R_aux, D1C_R_aux, SV_pos_R1, SV_vel_R1, iono_corr_R1, trop_corr_R1, SV_CLK_R1, SV_SNR_R1, SV_elev_R1, Number_of_satellites_R1(iGNSS), N_aux_R1, satPRN_R1 ] = RangeDopplerRINEXdecoder( iGNSS, Eph, C1C_R, L1C_R, D1C_R, S1C_R, time_GPS_R, nSatTot, pos_receiv_R, phiR, lamR, hR, iono, elevationMask, 0, goGNSS.FL1 );
    [ C2C_R_aux, L2C_R_aux, D2C_R_aux, SV_pos_R2, SV_vel_R2, iono_corr_R2, trop_corr_R2, SV_CLK_R2, SV_SNR_R2, SV_elev_R2, Number_of_satellites_R2(iGNSS), N_aux_R2, satPRN_R2 ] = RangeDopplerRINEXdecoder( iGNSS, Eph, C2C_R, L2C_R, D2C_R, S2C_R, time_GPS_R, nSatTot, pos_receiv_R, phiR, lamR, hR, iono, elevationMask, 0, goGNSS.FL2 );
    [ C1C_B_aux, L1C_B_aux, D1C_B_aux, SV_pos_B1, SV_vel_B1, iono_corr_B1, trop_corr_B1, SV_CLK_B1, SV_SNR_B1, SV_elev_B1, Number_of_satellites_B1(iGNSS), N_aux_B1, satPRN_B1 ] = RangeDopplerRINEXdecoder( iGNSS, Eph, C1C_B, L1C_B, D1C_B, S1C_B, time_GPS_B, nSatTot, pos_receiv_B, phiR, lamR, hR, iono, elevationMask, 0, goGNSS.FL1 );
    [ C2C_B_aux, L2C_B_aux, D2C_B_aux, SV_pos_B2, SV_vel_B2, iono_corr_B2, trop_corr_B2, SV_CLK_B2, SV_SNR_B2, SV_elev_B2, Number_of_satellites_B2(iGNSS), N_aux_B2, satPRN_B2 ] = RangeDopplerRINEXdecoder( iGNSS, Eph, C2C_B, L2C_B, D2C_B, S2C_B, time_GPS_B, nSatTot, pos_receiv_B, phiR, lamR, hR, iono, elevationMask, 0, goGNSS.FL2 );
    
    LS_Class.calculate_position(C1C_R_aux, iono_corr_R1, trop_corr_R1, goGNSS.V_LIGHT*SV_CLK_R1, SV_pos_R1(:,1)', SV_pos_R1(:,2)', SV_pos_R1(:,3)', SV_SNR_R1, SV_elev_R1);
    if norm( RTK.state_(1:RTK.sizeState_) -  initState ) < 10^-4
        RTK.state_(1:3) = LS_Class.get_Position();
    end
    
    % Cycle Slip Detection
    in_function_cycleSlip(useCycleSlipDetection);
    
  
          
%%  EKF - RTK   
    % 1) Prediction step
    tmp = RTK.prediction(dt);
    
    
    % 2) GPS L1) Selection of the GPS satellites in L1 observed by the base and rover simultaneously
    satLOSL1            = intersect(satPRN_R1,satPRN_B1); satUsedInd_R1 = ismember(satPRN_R1,satLOSL1); satUsedInd_B1 = ismember(satPRN_B1,satLOSL1);
    if length(satLOSL1)<2 % Check whether there is at least two satellites for this frequency
        DD_C1 = [];  DD_L1 = []; satRefPRNL1 = []; satPosL1 = []; satRefPosL1 = []; satPRNL1 = []; satUsedInd_R1 = []; delta_L1= [];
    else
        if isempty(intersect(satRefPRNL1,satPRN_R1))                              % Find the reference satellite for GPS L1
            [~,satRefL1_R]    = max(SV_elev_R1);
            satRefPRNL1     = satPRN_R1(satRefL1_R);
        end
        C1C_R_corr = C1C_R_aux - iono_corr_R1 - trop_corr_R1 + goGNSS.V_LIGHT*SV_CLK_R1;
        L1C_R_corr = WavelengthL1 * L1C_R_aux + iono_corr_R1 - trop_corr_R1 + goGNSS.V_LIGHT*SV_CLK_R1;
        C1C_B_corr = C1C_B_aux - iono_corr_B1 - trop_corr_B1 + goGNSS.V_LIGHT*SV_CLK_B1;
        L1C_B_corr = WavelengthL1 * L1C_B_aux + iono_corr_B1 - trop_corr_B1 + goGNSS.V_LIGHT*SV_CLK_B1;
        satRefL1_R          = find(satPRN_R1==satRefPRNL1);
        satRefL1_B          = find(satPRN_B1==satRefPRNL1);
        satPRNL1            = satPRN_R1(satUsedInd_R1);
        satRefPosL1         = SV_pos_R1(satRefL1_R,:);
        satPosL1            = SV_pos_R1(satUsedInd_R1,:);
        
        DDIonoL1            = iono_corr_R1(satUsedInd_R1) - iono_corr_R1(satRefL1_R) - (iono_corr_B1(satUsedInd_B1) - iono_corr_B1(satRefL1_R));
        DDTropL1            = trop_corr_R1(satUsedInd_R1) - trop_corr_R1(satRefL1_R) - (trop_corr_B1(satUsedInd_B1) - trop_corr_B1(satRefL1_R));
        DDSvCLKL1           = SV_CLK_R1(satUsedInd_R1) - SV_CLK_R1(satRefL1_R) - (SV_CLK_B1(satUsedInd_B1) - SV_CLK_B1(satRefL1_R));
        DDRangeL1           = C1C_R_corr(satUsedInd_R1) - C1C_R_corr(satRefL1_R) - (C1C_B_corr(satUsedInd_B1) - C1C_B_corr(satRefL1_R))- DDIonoL1 - DDTropL1 + DDSvCLKL1;
        DDPhaseL1           = L1C_R_corr(satUsedInd_R1) - L1C_R_corr(satRefL1_R) - (L1C_B_corr(satUsedInd_B1) - L1C_B_corr(satRefL1_R)) + DDIonoL1 - DDTropL1 + DDSvCLKL1;
        
        ReceiverPos  = RTK.getTrackedPosition();
        d_iRef_RL1          = sqrt((ReceiverPos(1)-SV_pos_R1(satUsedInd_R1,1)').^2 + (ReceiverPos(2)-SV_pos_R1(satUsedInd_R1,2)').^2 + (ReceiverPos(3)-SV_pos_R1(satUsedInd_R1,3)').^2 ) - sqrt((ReceiverPos(1)-SV_pos_R1(satRefL1_R,1)').^2 + (ReceiverPos(2)-SV_pos_R1(satRefL1_R,2)').^2 + (ReceiverPos(3)-SV_pos_R1(satRefL1_R,3)').^2 );
        d_iRef_BL1          = sqrt((basePosition(1)-SV_pos_B1(satUsedInd_B1,1)').^2 + (basePosition(2)-SV_pos_B1(satUsedInd_B1,2)').^2 + (basePosition(3)-SV_pos_B1(satUsedInd_B1,3)').^2 ) - sqrt((basePosition(1)-SV_pos_B1(satRefL1_B,1)').^2 + (basePosition(2)-SV_pos_B1(satRefL1_B,2)').^2 + (basePosition(3)-SV_pos_B1(satRefL1_B,3)').^2 );
        
        DDRangeL1_d         = DDRangeL1 + d_iRef_BL1 - d_iRef_RL1;
        DDPhaseL1_d         = DDPhaseL1 + d_iRef_BL1 - d_iRef_RL1;
        
        shortenObsC1_R      = C1C_R_corr(satUsedInd_R1)' - sqrt( (RTK.state_(1) - SV_pos_R1(satUsedInd_R1,1)).^2 + (RTK.state_(2) - SV_pos_R1(satUsedInd_R1,2)).^2 + (RTK.state_(3) - SV_pos_R1(satUsedInd_R1,3)).^2 ) - CLKoffset_R(end);
        shortenObsC1_B      = C1C_B_corr(satUsedInd_B1)' - sqrt( (basePosition(1) - SV_pos_B1(satUsedInd_B1,1)).^2 + (basePosition(2) - SV_pos_B1(satUsedInd_B1,2)).^2 + (basePosition(3) - SV_pos_B1(satUsedInd_B1,3)).^2 ) - CLKoffset_B(end);
        shortenObsL1_R      = L1C_R_corr(satUsedInd_R1)' - sqrt( (RTK.state_(1) - SV_pos_R1(satUsedInd_R1,1)).^2 + (RTK.state_(2) - SV_pos_R1(satUsedInd_R1,2)).^2 + (RTK.state_(3) - SV_pos_R1(satUsedInd_R1,3)).^2 ) - CLKoffset_R(end);
        shortenObsL1_B      = L1C_B_corr(satUsedInd_B1)' - sqrt( (basePosition(1) - SV_pos_B1(satUsedInd_B1,1)).^2 + (basePosition(2) - SV_pos_B1(satUsedInd_B1,2)).^2 + (basePosition(3) - SV_pos_B1(satUsedInd_B1,3)).^2 ) - CLKoffset_B(end);
        
        range_1_R=sqrt( (RTK.state_(1) - SV_pos_R1(satUsedInd_R1,1)).^2 + (RTK.state_(2) - SV_pos_R1(satUsedInd_R1,2)).^2 + (RTK.state_(3) - SV_pos_R1(satUsedInd_R1,3)).^2 );
        range_1_B=sqrt( (basePosition(1) - SV_pos_B1(satUsedInd_B1,1)).^2 + (basePosition(2) - SV_pos_B1(satUsedInd_B1,2)).^2 + (basePosition(3) - SV_pos_B1(satUsedInd_B1,3)).^2 );

        SD_C1               = shortenObsC1_R - shortenObsC1_B;
        SD_L1               = shortenObsL1_R - shortenObsL1_B;
        
        DD_C1               = shortenObsC1_R - shortenObsC1_R(satRefL1_R) - (shortenObsC1_B - shortenObsC1_B(satRefL1_R));
        DD_L1               = shortenObsL1_R - shortenObsL1_R(satRefL1_R) - (shortenObsL1_B - shortenObsL1_B(satRefL1_R));
    % shortenObsL1_B=ObsL1_B-delta_L1
        delta_L1=range_1_R-range_1_R(satRefL1_R)-range_1_B+range_1_B(satRefL1_R);
    end
    
    % 2) GPS L2) Selection of the GPS satellites in L2 observed by the base and rover simultaneously
    satLOSL2            = intersect(satPRN_R2,satPRN_B2); satUsedInd_R2 = ismember(satPRN_R2,satLOSL2); satUsedInd_B2 = ismember(satPRN_B2,satLOSL2);
    if length(satLOSL2)<2 % Check whether there is at least two satellites for this frequency
        DD_C2 = [];  DD_L2 = []; satRefPRNL2 = []; satPosL2 = []; satRefPosL2 = []; satPRNL2 = []; satUsedInd_R2 = []; delta_L2= [];
    else
        if isempty(intersect(satRefPRNL2,satPRN_R2))                              % Find the reference satellite for GPS L2
            [~,satRefL2_R]    = max(SV_elev_R2);
            satRefPRNL2     = satPRN_R2(satRefL2_R);
        end
        C2C_R_corr = C2C_R_aux - iono_corr_R2 - trop_corr_R2 + goGNSS.V_LIGHT*SV_CLK_R2;
        L2C_R_corr = WavelengthL2 * L2C_R_aux + iono_corr_R2 - trop_corr_R2 + goGNSS.V_LIGHT*SV_CLK_R2;
        C2C_B_corr = C2C_B_aux - iono_corr_B2 - trop_corr_B2 + goGNSS.V_LIGHT*SV_CLK_B2;
        L2C_B_corr = WavelengthL2 * L2C_B_aux + iono_corr_B2 - trop_corr_B2 + goGNSS.V_LIGHT*SV_CLK_B2;
        satRefL2_R          = find(satPRN_R2==satRefPRNL2);
        satRefL2_B          = find(satPRN_B2==satRefPRNL2);
        satPRNL2            = satPRN_R2(satUsedInd_R2);
        satRefPosL2         = SV_pos_R2(satRefL2_R,:);
        satPosL2            = SV_pos_R2(satUsedInd_R2,:);
        
        DDIonoL2            = iono_corr_R2(satUsedInd_R2) - iono_corr_R2(satRefL2_R) - (iono_corr_B2(satUsedInd_B2) - iono_corr_B2(satRefL2_R));
        DDTropL2            = trop_corr_R2(satUsedInd_R2) - trop_corr_R2(satRefL2_R) - (trop_corr_B2(satUsedInd_B2) - trop_corr_B2(satRefL2_R));
        DDSvCLKL2           = SV_CLK_R2(satUsedInd_R2) - SV_CLK_R2(satRefL2_R) - (SV_CLK_B2(satUsedInd_B2) - SV_CLK_B2(satRefL2_R));
        DDRangeL2           = C2C_R_corr(satUsedInd_R2) - C2C_R_corr(satRefL2_R) - (C2C_B_corr(satUsedInd_B2) - C2C_B_corr(satRefL2_R)) - DDIonoL2 - DDTropL2 + DDSvCLKL2;
        DDPhaseL2           = L2C_R_corr(satUsedInd_R2) - L2C_R_corr(satRefL2_R) - (L2C_B_corr(satUsedInd_B2) - L2C_B_corr(satRefL2_R)) + DDIonoL2 - DDTropL2 + DDSvCLKL2;
        
        ReceiverPos  = RTK.getTrackedPosition();
        d_iRef_RL2          = sqrt((ReceiverPos(1)-SV_pos_R2(satUsedInd_R2,1)').^2 + (ReceiverPos(2)-SV_pos_R2(satUsedInd_R2,2)').^2 + (ReceiverPos(3)-SV_pos_R2(satUsedInd_R2,3)').^2 ) - sqrt((ReceiverPos(1)-SV_pos_R2(satRefL2_R,1)').^2 + (ReceiverPos(2)-SV_pos_R2(satRefL2_R,2)').^2 + (ReceiverPos(3)-SV_pos_R2(satRefL2_R,3)').^2 );
        d_iRef_BL2          = sqrt((basePosition(1)-SV_pos_B1(satUsedInd_B2,1)').^2 + (basePosition(2)-SV_pos_B1(satUsedInd_B2,2)').^2 + (basePosition(3)-SV_pos_B1(satUsedInd_B2,3)').^2 ) - sqrt((basePosition(1)-SV_pos_B1(satRefL2_B,1)').^2 + (basePosition(2)-SV_pos_B1(satRefL2_B,2)').^2 + (basePosition(3)-SV_pos_B1(satRefL2_B,3)').^2 );
        
        DDRangeL2_d           = DDRangeL2 + d_iRef_BL2 - d_iRef_RL2;
        DDPhaseL2_d           = DDPhaseL2 + d_iRef_BL2 - d_iRef_RL2;
        
        shortenObsC2_R      = C2C_R_corr(satUsedInd_R2)' - sqrt( (RTK.state_(1) - SV_pos_R2(satUsedInd_R2,1)).^2 + (RTK.state_(2) - SV_pos_R2(satUsedInd_R2,2)).^2 + (RTK.state_(3) - SV_pos_R2(satUsedInd_R2,3)).^2 ) - CLKoffset_R(end);
        shortenObsC2_B      = C2C_B_corr(satUsedInd_B2)' - sqrt( (basePosition(1) - SV_pos_B2(satUsedInd_B2,1)).^2 + (basePosition(2) - SV_pos_B2(satUsedInd_B2,2)).^2 + (basePosition(3) - SV_pos_B1(satUsedInd_B2,3)).^2 ) - CLKoffset_B(end);
        shortenObsL2_R      = L2C_R_corr(satUsedInd_R2)' - sqrt( (RTK.state_(1) - SV_pos_R2(satUsedInd_R2,1)).^2 + (RTK.state_(2) - SV_pos_R2(satUsedInd_R2,2)).^2 + (RTK.state_(3) - SV_pos_R2(satUsedInd_R2,3)).^2 ) - CLKoffset_R(end);
        shortenObsL2_B      = L2C_B_corr(satUsedInd_B2)' - sqrt( (basePosition(1) - SV_pos_B2(satUsedInd_B2,1)).^2 + (basePosition(2) - SV_pos_B2(satUsedInd_B2,2)).^2 + (basePosition(3) - SV_pos_B1(satUsedInd_B2,3)).^2 ) - CLKoffset_B(end);
        
        range_2_R=sqrt( (RTK.state_(1) - SV_pos_R2(satUsedInd_R2,1)).^2 + (RTK.state_(2) - SV_pos_R2(satUsedInd_R2,2)).^2 + (RTK.state_(3) - SV_pos_R2(satUsedInd_R2,3)).^2 );
        range_2_B=sqrt( (basePosition(1) - SV_pos_B2(satUsedInd_B2,1)).^2 + (basePosition(2) - SV_pos_B2(satUsedInd_B2,2)).^2 + (basePosition(3) - SV_pos_B1(satUsedInd_B2,3)).^2 );


%         shortenObsC2_R      = C2C_R_corr(satUsedInd_R2)' - CLKoffset_R(end);
%         shortenObsC2_B      = C2C_B_corr(satUsedInd_B2)' - CLKoffset_B(end);
%         shortenObsL2_R      = L2C_R_corr(satUsedInd_R2)'  - CLKoffset_R(end);
%         shortenObsL2_B      = L2C_B_corr(satUsedInd_B2)'  - CLKoffset_B(end);
       
        
        SD_C2               = shortenObsC2_R - shortenObsC2_B;
        SD_L2               = shortenObsL2_R - shortenObsL2_B;
        
        DD_C2               = shortenObsC2_R - shortenObsC2_R(satRefL2_R) - (shortenObsC2_B - shortenObsC2_B(satRefL2_R));
        DD_L2               = shortenObsL2_R - shortenObsL2_R(satRefL2_R) - (shortenObsL2_B - shortenObsL2_B(satRefL2_R));
    % shortenObsL2_B=ObsL2_B-delta_L2
        delta_L2=range_2_R-range_2_R(satRefL2_R)-range_2_B+range_2_B(satRefL2_R);
    
    end
    
    DD_C1(outlier_ind(1:outlength_C1))=DD_C1(outlier_ind(1:outlength_C1))+outlier_C1(iGNSS,outlier_ind(1:outlength_C1))';
    DD_C2(outlier_ind(1:outlength_C2))=DD_C2(outlier_ind(1:outlength_C2))+outlier_C2(iGNSS,outlier_ind(1:outlength_C2))';
    DD_L1(outlier_ind(1:outlength_L1))=DD_L1(outlier_ind(1:outlength_L1))+outlier_L1(iGNSS,outlier_ind(1:outlength_L1))';
    DD_L2(outlier_ind(1:outlength_L2))=DD_L2(outlier_ind(1:outlength_L2))+outlier_L2(iGNSS,outlier_ind(1:outlength_L2))';
    
    
    
    % 3) Pile up all the data from different constellations/frequencies into single variables for the position of the satellites, pseudoranges, phases...
    satPos              = [satPosL1; satPosL2];
    satRefPos           = [repmat(satRefPosL1,length(satPRNL1),1); repmat(satRefPosL2,length(satPRNL2),1)];
    wavelengthVector    = [WavelengthL1*ones(length(satPRNL1),1); WavelengthL2*ones(length(satPRNL2),1)];
    typeObs             = [ones(length(satPRNL1),1); 2*ones(length(satPRNL2),1)];
    DDPhase             = [DDPhaseL1_d, DDPhaseL2_d];
    DDRange             = [DDRangeL1_d, DDRangeL2_d];
    satPRN              = [satPRNL1; satPRNL2 + 100];
    satRefPRN           = [satRefPRNL1; satRefPRNL2 + 100 ];                 % just add a +100 to see the PRN of the satellites of GPS L2 frequency
    
    DDRange             = [DD_C1; DD_C2];
    DDPhase             = [DD_L1; DD_L2];
    DDDelta             = [delta_L1; delta_L2];
    
    numberObservations  = [numberObservations; [ length(DD_C1), length(DD_C2), length(DD_C1)+length(DD_C2)] ];
    
    numberObSatellitesUsedL1(iGNSS) = length(DD_C1);
    numberObSatellitesUsedL2(iGNSS) = length(DD_C2);
    
    % 4) Variance model for the observations
    switch observationModel
        case 'elevation-based' 
            varianceObs         = [elevationVarianceModel( a_elev, b_elev, SV_elev_R1(satUsedInd_R1) ), elevationVarianceModel( a_elev, b_elev, SV_elev_R2(satUsedInd_R2) )];
        case 'CN0-based'
            varianceObs         = [cn0VarianceModel( a_CN0, b_CN0, SV_SNR_R1(satUsedInd_R1) ), cn0VarianceModel( a_CN0, b_CN0, SV_SNR_R2(satUsedInd_R2) )];
        otherwise 
            varianceObs         = [elevationVarianceModel( a_elev, b_elev, SV_elev_R1(satUsedInd_R1) ), elevationVarianceModel( a_elev, b_elev, SV_elev_R2(satUsedInd_R2) )];
    end
    RTK.R_              = diag([ScalingPhase*varianceObs, ScalingCode*varianceObs]);
    
    
     % 5) Ambiguity reset: Cycle slip detection! Instantaneous or Contineous mode:
    if useInstantaneousMode
        RTK.cycleSlipReset();
    end % End of the instanteous mode

    % 6) CORRECTION STEP
    if useVariationalFilter == true
        [ filterOutput]        = RTK.correctionNonINSRTK_Variational(satPRN, satRefPRN, satPos, satRefPos, DDPhase, DDRange, wavelengthVector, typeObs,iGNSS,DDDelta)';
        z_i_(iGNSS)=filterOutput(end);
        filterOutput=filterOutput(1,1:length(filterOutput)-1);
    else % Use the regular EKF
        [ filterOutput]        = RTK.correctionNonINSRTK(satPRN, satRefPRN, satPos, satRefPos, DDPhase, DDRange, wavelengthVector, typeObs,iGNSS,DDDelta)';
    end
    
    % 7) MLAMBDA Algorithm for ambiguity fixing
    % A) Integer Ambiguity estimation
    if length(satPRN) >= 2
        [DDAmb, r] = mlambda(RTK.QN_,RTK.DDAmb_,2);
    else
        r = [1,1];
        DDAmb = NaN(length([satPRNL1(satPRNL1~=satRefPRNL1);satPRNL2(satPRNL2~=satRefPRNL2)+33]),1);
    end
    % B) Fix ratio Test
    fixRatio(iGNSS) = r(2)/r(1); 
    % C) Fixed solution 
    if fixRatio(iGNSS)>Ratio2FixAmb                                             % Getting a fixed solution
        fixSolution = filterOutput(1:sizeState)' - RTK.QRN_ * inv(RTK.QN_)*(RTK.DDAmb_-DDAmb(:,1));
        if useFixHold == 1
            tmp = RTK.FixNHoldMode( DDAmb(:,1) );
        else
            RTK.state_(1:sizeState) = fixSolution';
        end
    end
    
    % Saving the estimation for the position, velocity and ambiguities
    nObs = length(fixRatio);
    Pos_ECEF_EKF(iGNSS,:) = RTK.state_(1:3)';
    [Pos_ELL_EKF(iGNSS,1), Pos_ELL_EKF(iGNSS,2), Pos_ELL_EKF(iGNSS,3)]    = cart2geod(Pos_ECEF_EKF(iGNSS,1), Pos_ECEF_EKF(iGNSS,2), Pos_ECEF_EKF(iGNSS,3)); 
    Pos_ELL_EKF(iGNSS,1:2) = Pos_ELL_EKF(iGNSS,1:2)*180/pi;
    Rotation_Matrix_ECEF2ENU = ECEF_2_ENU_rotation( Pos_ELL_EKF(iGNSS,:) );
    if ~exist('ENU_Reference_point')
        ENU_Reference_point = Pos_ECEF_EKF(iGNSS,:);
    end
    Pos_ENU_EKF(iGNSS,:) = [ Rotation_Matrix_ECEF2ENU * ( Pos_ECEF_EKF(end,:) - ENU_Reference_point)' ]';
    Vel_ECEF_EKF(iGNSS,:) = filterOutput(4:6)';
    PhaseAmbiguitiesFix(iGNSS,[satPRNL1(satPRNL1~=satRefPRNL1);satPRNL2(satPRNL2~=satRefPRNL2)+33]) = DDAmb(:,1);
    if ~isempty(satPRNL1) && ~isempty(satPRNL2)
        PhaseAmbiguities(iGNSS,[satPRNL1(satPRNL1~=satRefPRNL1);satPRNL2(satPRNL2~=satRefPRNL2)+33]) = RTK.DDAmb_;
    end
    if fixRatio(end) > Ratio2FixAmb
        Pos_ECEF_EKFfix(iGNSS,:) = RTK.state_(1:3)';
        [Pos_ELL_EKFfix(iGNSS,1), Pos_ELL_EKFfix(iGNSS,2), Pos_ELL_EKFfix(iGNSS,3)]    = cart2geod(Pos_ECEF_EKFfix(iGNSS,1), Pos_ECEF_EKFfix(iGNSS,2), Pos_ECEF_EKFfix(iGNSS,3)); 
        Pos_ELL_EKFfix(iGNSS,1:2) = Pos_ELL_EKFfix(iGNSS,1:2)*180/pi;
        Pos_ENU_EKFfix(iGNSS,:) = [ Rotation_Matrix_ECEF2ENU * ( Pos_ECEF_EKFfix(iGNSS,:) - ENU_Reference_point)' ]';
        Vel_ECEF_EKFfix(iGNSS,:) = RTK.state_(4:6)';
    end
    CovMatrixEKF = RTK.P_;
    CovMatrixENUEKF = Rotation_Matrix_ECEF2ENU' * CovMatrixEKF(1:3,1:3) * Rotation_Matrix_ECEF2ENU;
    Accuracy_Hor_EKF(iGNSS) = sqrt(  CovMatrixENUEKF(1,1) + CovMatrixENUEKF(2,2)  );


end



%% Visualization of the results
% Estimated EKF fix and non-fixed + map
figure; hold on; set(gcf,'color','w'); grid on;  box on;
ax = gca; ax.ColorOrderIndex = 2;
plot(Pos_ELL_EKF(1:nObs,2), Pos_ELL_EKF(1:nObs,1),'.')
ax = gca; ax.ColorOrderIndex = 5;
plot(Pos_ELL_EKFfix(1:nObs,2),Pos_ELL_EKFfix(1:nObs,1),'.')
axis off
% plot_google_map( 'scale',2,'maptype','roadmap','showlabels',0,'mapscale',1, 'APIKey', 'AIzaSyCWyF15LwrY5ftaMNr_PQICwKwXo9NLj0c'  )
leg = legend('Float solution', 'Fix solution');
set(leg, 'interpreter','latex','fontsize',12);

% ENU Solution
figure; hold on; set(gcf,'color','w'); grid on;  box on;
ax = gca; ax.ColorOrderIndex = 2;
plot(Pos_ENU_EKF(1:nObs,1), Pos_ENU_EKF(1:nObs,2),'.-')
ax = gca; ax.ColorOrderIndex = 5;
plot(Pos_ENU_EKFfix(1:nObs,1),Pos_ENU_EKFfix(1:nObs,2),'.-')
leg = legend('Float solution', 'Fix solution');
set(leg, 'interpreter','latex','fontsize',12);
xlabel('East [m]','Fontsize',13,'interpreter','latex')
ylabel('North [m]','Fontsize',13,'interpreter','latex')
set(gca,'FontSize',13,'TickLabelInterpreter','latex');

% ENU solution
figure; hold on; set(gcf,'color','w'); grid on;  box on;
subplot(3,1,1); hold on; set(gcf,'color','w'); grid on;  box on; axis tight;
ax = gca; ax.ColorOrderIndex = 2;
plot(Pos_ENU_EKF(:,1),'.-')
ax = gca; ax.ColorOrderIndex = 5;
plot(Pos_ENU_EKFfix(:,1),'.-')
leg = legend('Float solution', 'Fix solution');
set(leg, 'interpreter','latex','fontsize',12);
ylabel('East [m]','Fontsize',13,'interpreter','latex')
set(gca,'FontSize',13,'TickLabelInterpreter','latex');
subplot(3,1,2); hold on; set(gcf,'color','w'); grid on;  box on; axis tight;
ax = gca; ax.ColorOrderIndex = 2;
plot(Pos_ENU_EKF(:,2),'.-')
ax = gca; ax.ColorOrderIndex = 5;
plot(Pos_ENU_EKFfix(:,2),'.-')
ylabel('North [m]','Fontsize',13,'interpreter','latex')
set(gca,'FontSize',13,'TickLabelInterpreter','latex');
subplot(3,1,3); hold on; set(gcf,'color','w'); grid on;  box on; axis tight;
ax = gca; ax.ColorOrderIndex = 2;
plot(Pos_ENU_EKF(:,3),'.-')
ax = gca; ax.ColorOrderIndex = 5;
plot(Pos_ENU_EKFfix(:,3),'.-')
ylabel('Up [m]','Fontsize',13,'interpreter','latex')
xlabel('time [s]','Fontsize',13,'interpreter','latex')
set(gca,'FontSize',13,'TickLabelInterpreter','latex');


% Fix ratio of the RTK -> subplot mode
figure; hold on; set(gcf,'color','w'); grid on; axis tight; box on;
title('Fix ratio')
plot(Ratio2FixAmb*ones(nObs,1),'--','linewidth',2);
plot(fixRatio,'.-', 'Color',[0.929411768913269 0.694117665290833 0.125490203499794]);
leg = legend('Threshold Value','Fix Ratio');
set(leg, 'interpreter','latex','fontsize',12);
xlabel('time [s]','Fontsize',13,'interpreter','latex')
ylabel('ratio','Fontsize',13,'interpreter','latex')
set(gca,'FontSize',13,'TickLabelInterpreter','latex');

% Estimated velocity
figure; 
title('Velocity Estimation')
subplot(3,1,1); hold on; set(gcf,'color','w'); grid on; axis tight; box on;
title('Velocity North-East-Down')
ax = gca; ax.ColorOrderIndex = 2;
plot(Vel_ECEF_EKF(:,1),'.-')
ax = gca; ax.ColorOrderIndex = 5;
plot(Vel_ECEF_EKFfix(:,1),'.-')
leg = legend('Float solution', 'Fix solution');
set(leg, 'interpreter','latex','fontsize',12);
ylabel('$v_N$ [m/s]','Fontsize',13,'interpreter','latex')
set(gca,'FontSize',13,'TickLabelInterpreter','latex');
subplot(3,1,2); hold on; set(gcf,'color','w'); grid on; axis tight; box on;
ax = gca; ax.ColorOrderIndex = 2;
plot(Vel_ECEF_EKF(:,2),'.-')
ax = gca; ax.ColorOrderIndex = 5;
plot(Vel_ECEF_EKFfix(:,2),'.-')
ylabel('$v_E$ [m/s]','Fontsize',13,'interpreter','latex')
set(gca,'FontSize',13,'TickLabelInterpreter','latex');
subplot(3,1,3); hold on; set(gcf,'color','w'); grid on; axis tight; box on;
ax = gca; ax.ColorOrderIndex = 2;
plot(Vel_ECEF_EKF(:,3),'.-')
ax = gca; ax.ColorOrderIndex = 5;
plot(Vel_ECEF_EKFfix(:,3),'.-')
xlabel('time [s]','Fontsize',13,'interpreter','latex')
ylabel('$v_D$ [m/s]','Fontsize',13,'interpreter','latex')
set(gca,'FontSize',13,'TickLabelInterpreter','latex');

% Double-difference ambiguities of the satellites
figure; hold on; set(gcf,'color','w'); grid on; axis tight; box on;
plot(PhaseAmbiguitiesFix,'--')
ax = gca; ax.ColorOrderIndex = 1;
plot(PhaseAmbiguities,'-')
xlabel('time [s]','Fontsize',13,'interpreter','latex')
ylabel('DD Amb [cycles]','Fontsize',13,'interpreter','latex')
title('Estimation of ambiguities')
set(gca,'FontSize',13,'TickLabelInterpreter','latex');

% Representation of uncertainty
x = 1:nObs;
y1 = zeros(1,nObs);
y2 = Accuracy_Hor_EKF';
X=[x,fliplr(x)];                %#create continuous x value array for plotting
Y=[y1,fliplr(y2(1:nObs))];              %#create y values for out and then back
figure; hold on; set(gcf,'color','w'); grid on; axis tight; box on;
fill(X,Y,'k', 'FaceColor',[0.862745106220245 0.862745106220245 0.862745106220245]);
ylim([0,4.5*median(y2)])
leg = legend('Accuracy');
set(leg, 'interpreter','latex','fontsize',12);
title('Formal accuracy of the horizontal position')
xlabel('time [s]','Fontsize',13,'interpreter','latex')
ylabel('[m]','Fontsize',13,'interpreter','latex')
set(gca,'FontSize',13,'TickLabelInterpreter','latex');


% Number of satellites observed VS satellites used 
figure; hold on; set(gcf,'color','w'); grid on; axis tight; box on;
subplot(2,1,1); hold on; grid on; axis tight; box on;
title('Satellites observed L1')
plot(Number_of_satellites_R1,'.-')
plot(numberObSatellitesUsedL1,'-')
set(gca,'FontSize',13,'TickLabelInterpreter','latex');
ylabel('satellites','Fontsize',13,'interpreter','latex')
subplot(2,1,2); hold on; grid on; axis tight; box on;
title('Satellites observed L2')
plot(Number_of_satellites_R2,'.-')
plot(numberObSatellitesUsedL2,'-')
xlabel('time [s]','Fontsize',13,'interpreter','latex')
ylabel('satellites','Fontsize',13,'interpreter','latex')

set(gca,'FontSize',13,'TickLabelInterpreter','latex');

figure
 for i=1:iGNSS
covar(i)=sqrt(dif(:,i)'*dif(:,i)/size(dif,1));
end
plot(1:iGNSS,covar)








%% Additional functions
function [] =  in_function_cycleSlip ( useCycleSlipDetection )
if useCycleSlipDetection
    if useCycleSlipDetectionFreqComb
        slideWindow = 10;
        polifitOrder = 2;
        ThresholdCycleSlip = 1.5 * ( WavelengthL2 - WavelengthL1 );
        for iSat = 1:nSatTot
            if L1C_R(iSat,iGNSS)==0 || L2C_R(iSat,iGNSS)==0 % there is an observation missing
                cycleSlipDetected(iSat) = 1;
                freqCombi_memory{iSat} = [];
            else
                phaseFreqDiff = WavelengthL1 * L1C_R(iSat,iGNSS) - WavelengthL2 * L2C_R(iSat,iGNSS);
                lengthSlipMemory = length(freqCombi_memory{iSat});
                if   lengthSlipMemory < slideWindow
                    freqCombi_memory{iSat} = [freqCombi_memory{iSat}, phaseFreqDiff];
                    cycleSlipDetected(iSat) = 0;
                elseif length(freqCombi_memory{iSat}) >= slideWindow
                    p = polyfit([1:lengthSlipMemory], freqCombi_memory{iSat},polifitOrder);
                    preditPhaseFreqDiff = polyval( p, lengthSlipMemory+1 );
                    if abs(preditPhaseFreqDiff - phaseFreqDiff) < ThresholdCycleSlip
                        freqCombi_memory{iSat} = [freqCombi_memory{iSat}, phaseFreqDiff];
                        cycleSlipDetected(iSat) = 0;
                    else
                        %                         figure; hold on; grid on;
                        %                         plot(phaseFreqDiff_memory{iSat},'o')
                        %                         plot( polyval( p, [1:lengthSlipMemory+1] ))
                        %                         plot( lengthSlipMemory+1, polyval( p, slideWindow+1 ), '*')
                        %                         plot(lengthSlipMemory+1,phaseFreqDiff,'s','markerfacecolor','r')
                        %                         legend('Previous phase diff', 'Polifit','Prediction','Current phase diff')
                        freqCombi_memory{iSat} = phaseFreqDiff;
                        cycleSlipDetected(iSat) = 1;
                    end % check whether the current phase difference agrees with the polifit from the memory
                end % Check whether you have enough observations on the memory on the ith satellite
            end % check whether you have measurements
        end % satellite loop
        badSat = find(cycleSlipDetected==1);
        [meas2delete, meas2deleteIndx] = intersect(satPRN_R1,badSat);
        C1C_R_aux(meas2deleteIndx) = []; L1C_R_aux(meas2deleteIndx) = []; D1C_R_aux(meas2deleteIndx) = []; SV_pos_R1(meas2deleteIndx,:) = []; SV_vel_R1(meas2deleteIndx,:) = [];  iono_corr_R1(meas2deleteIndx) = []; trop_corr_R1(meas2deleteIndx) = []; SV_CLK_R1(meas2deleteIndx) = []; SV_SNR_R1(meas2deleteIndx) = []; SV_elev_R1(meas2deleteIndx) = []; satPRN_R1(meas2deleteIndx) = [];  %Number_of_satellites_R1(iGNSS) = length(satPRN_R1);  N_aux_R1 = 1: length(satPRN_R1);
        [meas2delete, meas2deleteIndx] = intersect(satPRN_R2,badSat);
        C2C_R_aux(meas2deleteIndx) = []; L2C_R_aux(meas2deleteIndx) = []; D2C_R_aux(meas2deleteIndx) = []; SV_pos_R2(meas2deleteIndx,:) = []; SV_vel_R2(meas2deleteIndx,:) = [];  iono_corr_R2(meas2deleteIndx) = []; trop_corr_R2(meas2deleteIndx) = []; SV_CLK_R2(meas2deleteIndx) = []; SV_SNR_R2(meas2deleteIndx) = []; SV_elev_R2(meas2deleteIndx) = []; satPRN_R2(meas2deleteIndx) = [];%  Number_of_satellites_R2(iGNSS) = length(satPRN_R2);  N_aux_R2 = 1: length(satPRN_R2);
        [meas2delete, meas2deleteIndx] = intersect(satPRN_B1,badSat);
        C1C_B_aux(meas2deleteIndx) = []; L1C_B_aux(meas2deleteIndx) = []; D1C_B_aux(meas2deleteIndx) = []; SV_pos_B1(meas2deleteIndx,:) = []; SV_vel_B1(meas2deleteIndx,:) = [];  iono_corr_B1(meas2deleteIndx) = []; trop_corr_B1(meas2deleteIndx) = []; SV_CLK_B1(meas2deleteIndx) = []; SV_SNR_B1(meas2deleteIndx) = []; SV_elev_B1(meas2deleteIndx) = []; satPRN_B1(meas2deleteIndx) = [];  %Number_of_satellites_B1(iGNSS) = length(satPRN_B1);  N_aux_B1 = 1: length(satPRN_B1);
        [meas2delete, meas2deleteIndx] = intersect(satPRN_B2,badSat);
        C2C_B_aux(meas2deleteIndx) = []; L2C_B_aux(meas2deleteIndx) = []; D2C_B_aux(meas2deleteIndx) = []; SV_pos_B2(meas2deleteIndx,:) = []; SV_vel_B2(meas2deleteIndx,:) = [];  iono_corr_B2(meas2deleteIndx) = []; trop_corr_B2(meas2deleteIndx) = []; SV_CLK_B2(meas2deleteIndx) = []; SV_SNR_B2(meas2deleteIndx) = []; SV_elev_B2(meas2deleteIndx) = []; satPRN_B2(meas2deleteIndx) = [];  %Number_of_satellites_B2(iGNSS) = length(satPRN_B2);  N_aux_B2 = 1: length(satPRN_B2);
        
    else % Use phase-code combination
        slideWindow = 30;
        polifitOrder = 2;
        ThresholdCycleSlip = 1.5 * ( WavelengthL2 - WavelengthL1 );
        ThresholdCycleSlip = 3;% * 1.798192626953125e+03;
        for iSat = 1:nSatTot
            
            if L1C_R(iSat,iGNSS)==0 || C1C_R(iSat,iGNSS)==0
                cycleSlipL1Detected(iSat) = 1;
            else
                phaseCodeL1Combi =WavelengthL1*L1C_R(iSat,iGNSS) - C1C_R(iSat,iGNSS);
                if length(phaseCodeL1Combi_memory{iSat}) < slideWindow
                    phaseCodeL1Combi_memory{iSat} = [phaseCodeL1Combi_memory{iSat}, phaseCodeL1Combi];
                    cycleSlipL1Detected(iSat) = 0;
                else
                    p = polyfit([1:length(phaseCodeL1Combi_memory{iSat})], phaseCodeL1Combi_memory{iSat},polifitOrder);
                    preditPhaseCodeCombi = polyval( p, length(phaseCodeL1Combi_memory{iSat})+1 );
                    if abs(preditPhaseCodeCombi - phaseCodeL1Combi) < ThresholdCycleSlip
                        phaseCodeL1Combi_memory{iSat} = [phaseCodeL1Combi_memory{iSat}, phaseCodeL1Combi];
                        cycleSlipL1Detected(iSat) = 0;
                    else
                        %                             figure; hold on; grid on;
                        %                             plot(phaseCodeL1Combi_memory{iSat},'o')
                        %                             plot( polyval( p, [1: length(phaseCodeL1Combi_memory{iSat})+1] ))
                        %                             plot(  length(phaseCodeL1Combi_memory{iSat})+1, polyval( p,  length(phaseCodeL1Combi_memory{iSat})+1 ), '*')
                        %                             plot( length(phaseCodeL1Combi_memory{iSat})+1,phaseCodeL1Combi,'s','markerfacecolor','r')
                        %                             legend('Previous phase diff', 'Polifit','Prediction','Current phase diff')
                        phaseCodeL1Combi_memory{iSat} = phaseCodeL1Combi;
                        cycleSlipL1Detected(iSat) = 1;
                    end
                end
            end
            
            if L2C_R(iSat,iGNSS)==0 || C2C_R(iSat,iGNSS)==0
                cycleSlipL2Detected(iSat) = 1;
            else
                phaseCodeL2Combi =WavelengthL2*L2C_R(iSat,iGNSS) - C2C_R(iSat,iGNSS);
                if length(phaseCodeL2Combi_memory{iSat}) < slideWindow
                    phaseCodeL2Combi_memory{iSat} = [phaseCodeL2Combi_memory{iSat}, phaseCodeL2Combi];
                    cycleSlipL2Detected(iSat) = 0;
                else
                    p = polyfit([1:length(phaseCodeL2Combi_memory{iSat})], phaseCodeL2Combi_memory{iSat},polifitOrder);
                    preditPhaseCodeCombi = polyval( p, length(phaseCodeL2Combi_memory{iSat})+1 );
                    if abs(preditPhaseCodeCombi - phaseCodeL2Combi) < ThresholdCycleSlip
                        phaseCodeL2Combi_memory{iSat} = [phaseCodeL2Combi_memory{iSat}, phaseCodeL2Combi];
                        cycleSlipL2Detected(iSat) = 0;
                    else
                        %                             figure; hold on; grid on;
                        %                             plot(phaseCodeL2Combi_memory{iSat},'o')
                        %                             plot( polyval( p, [1: length(phaseCodeL2Combi_memory{iSat})+1] ))
                        %                             plot(  length(phaseCodeL2Combi_memory{iSat})+1, polyval( p,  length(phaseCodeL2Combi_memory{iSat})+1 ), '*')
                        %                             plot( length(phaseCodeL2Combi_memory{iSat})+1,phaseCodeL2Combi,'s','markerfacecolor','r')
                        %                             legend('Previous phase diff', 'Polifit','Prediction','Current phase diff')
                        phaseCodeL2Combi_memory{iSat} = phaseCodeL2Combi;
                        cycleSlipL2Detected(iSat) = 1;
                    end
                end
            end
            
        end
        
        badSatL1 = find(cycleSlipL1Detected==1);
        [meas2deleteL1, meas2deleteIndxL1] = intersect(satPRN_R1,badSatL1);
        C1C_R_aux(meas2deleteIndxL1) = []; L1C_R_aux(meas2deleteIndxL1) = []; D1C_R_aux(meas2deleteIndxL1) = []; SV_pos_R1(meas2deleteIndxL1,:) = []; SV_vel_R1(meas2deleteIndxL1,:) = [];  iono_corr_R1(meas2deleteIndxL1) = []; trop_corr_R1(meas2deleteIndxL1) = []; SV_CLK_R1(meas2deleteIndxL1) = []; SV_SNR_R1(meas2deleteIndxL1) = []; SV_elev_R1(meas2deleteIndxL1) = []; satPRN_R1(meas2deleteIndxL1) = [];  %Number_of_satellites_R1(iGNSS) = length(satPRN_R1);  N_aux_R1 = 1: length(satPRN_R1);
        badSatL2 = find(cycleSlipL2Detected==1);
        [meas2deleteL2, meas2deleteIndxL2] = intersect(satPRN_R2,badSatL2);
        C2C_R_aux(meas2deleteIndxL2) = []; L2C_R_aux(meas2deleteIndxL2) = []; D2C_R_aux(meas2deleteIndxL2) = []; SV_pos_R2(meas2deleteIndxL2,:) = []; SV_vel_R2(meas2deleteIndxL2,:) = [];  iono_corr_R2(meas2deleteIndxL2) = []; trop_corr_R2(meas2deleteIndxL2) = []; SV_CLK_R2(meas2deleteIndxL2) = []; SV_SNR_R2(meas2deleteIndxL2) = []; SV_elev_R2(meas2deleteIndxL2) = []; satPRN_R2(meas2deleteIndxL2) = []; % Number_of_satellites_R2(iGNSS) = length(satPRN_R2);  N_aux_R2 = 1: length(satPRN_R2);
        cycleSlipCounterL1 = cycleSlipCounterL1 + length(meas2deleteIndxL1);
        cycleSlipCounterL2 = cycleSlipCounterL2 + length(meas2deleteIndxL2);
    end
    
    
end

end
