classdef InertialGNSS_EKF < TutorialBasePack.AbstractKalmanFilter
    properties
        F_ % Jacobian of f (derivative w.r.t. state)
        C_ % Jacobian of f (derivative w.r.t. noise)
        H_ % Jacobian of h 
%         sizeState_ % number of elements which are not ambiguities in the state (in case you have position, velocity, attitude, etc etc... or only position, etc)
        errorState_ % error state, composed by: [ pos[1:3], vel[4:6], att[7:9], b_acce[10:12], b_omega[13:15] ]
        
        % Baseline between the GNSS antenna and the IMU 
        baselineIMU2GNSS_ = [0,0,0]';
        angularRate_    = [0,0,0]'; % Angular rate of the vehicle. Useful for applying the lever-arm of the GNSS measurements when the tracked position lays corresponds to the IMU y
        
        % General parameter values
        FL1 = 1575.420 * 10^6;  % GPS [Hz]
        FL2 = 1227.600 * 10^6;  %
        FL5 = 1176.450 * 10^6;  %
        
        FR1_base  = 1602.000 * 10^6;  % GLONASS [Hz]
        FR2_base  = 1246.000 * 10^6;  %
       
        FE1  = 1575.420 * 10^6;     % Galileo [Hz]
        FE5a = 1176.450 * 10^6;     %
        
        V_LIGHT = 299792458;                  % velocity of light in the void [m/s]
        
        % Non-Inertial Noise declaration
        Qvel_   = 0.5;      % [m/s^2]
        Qatt_   = 0.0175;   % [rad/s]
        Qbacc_  = 0.005;    % [m/s^3]
        Qbgyr_  = 0.001;    % [rad/s^2]    
        
        
        
        
    end
    
    methods 
        function obj = InertialGNSS_EKF(varargin)
            
            superclassargs = {[]};
            idxSub = [];         
            for n=1:length(varargin)
                
                if iscell(varargin{n})
                    continue
                end
                
                if strcmp(varargin{n},'-masterObj') && length(varargin) > n
                    superclassargs = [superclassargs varargin(n) varargin(n+1)];
                    idxSub = [idxSub n n+1];
                elseif strcmp(varargin{n},'-plot') && length(varargin) > n
                    superclassargs = [superclassargs varargin(n) varargin(n+1)];
                    idxSub = [idxSub n n+1];
                end
            end
            
            superclassargs(cellfun(@isempty,superclassargs)) = [];
            
            obj@TutorialBasePack.AbstractKalmanFilter(superclassargs{:});
            
            varargin(idxSub) = [];
            varargin(cellfun(@isempty,varargin)) = [];
            
            while ~isempty(varargin)
                if length(varargin) > 1
                    if strcmp(varargin{1},'-dynModel') && ...
                            isa(varargin{2},'function_handle')
                        obj.f_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-measModel') && ...
                            isa(varargin{2},'function_handle')
                        obj.h_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-FJ') && ...
                            isa(varargin{2},'function_handle')
                        obj.F_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-CJ') && ...
                            isa(varargin{2},'function_handle')
                        obj.C_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-HJ') && ...
                            ismatrix(varargin{2})
                        obj.H_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-R') && ...
                            ismatrix(varargin{2})
                        obj.R_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-Q') && ...
                            ismatrix(varargin{2})
                        obj.Q_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-P') && ...
                            ismatrix(varargin{2})
                        obj.P_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-state') && ...
                            isvector(varargin{2})
                        obj.state_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-basePosition') && ...
                            isvector(varargin{2})
                        obj.basePosition_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-sizeState') && ...
                            isvector(varargin{2})
                        obj.sizeState_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-errorState') && ...
                            isvector(varargin{2})
                        obj.errorState_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-labelsTitles') && ...
                            iscell(varargin{2})
                        obj.cTitleLabels_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-labelsYAxes') && ...
                            iscell(varargin{2})
                        obj.cYlabelLabels_ = varargin{2};
                        varargin{2} = [];
                    else
                        myExc = MException('ExtendedKalmanFilter:Constructor', ...
                        sprintf('unknown proporty %s',varargin{1}));
                        throw(myExc);
                    end
                elseif strcmp(varargin{1},'-h') || strcmp(varargin{1},'--help')
                    fprintf('[INFO] help to %s requested:\n\n',obj.getNameFromStack)
                    obj.usage();
                    return
                elseif ~isempty(varargin{1})
                    myExc = MException('ExtendedKalmanFilter:Constructor', ...
                        sprintf('unknown proporty %s',varargin{1}));
                    throw(myExc);
                end
                varargin{1} = [];
                varargin(cellfun(@isempty,varargin)) = [];
            end
            
            if ~isempty(obj.masterObj_)
                addlistener(obj.masterObj_,'observation',...
                    @(src,evnt)obj.eventListener__onObservation(obj,src,evnt));
                
                addlistener(obj.masterObj_,'prediction',...
                    @(src,evnt)obj.eventListener__onPrediction(obj,src,evnt));
            end
            
            obj.init();
        end
        
        function init(obj)
            % init state covariance
                               
            if isempty(obj.P_)
                obj.P_ = 0.5 * eye(length(obj.state_));
            end            
        end
        
        function out = sizeState(obj) % number of elements which are not ambiguities in the state (in case you have position, velocity, attitude, etc etc... or only position, etc)
            out = length(obj.state_);
        end
        
        function resetErrorState(obj)
           obj.errorState_  = zeros(size(obj.errorState_));
        end
        
        function updateTotalState(obj)
           obj.state_ = obj.state_ + obj.errorState_;
        end
        
        function varargout = predictionSola(obj,acc_meas,omega_meas,dt)
                       
            % Calibration of the measurements and other useful operations
            bias_acce   = obj.state_(10:12); bias_gyro = obj.state_(13:15);
            omega_corr  = omega_meas - bias_acce;           % Angular rate corrected using the bias of the gyroscope
            acc_corr    = acc_meas - bias_gyro;             % Acceleration correcting using the bias of the accelerometer
            R_eb        = Euler2DCM(obj.state_(7:9));                       % Rotation from the body to the ECEF frame (using the absolute "total" attitude from the state)
            skew_acc    = skewMatrix(acc_corr);                             % Skew matrix with the corrected measured acceleration
            R_omega     = expm(-0.5*skewMatrix(omega_corr).*dt);
            
            if isempty(obj.angularRate_), obj.angularRate_ = omega_corr; else obj.angularRate_ = (obj.angularRate_ + omega_corr)/2; end % Saving into memory what is the current angular rate of the body --> Necessary for the lever-arm of the GNSS measurements
            
            % Update the error-state -> [ pos[1:3], vel[4:6], att[7:9], b_acce[10:12], b_omega[13:15] ]
            % Disclaimer: this step can be skipped as actually the mean of the error-state will NOT change over time. However, the uncertainty of the system is taken into account 
            obj.errorState_(1:3) = obj.errorState_(1:3) + dt*obj.errorState_(4:6); % Update of the error-position
            obj.errorState_(4:6) = obj.errorState_(4:6) + dt*(-R_eb * skew_acc * obj.errorState_(7:9) - R_eb * obj.errorState_(10:12) ); % Update of the error-velocity
            obj.errorState_(7:9) = skewMatrix(-omega_corr*dt) * obj.errorState_(7:9) - obj.errorState_(13:15)*dt; % Update of the error-attitude
            obj.errorState_(10:12) = obj.errorState_(10:12);                % Update of the error-bias of the accelerometer
            obj.errorState_(13:15) = obj.errorState_(13:15);                % Update of the error-bias of the gyroscope
            
            % System's transition matrix - Jacobian matrix
            F = zeros(15);
            F(1:3,:)    = [eye(3),      dt*eye(3),         zeros(3),               zeros(3),            zeros(3) ];   % Jacobian of the position
            F(4:6,:)    = [zeros(3),    eye(3),            -R_eb*skew_acc*dt,      -R_eb*dt,            zeros(3) ];   % Jacobian of the error-velocity
            F(7:9,:)    = [zeros(3),    zeros(3),           R_omega,               zeros(3),            -eye(3)*dt ]; % Jacobian of the error-attitude
            F(10:12,:)  = [zeros(3),    zeros(3),          zeros(3),               eye(3),              zeros(3) ];   % Jacobian of the error-bias of the accelerometer   
            F(13:15,:)  = [zeros(3),    zeros(3),          zeros(3),               zeros(3),            eye(3)   ];   % Jacobian of the error.bias of the gyroscope
            
            % Jacobian of the noise
            % Q = diag([Noise_acce, Noise_gyro, Noise_b_acce, Noise_b_gyro])
            Jq          = zeros(15,12);
            Jq(1:3,:)   = [zeros(3),        zeros(3),       zeros(3),       zeros(3)];
            Jq(4:6,:)   = [eye(3)*dt,       zeros(3),       zeros(3),       zeros(3)];
            Jq(7:9,:)   = [zeros(3),        eye(3)*dt,      zeros(3),       zeros(3)];
            Jq(10:12,:) = [zeros(3),        zeros(3),       eye(3)*dt,      zeros(3)];
            Jq(13:15,:) = [zeros(3),        zeros(3),       zeros(3),       dt*eye(3)];
                                        
            % Update the covariance matrix
            obj.P_      = F*obj.P_*F' + Jq*obj.Q_*Jq';
                       
            if nargout == 1
                varargout{1} = obj.errorState_;
            else
                varargout = [];
            end
        end
        
        function varargout = predictionNonInertial(obj,dt)
                                  
            % System's transition matrix - Jacobian matrix
            F           = eye(length(obj.state_));
            F(1:3,4:6)  = eye(3)*dt;
            
            obj.state_ = F * obj.state_;
            
            % Jacobian of the noise  -->   Q = diag([Noise_vel, Noise_attitude, Noise_b_acce, Noise_b_gyro, Noise_others])
            Q           = blkdiag(obj.Qvel_*eye(3), obj.Qatt_*eye(3), obj.Qbacc_*eye(3), obj.Qbgyr_*eye(3), zeros(length(obj.state_)-15) );            
            Jq          = zeros([length(obj.state_),length(obj.state_)-3]);
            Jq(1:3,1:3) = eye(3)*dt^2/2;
            Jq(4:end,:) = eye(length(obj.state_)-3)*dt;
                              
            % Update the covariance matrix
            obj.P_      = F*obj.P_*F' + Jq*Q*Jq';
                       
            if nargout == 1
                varargout{1} = obj.errorState_;
            else
                varargout = [];
            end
        end
        
        function varargout = predictionGroves(obj,acc_meas,omega_meas,dt)
                       
            % Error State ->  [delta_pos [1-3], delta_vel[4-6], delta_att[7-9], delta_biasAcce[10-12], delta_biasGyro[13-15], delta_gravity[16-18]]
            % Kinematics             
            omega_corr          = omega_meas - obj.state_(10:12);           % Angular rate corrected using the bias of the gyroscope
            acc_corr            = acc_meas - obj.state_(13:15);             % Acceleration correcting using the bias of the accelerometer
            R_eb        = Euler2DCM(obj.state_(7:9));                       % Rotation from the body to the ECEF frame (using the absolute "total" attitude from the state)
            skew_acc    = skewMatrix(acc_corr);                             % Skew matrix with the corrected measured acceleration
            R_omega     = Euler2DCM(omega_corr*dt)';                        % Skew matrix for the integration of the "corrected" measured angular rate
            
            % Update the state
            obj.errorState_(1:3) = obj.errorState_(1:3) + dt*obj.errorState_(4:6);
            obj.errorState_(4:6) = obj.errorState_(4:6) + (-R_omega * skew_acc * obj.errorState_(7:9) - R_eb * obj.errorState_(10:12) + obj.errorState_(16:18) )*dt;
            obj.errorState_(7:9) = R_omega * obj.errorState_(7:9) - obj.errorState_(13:15)*dt;
            obj.errorState_(10:12) = obj.errorState_(10:12);
            obj.errorState_(13:15) = obj.errorState_(13:15);
            obj.errorState_(16:18) = obj.errorState_(16:18);
            % System's transition matrix
            F = zeros(18);
            F(1:3,:)    = [eye(3),      dt*eye(3),         zeros(3),               zeros(3),            zeros(3),       zeros(3)];   % Jacobian of the position
            F(4:6,:)    = [zeros(3),    eye(3),            -R_eb*skew_acc*dt,      -R_eb*dt,            zeros(3),       eye(3)*dt];
            F(7:9,:)    = [zeros(3),    zeros(3),          R_omega,                zeros(3),            -eye(3)*dt,     zeros(3)];
            F(10:12,:)  = [zeros(3),    zeros(3),          zeros(3),               eye(3),              zeros(3),       zeros(3)]; 
            F(13:15,:)  = [zeros(3),    zeros(3),          zeros(3),               zeros(3),            eye(3),         zeros(3)];
            F(16:18,:)  = [zeros(3),    zeros(3),          zeros(3),               zeros(3),            zeros(3),       eye(3)];
            % Jacobian of the noise
            % Q = diag([Sigma_omega, Sigma_acce, Sigma_biasAcce, Sigma_biasGyro])
            Jq          = zeros(18,12);
            Jq(1:3,:)   = [zeros(3),        zeros(3),       zeros(3),       zeros(3)];
            Jq(4:6,:)   = [eye(3)*dt^2,     zeros(3),       zeros(3),       zeros(3)];
            Jq(7:9,:)   = [zeros(3),        eye(3)*dt^2,    zeros(3),       zeros(3)];
            Jq(10:12,:) = [zeros(3),        zeros(3),       eye(3)*dt,      zeros(3)];
            Jq(13:15,:) = [zeros(3),        zeros(3),       zeros(3),       dt*eye(3)];
            Jq(16:18,:) = [zeros(3),        zeros(3),       zeros(3),       zeros(3)];
                                        
            % Update the covariance matrix
            obj.P_      = F*obj.P_*F' + Jq*obj.Q_*Jq';
%             obj.errorP_      = F*obj.errorP_*F' + Jq*obj.Q_*Jq';
                       
            if nargout == 1
                varargout{1} = obj.errorState_;
            else
                varargout = [];
            end
        end
        
        function varargout = correctionLC(obj, position, position_P)
            % Update the Jacobian matrices
            obj.H_ = zeros(3,15);
            obj.H_(1:3,1:3) = eye(3);
            obj.R_ = position_P;
            % residual covariance
            obj.S_ = obj.H_*obj.P_*obj.H_' + obj.R_;
            % Kalman Gain
            obj.K_ = obj.P_*obj.H_'/obj.S_;
            % compute innovation
            obj.y_ = (position - obj.state_(1:3));
            % update the state estimate.
            obj.state_ = obj.state_ + obj.K_ * obj.y_;
            % update the state covariance.
            obj.P_ = (eye(size(obj.K_,1))-obj.K_*obj.H_)*obj.P_;
            
            if nargout == 1
                varargout{1} = obj.errorState_;
            else
                varargout = [];
            end
        end
        
        function varargout = correctionTC(obj, varargin)
            % Extraction of the data and elimination of the empty
            % information slots ---> Remember that as you need to call for the object, always nargin = the number of inputs you are giving while calling the method +1
            if nargin < 3 || nargin==4 || nargin==5  % there are no ionospheric, tropospheric or SV clockoffset corrections
                disp('Error: Invalid number of arguments')
                return;
            elseif nargin == 3
                C1C          = varargin{1}; C1C(isnan(C1C)) = [];
                SV_Pos       = varargin{2}; SV_Pos(isnan(SV_Pos)) = [];
            elseif nargin < 10  % we are using only the code pseudoranges to update the position
                useDopplerUpdate = 0;
                C1C          = varargin{1}; C1C(isnan(C1C)) = [];
                iono_corr    = varargin{2}; iono_corr(isnan(iono_corr)) = [];
                trop_corr    = varargin{3}; trop_corr(isnan(trop_corr)) = [];
                SV_Clk       = varargin{4}; SV_Clk(isnan(SV_Clk)) = [];
                SV_Pos       = varargin{5}; SV_Pos(isnan(SV_Pos)) = [];
                SV_CN0       = varargin{6}; SV_CN0(isnan(SV_CN0)) = [];
                SV_Elev      = varargin{7}; SV_Elev(isnan(SV_Elev)) = [];
            elseif nargin == 10 % the Doppler shift is also being used to correct the velocity
                useDopplerUpdate = 1;
                C1C          = varargin{1}; C1C(isnan(C1C)) = [];
                iono_corr    = varargin{2}; iono_corr(isnan(iono_corr)) = [];
                trop_corr    = varargin{3}; trop_corr(isnan(trop_corr)) = [];
                SV_Clk       = varargin{4}; SV_Clk(isnan(SV_Clk)) = [];
                SV_Pos       = varargin{5}; SV_Pos(isnan(SV_Pos)) = [];
                SV_CN0       = varargin{6}; SV_CN0(isnan(SV_CN0)) = [];
                SV_Elev      = varargin{7}; SV_Elev(isnan(SV_Elev)) = [];
                D1C          = varargin{8}; D1C(isnan(D1C)) = [];
                SV_Vel       = varargin{9}; SV_Vel(isnan(SV_Vel)) = [];
            end
            
            % Check the number of measurements. As long as we have >1 the state correction is more or less effective
            NumberOfSatellite = length(C1C);
            C1C = C1C - iono_corr - trop_corr + SV_Clk;
            if NumberOfSatellite < 4, 
                disp('Not enough satellites'); 
                varargout{1} = [];
                return;
            end
            
            % External SPP to track the receiver clock offset, avoiding convergence issues whenever the receiver internally resets its clock offset
            [~, CLKOffset_aux, ~, ~, ~, ~, ~] = Snapshot_Position ([0,0,0]', ones(NumberOfSatellite,1), C1C, SV_Pos(:,1), SV_Pos(:,2), SV_Pos(:,3), 20, 10^(-5), []);
            if abs(CLKOffset_aux - obj.state_(16)) > 10^4,      obj.state_(16) =  CLKOffset_aux;    end % Overwrite the current estimated clock offset if a clock offset reset is perceived
                     
            % Noise model for the measurements
            R_C1C       = diag([ 100*(0.13 + 0.56*exp( -SV_Elev/0.1745 ) ).^2 ]);
            if useDopplerUpdate, R_D1C       = diag([ (0.013 + 0.056*exp( -SV_Elev/0.1745 ) ).^2 ]); end
            
            
            Rec_ClkOffset = obj.state_(16);
            Rec_ClkOffsetRate = obj.state_(17);
            
            % Lever-arm from IMU -> GNSS antenna
            positionGNSS = obj.state_(1:3) + Euler2DCM(obj.state_(7:9))*obj.baselineIMU2GNSS_; % Lever-arm to get the position of the antenna from the position of the IMU + the baseline between both  REF: [Groves Book] Pag. 393
            if useDopplerUpdate, 
                GNSS_vel = obj.state_(4:6) + Euler2DCM(obj.state_(7:9))* cross(obj.angularRate_, obj.baselineIMU2GNSS_) ; 
            end % Lever-arm to get the velocity of the antenna from the position of the IMU + the baseline between both REF: [Groves Book] Pag. 393 
            
            if useDopplerUpdate
                % Initialization of the variables
                h_range = zeros(NumberOfSatellite,1);
                h_doppler = zeros(NumberOfSatellite,1);
                z_range = zeros(NumberOfSatellite,1);
                z_doppler = zeros(NumberOfSatellite,1);
                H_range = zeros(NumberOfSatellite, obj.sizeState());
                H_doppler = zeros(NumberOfSatellite, obj.sizeState());
                for iSV = 1:NumberOfSatellite
                    % Measurements
                    z_range(iSV) = C1C(iSV); % Measurements for the code measurements
                    z_doppler(iSV) = obj.V_LIGHT*D1C(iSV)/obj.FL1; % Measurements for the Doppler shift
                    
                    % Geometric distance between the satellite position and the target receiver predicted position
                    geoRange(iSV) = norm(SV_Pos(iSV,:) - positionGNSS');
                    uVector(iSV,:) = (SV_Pos(iSV,:) - positionGNSS')./geoRange(iSV);
                    
                    % Observation model
                    h_range(iSV) = geoRange(iSV) + Rec_ClkOffset;
                    h_doppler(iSV) = uVector(iSV,:)*(SV_Vel(iSV,:)-GNSS_vel')' + Rec_ClkOffsetRate;
                    
                    % Vectors with the differences of velocity and position of satellite and receiver
                    Pos_diff = [ SV_Pos(iSV,1)-positionGNSS(1), SV_Pos(iSV,2)-positionGNSS(2), SV_Pos(iSV,3)-positionGNSS(3) ];
                    Vel_diff = [ SV_Vel(iSV,1)-GNSS_vel(1), SV_Vel(iSV,2)-GNSS_vel(2), SV_Vel(iSV,3)-GNSS_vel(3) ];
                    
                    % Jacobian matrix of the observation model
                    H_range(iSV,1:3) = -uVector(iSV,:);
                    H_range(iSV,16) = 1;
                    
                    H_doppler(iSV,4:6) = -uVector(iSV,:);
                    H_doppler(iSV,17) = 1;
                    
                    % When estimating the velocity, there is a weak dependance on the position, therefore generally that is eliminated
                    %                 H_doppler(iSV,1) = -Vel_diff(1)/geoRange(iSV) + Vel_diff(1)*Pos_diff(1)^2/geoRange(iSV)^3 + Pos_diff(2:3)*Vel_diff(2:3)'     * ( Pos_diff(1)/geoRange(iSV)^3 );
                    %                 H_doppler(iSV,2) = -Vel_diff(2)/geoRange(iSV) + Vel_diff(2)*Pos_diff(2)^2/geoRange(iSV)^3 + Pos_diff([1,3])*Vel_diff([1,3])' * ( Pos_diff(2)/geoRange(iSV)^3 );
                    %                 H_doppler(iSV,3) = -Vel_diff(3)/geoRange(iSV) + Vel_diff(3)*Pos_diff(3)^2/geoRange(iSV)^3 + Pos_diff(1:2)*Vel_diff(1:2)'     * ( Pos_diff(3)/geoRange(iSV)^3 );
                end
                obj.H_ = [H_range;H_doppler];
                obj.R_   = blkdiag( R_C1C, R_D1C );
                
            else
                % Initialization of the variables
                h_range = zeros(NumberOfSatellite,1);
                z_range = zeros(NumberOfSatellite,1);
                H_range = zeros(NumberOfSatellite, obj.sizeState());
                for iSV = 1:NumberOfSatellite
                    % Measurements
                    z_range(iSV) = C1C(iSV); % Measurements for the code measurements                    
                    % Geometric distance between the satellite position and the target receiver predicted position
                    geoRange(iSV) = norm(SV_Pos(iSV,:) - positionGNSS');
                    uVector(iSV,:) = (SV_Pos(iSV,:) - positionGNSS')./geoRange(iSV);
                    % Observation model
                    h_range(iSV) = geoRange(iSV) + Rec_ClkOffset;
                    % Jacobian matrix of the observation model
                    H_range(iSV,1:3) = uVector(iSV,:);
                    H_range(iSV,16) = 1;
                end
                obj.H_ = [H_range];
                obj.R_   = blkdiag( R_C1C);
                
            end
            
            % Correction Step
            obj.S_ = obj.H_*obj.P_*obj.H_' + obj.R_;           % Innovation covariance
            obj.K_ = obj.P_*obj.H_'/obj.S_;                    % Kalman Gain
            y_range = (z_range - h_range);                            % Innovation for the pseudorange measurements
            if useDopplerUpdate, y_doppler = (z_doppler - h_doppler);   end                         % Innovation for the Doppler shift measurements
            if useDopplerUpdate, 
                obj.y_  = [y_range; y_doppler];
            else
                obj.y_  = [y_range];
            end
            obj.state_          = obj.state_ + obj.K_ * obj.y_;             % update the mean of the state
            obj.P_              = (eye(size(obj.P_,1)) - obj.K_ * obj.H_) * obj.P_;  % update the state covariance.
            
            if nargout == 1
                varargout{1} = obj.errorState_;
            else
                varargout = [];
            end
        end
        
        
        function varargout = correction(obj,sObs)       
            % residual covariance
            obj.S_ = obj.H_*obj.P_*obj.H_' + obj.R_;
            
            % Kalman Gain
            obj.K_ = obj.P_*obj.H_'/obj.S_;
                        
            % compute innovation
            obj.y_ = (sObs - obj.h_(obj.state_));
            % update the state estimate.
            obj.state_ = obj.state_ + obj.K_ * obj.y_;
            % update the state covariance.
            obj.P_ = (eye(size(obj.K_,1))-obj.K_*obj.H_)*obj.P_; 
            
            if nargout == 1
                varargout{1} = obj.state_;
            else
                varargout = [];
            end           
        end
        
        function varargout = prediction(obj)       
            
            if nargout == 1
                varargout{1} = obj.state_;
            else
                varargout = [];
            end           
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% LISTENER IMPLEMENTATION %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Static)
        function eventListener__onObservation(obj,~,evt)
            %%% event is triggered on new observation
            
            if obj.status_ ~= TutorialBasePack.BaseStates.Idle && ...
                    obj.status_ ~= TutorialBasePack.BaseStates.Finished
                fprintf('[INFO] %s is currently busy, waiting to continue at %s \n',class(obj),obj.getNameFromStack);
                waitfor(obj.getID(),'status_',TutorialBasePack.BaseStates.Finished)
            end
            
            %%% set internal status to active
            obj.status_ = TutorialBasePack.BaseStates.Active;
            
            %%% update state with new measurement
            obj.correction(evt.obs);

            %%% set internal state back to finished
            obj.status_ = TutorialBasePack.BaseStates.Finished;
            
            if obj.plot_
                obj.plotIntermediates('update',evt.dt)
            end
        end
        
        function eventListener__onPrediction(obj,~,evt)
            %%% event is triggered on new observation
            
            fprintf('event listener activated at %s\n',...
                datestr(now,'dd-mmm-yyyy HH:MM:SS.FFF'));
            
            if obj.status_ ~= TutorialBasePack.BaseStates.Idle && ...
                    obj.status_ ~= TutorialBasePack.BaseStates.Finished
                fprintf('[INFO] %s is currently busy, waiting to continue at %s \n',class(obj),obj.getNameFromStack);
                waitfor(obj.getID(),'status_',TutorialBasePack.BaseStates.Finished)
            end
            
            %%% set internal status to active
            obj.status_ = TutorialBasePack.BaseStates.Active;
            
            %%% update state with new measurement
            obj.prediction(evt.dt);
                        
            %%% set internal state back to finished
            obj.status_ = TutorialBasePack.BaseStates.Finished;
            
            if obj.plot_
                obj.plotIntermediates('prediction',evt.dt)
            end
        end   
        
    end
end