classdef RTK_Release < TutorialBasePack.AbstractKalmanFilter
    properties
        F_ % Jacobian of f (derivative w.r.t. state)
        C_ % Jacobian of f (derivative w.r.t. noise)
        H_ % Jacobian of h 
        satellitesLOS_ = []; % vector with the satellites PRN in LOS at the moment
        satRef_ = []; % reference satellite for the doble difference approach for RTK
        basePosition_ % position of the base or reference station, from which the position is known
        sizeState_ = 6; % number of elements which are not ambiguities in the state (in case you have position, velocity, attitude, etc etc... or only position, etc)
        WavelengthL1_ = 0.19029367;
        WavelengthL2_ = 0.24421021;
        sizeSystemNoise_ = 3;
        DDAmb_ 
        SDAmb_ 
        D_
        QRN_
        QNR_
        QN_
        initialAmbCovariance_ = 10^2;
        leverArmIMU2GNSS_ = [0 0 0]';
        eulerAngles_ = [];
        Qamb_ = (10^-16);
%         Qamb_ = 0;
        Qvel_     = 0.1^2; %(  m/s^2 )^2
        Qgyr_  = 1000*(1/3600*1/3600)*0.005*(pi/180)^2*diag([0.03*0.03 0.03*0.03 0.03*0.03]);    % Accelerometer measurement noise covariance (interpreted as a control noise)
        Qacc_  = 9.81*9.81*0.005*(10e-12)*diag([320*320 320*320 320*320]);            % Accelerometer measurement noise covariance (interpreted as a control noise)
        QbiasGyr_ = 0.001*(1/3600*1/3600)*1/0.005*(pi/180)^2*diag([0.01*0.01 0.01*0.01 0.01*0.01]);% Process noise for estimated gyroscope bias (should be scaled with sqrt(dt))
        QbiasAcc_  = 0.001*9.81*9.81*1/0.005*(10e-12)*diag([50*50 50*50 50*50]);            % Accelerometer measurement noise covariance (interpreted as a control noise)
        resetDDCrossVarianceP_ = 0;
        debugMode_ = 0;
        mode_
        acceleration_ = zeros(3,1);
        angularRate_ = zeros(3,1);
    end
    
    methods 
        function obj = RTK_Release(varargin)
            
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
                    elseif strcmp(varargin{1},'-mode') && ...
                            isa(varargin{2},'function_handle')
                        obj.mode_ = varargin{2};
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
                    elseif strcmp(varargin{1},'-leverArmIMU2GNSS') && ...
                            isvector(varargin{2})
                        obj.leverArmIMU2GNSS_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-sizeState') && ...
                            isvector(varargin{2})
                        obj.sizeState_ = varargin{2};
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
        
        %% Getter for the estimated position 
        function position = getTrackedPosition(obj)
            if obj.sizeState_ == 16 % In case we are using a INS system in which we are using orientation quaternion
                position = obj.state_(1:3) +  quat_rotate( obj.state_(7:10), obj.leverArmIMU2GNSS_ );
            else
                position = obj.state_(1:3);
            end
        end
        
        %% Prediction Step for the non-inertial system
        function varargout = prediction(obj,dt)           
            % Update the mean of the state (this can also be done using matrix multiplication)
            obj.state_ = obj.f_(dt,obj.state_);
            
            % Set the noise for the prediction
            obj.setSystemNoise( dt^2*blkdiag(obj.Qvel_*eye(3), obj.Qamb_*eye(length(obj.satellitesLOS_))) );
            
            % prediction of state covariance
            Jf = obj.F_(dt,obj.state_);
            Jq = obj.C_(dt,obj.state_);
            obj.P_ = Jf * obj.P_ * Jf' + Jq * obj.Q_ * Jq';
           
            if nargout == 1
                varargout{1} = obj.state_;
            else
                varargout = [];
            end
        end
        
        %% Original Correction Function
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
        
        %% Prediction Step for INS
        function varargout = predictionInertial(obj, acceMeas, omegaMeas, dt)
            % The state is constructed as following:
            % x = [position, velocity, attitudeQuaternion, biasAccelerometer, biasGyroscope, others...]
            
            % Saving the components of the state into variables with the corresponding name [not needed, but good to directly refer to those on the rest of the equations]
            posOld = obj.state_(1:3);
            velOld = obj.state_(4:6);
            quatOld = obj.state_(7:10);
            biasAcceOld = obj.state_(11:13);
            biasGyroOld = obj.state_(14:16);
            sizeNoINS = length(obj.state_) - obj.sizeState_;               % Number of elements in the state which are not related to Inertial navigation (GNSS clock offset, phase ambiguities and others)
            
            % Forcing sensor data to be a column vector
            acceMeas = reshape(acceMeas,3,1);
            omegaMeas = reshape(omegaMeas,3,1);
            
            % Correcting the sensors' bias
            acceCorr = acceMeas - obj.state_(11:13);
            omegaCorr = omegaMeas - obj.state_(14:16);
            gravity = xyz2grav(obj.state_(1),obj.state_(2),obj.state_(3));  % ECEF gravity using a model for it
            
            % Rotation matrix Body -> ECEF
            R_b2e = quat_2_DCM_nb(quatOld); 	% Rotation matrix associated to the previous quaternion
            
            
            %%%%%%%%%%%
            %%%%%%%%%%%
            theta = norm(0.5*dt*omegaCorr);
            Omega = [ 0, omegaCorr'; -omegaCorr, skewMatrix(omegaCorr) ] ;
            
            omega_x = [ 0 1 0 0; -1 0 0 0; 0 0 0 -1; 0 0 1 0  ];
            omega_y = [ 0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0 ];
            omega_z = [ 0 0 0 1; 0 0 -1 0; 0 1 0 0; -1 0 0 0 ];
            
            term1_x = eye(4) * sin( theta ) * biasGyroOld(1)/theta;
            term2_x = 1/2 * dt * omega_x * sin (theta)/theta;
            term3_x = 1/2*dt^2 * Omega * ( theta*cos(theta)-sin(theta) )/theta^3;
            term_x = -(term1_x + term2_x + term3_x);
            
            term1_y = eye(4) * sin( theta ) * biasGyroOld(2)/theta;
            term2_y = 1/2 * dt * omega_y * sin (theta)/theta;
            term3_y = 1/2*dt^2 * Omega * ( theta*cos(theta)-sin(theta) )/theta^3;
            term_y = -(term1_y + term2_y + term3_y);
            
            term1_z = eye(4) * sin( theta ) * biasGyroOld(3)/theta;
            term2_z = 1/2 * dt * omega_z * sin (theta)/theta;
            term3_z = 1/2*dt^2 * Omega * ( theta*cos(theta)-sin(theta) )/theta^3;
            term_z = -(term1_z + term2_z + term3_z);
            
            F_q_bw = [ term_x*quatOld, term_y*quatOld, term_z*quatOld ];
            %%%%%%%%%%%
            %%%%%%%%%%%
            
            % Derivative of the measured acceleration set into the ECEF frame ( q\otimes acceCorr \otimes q^* ) with respect to the quaternion. Reference: [Joan Sola 2017] Quaternion Dynamics for Error State KF, page 44  
            qw    = quatOld(1);                                 % Real part of the quaternion
            qu    = quatOld(2:4);                               % Imaginary part of the quaternion
            DiffAcceECEFdiffQuat = 2 * [ qw*acceMeas + cross(qu, acceMeas), qu'*acceMeas*eye(3) + qu*acceMeas' - acceMeas*qu' - qw*skewMatrix(acceMeas)  ];
            DiffAcceECEFdiffQuat = 1 * [ qw*acceMeas + cross(qu, acceMeas), qu'*acceMeas*eye(3) + qu*acceMeas' - acceMeas*qu' - qw*skewMatrix(acceMeas)  ];
            Fquat = DiffAcceECEFdiffQuat;
            
            % The update of the quaternion: q\otimes \delta q, where \delta q comes from the zeroth order integration of the angular rate, can also be expressed in a matrix multiplication form:  [\delta q]_R q. Look at quaternion multiplication properties,  [Joan Sola 2017] Page 7
            deltaQuat = [ cos(norm(omegaCorr)*dt/2);  omegaCorr/norm(omegaCorr)*sin(norm(omegaCorr)*dt/2) ]; % Matrix for q \otimes \delta q = [delta q]_R q . For that [delta q]_R
            deltaQuatMatrixR = deltaQuat(1)*eye(4) + [0, -deltaQuat(2:4)'; deltaQuat(2:4) -skewMatrix(deltaQuat(2:4))]; %  
             
             % Jacobian of the process model
%              Fx = [  eye(3),         dt*eye(3),      dt^2/2*Fquat,      -dt^2/2*R_b2e,     zeros(3),        zeros(3,sizeNoINS); ...
%                      zeros(3),       eye(3),         dt*Fquat,           -dt*R_b2e,         zeros(3),        zeros(3,sizeNoINS); ...
%                      zeros(4,3),     zeros(4,3),     deltaQuatMatrixR,   zeros(4,3),    zeros(4,3),      zeros(4,sizeNoINS); ...
%                      zeros(3),       zeros(3),       zeros(3,4),         eye(3),        zeros(3),        zeros(3,sizeNoINS); ...
%                      zeros(3),       zeros(3),       zeros(3,4),         zeros(3),      eye(3),          zeros(3,sizeNoINS); ...
%                      zeros(sizeNoINS,16),                                                                eye(sizeNoINS)      ];
                 
                 Fx = [  eye(3),         dt*eye(3),      dt^2/2*Fquat,      -dt^2/2*R_b2e,     zeros(3),        zeros(3,sizeNoINS); ...
                     zeros(3),       eye(3),         dt*Fquat,           -dt*R_b2e,         zeros(3),        zeros(3,sizeNoINS); ...
                     zeros(4,3),     zeros(4,3),     deltaQuatMatrixR,   zeros(4,3),    F_q_bw,      zeros(4,sizeNoINS); ...
                     zeros(3),       zeros(3),       zeros(3,4),         eye(3),        zeros(3),        zeros(3,sizeNoINS); ...
                     zeros(3),       zeros(3),       zeros(3,4),         zeros(3),      eye(3),          zeros(3,sizeNoINS); ...
                     zeros(sizeNoINS,16),                                                                eye(sizeNoINS)      ];
                 
            B  = [  dt^2/2 * eye(3); ...
                    dt * eye(3);     ...
                    zeros(4,3);      ...
                    zeros(3);        ...
                    zeros(3);        ...
                    zeros(sizeNoINS,3) ];
            
            u  = [ gravity ];
    
            % Jacobian of the process noise
%             Fq = [  dt^3/2*R_b2e,       zeros(3),      -dt^3/2*R_b2e,      zeros(3),       zeros(3,sizeNoINS); ...
%                      dt^3*R_b2e,         zeros(3),      -dt^2*R_b2e,        -zeros(3),      zeros(3,sizeNoINS); ...
%                      zeros(4,3),     zeros(4,3),     zeros(4,3),    zeros(4,3),     zeros(4,sizeNoINS); ...
%                      zeros(3),       zeros(3),       dt*eye(3),     zeros(3),       zeros(3,sizeNoINS); ...
%                      zeros(3),       zeros(3),       zeros(3),      dt*eye(3),      zeros(3,sizeNoINS); ...
%                      zeros(sizeNoINS,12),                                           dt*eye(sizeNoINS)      ];
                 
                 Fq = [  dt^3/2*R_b2e,       zeros(3),      -dt^3/2*R_b2e,      zeros(3),       zeros(3,sizeNoINS); ...
                     dt^3*R_b2e,         zeros(3),      -dt^2*R_b2e,        -zeros(3),      zeros(3,sizeNoINS); ...
                     zeros(4,3),     dt*F_q_bw,     zeros(4,3),    -dt*F_q_bw,     zeros(4,sizeNoINS); ...
                     zeros(3),       zeros(3),       dt*eye(3),     zeros(3),       zeros(3,sizeNoINS); ...
                     zeros(3),       zeros(3),       zeros(3),      dt*eye(3),      zeros(3,sizeNoINS); ...
                     zeros(sizeNoINS,12),                                           dt*eye(sizeNoINS)      ];
             
              newStateMatrix = Fx*obj.state_ + B*u;   
                 
              % Update of the state using directly the nonlinear equations
             accelerationECEF = R_b2e * (acceMeas - obj.state_(11:13) ) + gravity;
             newPosition = obj.state_(1:3) + dt*obj.state_(4:6) + dt^2/2 * accelerationECEF;
             newVelocity = obj.state_(4:6) + dt * accelerationECEF;
             newQuaternion = quat_mult( quatOld, deltaQuat );
             newBiasAcce = biasAcceOld;
             newBiasGyro = biasGyroOld;
             newStateEquation = [newPosition; newVelocity; newQuaternion; newBiasAcce; newBiasGyro; obj.state_(obj.sizeState_+1:end)];
             
             % Mean of the state update
             obj.state_ = newStateMatrix;
             obj.state_ = newStateEquation;
             obj.state_(7:10) = quat_normalize(obj.state_(7:10)); % Normalization of the quaternion, just to assure that it represents a proper rotation. 
             obj.eulerAngles_ = quat2euler(obj.state_(7:10));
                          
             % Set the noise for the prediction
             obj.setSystemNoise( dt^2*blkdiag(obj.Qacc_*eye(3), obj.Qgyr_*eye(3), obj.QbiasAcc_*eye(3), obj.QbiasGyr_*eye(3), obj.Qamb_*eye(length(obj.satellitesLOS_))) );
             
             % Covariance Update
             newPmatrix = Fx * obj.P_ * Fx' + Fq * obj.Q_ * Fq';
             obj.P_     = newPmatrix; 
             
            if nargout == 1
                varargout{1} = obj.state_;
            else
                varargout = [];
            end
        end
        
        %% Prediction Step for INS
        function varargout = NEWpredictionInertial(obj, acceMeas, omegaMeas, dt)
            % The state is constructed as following:
            % x = [position, velocity, attitudeQuaternion, biasAccelerometer, biasGyroscope, others...]
            
            % Forcing sensor data to be a column vector
            acceMeas = reshape(acceMeas,3,1);
            omegaMeas = reshape(omegaMeas,3,1);
            
            % Saving the components of the state into variables with the corresponding name [not needed, but good to directly refer to those on the rest of the equations]
            posOld = obj.state_(1:3);
            velOld = obj.state_(4:6);
            quatOld = obj.state_(7:10);
            qw = quatOld(1); % Real part of the quaternion
            qu = quatOld(2:4); % Imaginary part of the quaternion ---> to be used in Fq =  \diff quat / \diff quat
            biasAcceOld = obj.state_(11:13);
            biasGyroOld = obj.state_(14:16);
            sizeNoINS = length(obj.state_) - obj.sizeState_;               % Number of elements in the state which are not related to Inertial navigation (GNSS clock offset, phase ambiguities and others)
            C_b2e = quat_2_DCM(quatOld);    % Rotation from Body to ECEF
            
            % Earth rotation, gravity and centrifugal and Coriolis accelerations
            EarthRate = 7.292115*10^(-5);
            EarthRateVector = [0 0 EarthRate]';
            gravity = xyz2grav(obj.state_(1),obj.state_(2),obj.state_(3));  % ECEF gravity using a model for it
            centrifugalAcc = - cross(  EarthRateVector, cross(EarthRateVector, posOld)  );
            CoriolisAcc = - 2 * cross( EarthRateVector, velOld );
            
            % Correction of the angular rate
            omegaCorr = 0.5*(omegaMeas + obj.angularRate_) - biasGyroOld - C_b2e' * EarthRateVector ;
            % Correction of the acceleration
            acceCorrAux = C_b2e * acceMeas - biasAcceOld + gravity + centrifugalAcc + CoriolisAcc;
            acceCorr = 1/2 * ( acceCorrAux + obj.acceleration_ );
           
            % Saving the measured angular rate and acceleration to smooth out the following inertial measurements
            obj.angularRate_ = omegaMeas;
            obj.acceleration_ = acceCorrAux;
            
            
            % Non linear equations for the prediction model
            posNew = posOld + dt*velOld + 0.5*dt^2*acceCorr; % Update the position 
            velNew = velOld + dt * acceCorr;          % Update the velocity
            % Quaternion update
            theta = norm(omegaCorr)*dt*0.5;
            uVector = omegaCorr / norm(omegaCorr);
            deltaQuat = [ cos(theta); uVector * sin(theta) ];
            deltaQuatMatrix = quat_productMatrix(deltaQuat, 'right');
            quatNew = deltaQuatMatrix * quatOld;
            % Update the biases of the inertial sensors
            biasAcceNew = biasAcceOld;
            biasGyroNew = biasGyroOld;
            
            % Jacobian of the system dynamics
            Fx = zeros(size(obj.state_));
            
            % Jacobian of the position
            F_pos_pos = eye(3) + 0.5*dt*dt * EarthRate^2 * [ 1,0,0; 0,1,0; 0,0,0 ];
            F_pos_vel  = dt*eye(3) + dt*dt * EarthRate * [ 0,1,0; -1,0,0; 0,0,0 ];
            F_pos_quat = 0.5*dt*dt * [ qw*acceMeas + cross(qu, acceMeas), qu'*acceMeas*eye(3) + qu*acceMeas' - acceMeas*qu' - qw*skewMatrix(acceMeas)  ];
            F_pos_bAcc = -0.5*dt*dt* C_b2e; 
            F_pos_bGyr = zeros(3);
            F_pos_other = zeros(3, sizeNoINS);
            
            % Jacobian of the velocity
            F_vel_pos = 0.5*dt*EarthRate^2*[1,0,0; 0,1,0; 0,0,0];
            F_vel_vel  = eye(3) + dt * EarthRate * [ 0,1,0; -1,0,0; 0,0,0 ];
            F_vel_quat = dt * [ qw*acceMeas + cross(qu, acceMeas), qu'*acceMeas*eye(3) + qu*acceMeas' - acceMeas*qu' - qw*skewMatrix(acceMeas)  ];
            F_vel_bAcc = -dt * C_b2e; 
            F_vel_bGyr = zeros(3);
            F_vel_other = zeros(3, sizeNoINS);
            
            % Jacobian of the quaternion
            F_quat_pos = zeros(4,3);
            F_quat_vel = zeros(4,3);
            F_quat_quat = deltaQuatMatrix;
            F_quat_bAcc = zeros(4,3);
            % Jacobian of the quaternion with respect to the bias of the gyroscope
            theta = norm(0.5*dt*omegaCorr);
            Omega = zeros(4);
            Omega = [ 0, omegaCorr'; -omegaCorr, skewMatrix(omegaCorr) ] ;
            omega_x = [ 0 1 0 0; -1 0 0 0; 0 0 0 -1; 0 0 1 0  ];
            omega_y = [ 0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0 ];
            omega_z = [ 0 0 0 1; 0 0 -1 0; 0 1 0 0; -1 0 0 0 ];            
            term1_x = eye(4) * sin( theta ) * biasGyroOld(1)/theta;
            term2_x = 1/2 * dt * omega_x * sin (theta)/theta;
            term3_x = 1/2*dt^2 * Omega * ( theta*cos(theta)-sin(theta) )/theta^3;
            term_x = -(term1_x + term2_x + term3_x);
            term1_y = eye(4) * sin( theta ) * biasGyroOld(2)/theta;
            term2_y = 1/2 * dt * omega_y * sin (theta)/theta;
            term3_y = 1/2*dt^2 * Omega * ( theta*cos(theta)-sin(theta) )/theta^3;
            term_y = -(term1_y + term2_y + term3_y);            
            term1_z = eye(4) * sin( theta ) * biasGyroOld(3)/theta;
            term2_z = 1/2 * dt * omega_z * sin (theta)/theta;
            term3_z = 1/2*dt^2 * Omega * ( theta*cos(theta)-sin(theta) )/theta^3;
            term_z = -(term1_z + term2_z + term3_z);
            F_quat_bGyr = [ term_x*quatOld, term_y*quatOld, term_z*quatOld ];
            F_quat_other = zeros(4,sizeNoINS);
             % Jacobian of the accelerometer bias
             F_bAcc = [ zeros(3, 10), eye(3), zeros(3), zeros(3, sizeNoINS) ];
             % Jacobian of the gyroscope bias
             F_bGyr = [ zeros(3, 13), eye(3),  zeros(3, sizeNoINS) ];
             % Jacobian of the other elements (phase ambiguities, clockoffsets, others...)
             F_other = [ zeros(sizeNoINS, obj.sizeState_), eye(sizeNoINS) ];
             
            Fx = [ F_pos_pos, F_pos_vel, F_pos_quat, F_pos_bAcc, F_pos_bGyr, F_pos_other; ...
                      F_vel_pos, F_vel_vel, F_vel_quat, F_vel_bAcc, F_vel_bGyr, F_vel_other; ...
                      F_quat_pos, F_quat_vel, F_quat_quat, F_quat_bAcc, F_quat_bGyr, F_quat_other;
                      F_bAcc;...
                      F_bGyr;...
                      F_other ];
            
            % Control input and its relationship to the system dynamics     
                 u  = [ gravity + centrifugalAcc + CoriolisAcc ];
                 B  = [  dt^2/2 * eye(3); ...
                    dt * eye(3);     ...
                    zeros(4,3);      ...
                    zeros(3);        ...
                    zeros(3);        ...
                    zeros(sizeNoINS,3) ];
            
            % Jacobian of the process noise                
            Fq_pos_acceNoise = 0.5*dt*dt * C_b2e;
            Fq_pos_gyroNoise = zeros(3);
            Fq_pos_biasAcceNoise = 0.5*dt*dt * C_b2e;
            Fq_pos_biasGyroNoise = zeros(3);
            Fq_pos_other = zeros(3,sizeNoINS);
            
            Fq_vel_acceNoise = dt * C_b2e;
            Fq_vel_gyroNoise = zeros(3);
            Fq_vel_biasAcceNoise = dt * dt * C_b2e;
            Fq_vel_biasGyroNoise = zeros(3);
            Fq_vel_other = zeros(3,sizeNoINS);
            
            Fq_quat_acceNoise = zeros(4,3);
            Fq_quat_gyroNoise = F_quat_bGyr;
            Fq_quat_biasAcceNoise = zeros(4,3);
            Fq_quat_biasGyroNoise = F_quat_bGyr;
            Fq_quat_other = zeros(4,sizeNoINS);
            
            Fq_bAcc_acceNoise = zeros(3);
            Fq_bAcc_gyroNoise = zeros(3);
            Fq_bAcc_biasAcceNoise = eye(3);
            Fq_bAcc_biasGyroNoise = zeros(3);
            Fq_bAcc_other = zeros(3,sizeNoINS);
            
            Fq_bGyr_acceNoise = zeros(3);
            Fq_bGyr_gyroNoise = zeros(3);
            Fq_bGyr_biasAcceNoise = zeros(3);
            Fq_bGyr_biasGyroNoise = eye(3);
            Fq_bGyr_other = zeros(3,sizeNoINS);
            
            Fq_other_inertialNoise = zeros(sizeNoINS, 12);
            Fq_other_other = eye(sizeNoINS);
            
             Fq = [ Fq_pos_acceNoise,       Fq_pos_gyroNoise,       Fq_pos_biasAcceNoise,       Fq_pos_biasGyroNoise,       Fq_pos_other; ...
                        Fq_vel_acceNoise,           Fq_vel_gyroNoise,       Fq_vel_biasAcceNoise,       Fq_vel_biasGyroNoise,       Fq_vel_other; ...
                        Fq_quat_acceNoise,      Fq_quat_gyroNoise,      Fq_quat_biasAcceNoise,      Fq_quat_biasGyroNoise,      Fq_quat_other; ...
                        Fq_bAcc_acceNoise,      Fq_bAcc_gyroNoise,      Fq_bAcc_biasAcceNoise,      Fq_bAcc_biasGyroNoise,      Fq_bAcc_other;...
                        Fq_bGyr_acceNoise,      Fq_bGyr_gyroNoise,      Fq_bGyr_biasAcceNoise,      Fq_bGyr_biasGyroNoise,      Fq_bGyr_other;...
                        Fq_other_inertialNoise, Fq_other_other                        ];
             
              newStateMatrix = Fx*obj.state_ + B*u;   % One can compare that the update using the non-linear equations and the update using the matrix form are equivalent
                 
             % Mean of the state update
             obj.state_(1:obj.sizeState_) = [ posNew; velNew; quatNew; biasAcceNew; biasGyroNew ];
             obj.state_(7:10) = quat_normalize(obj.state_(7:10)); % Normalization of the quaternion, just to assure that it represents a proper rotation. 
             
             % Covariance Update
             obj.setSystemNoise( diag(  dt * [ diag(obj.Qacc_); diag(obj.Qgyr_); diag(obj.QbiasAcc_); diag(obj.QbiasGyr_); obj.Qamb_*ones(length(obj.satellitesLOS_),1) ] ) ); % Set the noise for the prediction
             newPmatrix = Fx * obj.P_ * Fx' + Fq * obj.Q_ * Fq';
             obj.P_     = newPmatrix; 
             
             % Estimation of the Euler Angles  
             C_b2e = quat_2_DCM(quatNew);
             C_e2n = ECEF_2_NED_rotation(posNew);
             C_b2n = C_e2n * C_b2e;
             obj.eulerAngles_ = DCM_2_euler(C_b2n)*180/pi;
             
            if nargout == 1
                varargout{1} = obj.state_;
            else
                varargout = [];
            end
        end
        
        %% Prediction Step for INS
        function varargout = LASTpredictionInertial(obj, acceMeas, omegaMeas, dt)
            % The state is constructed as following:
            % x = [position, velocity, attitudeQuaternion, biasAccelerometer, biasGyroscope, others...]
            
            % Forcing sensor data to be a column vector
            acceMeas = reshape(acceMeas,3,1);
            omegaMeas = reshape(omegaMeas,3,1);
            
            % Saving the components of the state into variables with the corresponding name [not needed, but good to directly refer to those on the rest of the equations]
            posOld = obj.state_(1:3);
            velOld = obj.state_(4:6);
            quatOld = obj.state_(7:10);
            qw = quatOld(1); % Real part of the quaternion
            qu = quatOld(2:4); % Imaginary part of the quaternion 
            biasAcceOld = obj.state_(11:13);
            biasGyroOld = obj.state_(14:16);
            sizeNoINS = length(obj.state_) - obj.sizeState_;               % Number of elements in the state which are not related to Inertial navigation (GNSS clock offset, phase ambiguities and others)
            C_b2e = quat_2_DCM(quatOld);    % Rotation from Body to ECEF
            
            % Earth rotation, gravity and centrifugal and Coriolis accelerations
            EarthRate = 7.292115*10^(-5);
            EarthRateVector = [0 0 EarthRate]';
            gravity = xyz2grav(obj.state_(1),obj.state_(2),obj.state_(3));  % ECEF gravity using a model for it
            centrifugalAcc = - cross(  EarthRateVector, cross(EarthRateVector, posOld)  );
            CoriolisAcc = - 2 * cross( EarthRateVector, velOld );
            
            % Correction of the IMU readings
            omegaCorr = omegaMeas - biasGyroOld - C_b2e' * EarthRateVector ;
            acceCorr = C_b2e*(acceMeas-biasAcceOld) + gravity + CoriolisAcc;
            obj.angularRate_ = omegaMeas;
            obj.acceleration_ = acceCorr;
            
            % Non linear equations for the prediction model
            posNew = posOld + dt*velOld + 0.5*dt^2*acceCorr; % Update the position 
            velNew = velOld + dt * acceCorr;          % Update the velocity
            theta = norm(omegaCorr)*dt*0.5; % Quaternion update
            uVector = omegaCorr / norm(omegaCorr);
            deltaQuat = [ cos(theta); uVector * sin(theta) ];
            deltaQuatMatrix = quat_productMatrix(deltaQuat, 'right');
            quatNew = deltaQuatMatrix * quatOld;
            biasAcceNew = biasAcceOld; % Update the biases
            biasGyroNew = biasGyroOld;
            
            % Jacobian of the system dynamics
            Fx = zeros(size(obj.state_));
            
            % Jacobian of the position
            F_pos_pos = eye(3) + 0.5*dt*dt*obj.Jacobian_grav_wrt_pos();
            F_pos_vel  = dt*eye(3);
            F_pos_quat = 0.5*dt*dt * obj.Jacobian_vector_wrt_quaternion(acceMeas);
            F_pos_bAcc = -0.5*dt*dt* C_b2e; 
            F_pos_bGyr = zeros(3);
            F_pos_other = zeros(3, sizeNoINS);
            % Jacobian of the velocity
            F_vel_pos = dt*obj.Jacobian_grav_wrt_pos;
            F_vel_vel  = eye(3) + dt * EarthRate * [ 0,1,0; -1,0,0; 0,0,0 ];
            F_vel_quat = dt * obj.Jacobian_vector_wrt_quaternion(acceMeas);
            F_vel_bAcc = -dt * C_b2e; 
            F_vel_bGyr = zeros(3);
            F_vel_other = zeros(3, sizeNoINS);
            % Jacobian of the quaternion
            F_quat_pos = zeros(4,3);
            F_quat_vel = zeros(4,3);
            F_quat_quat = deltaQuatMatrix;
            F_quat_bAcc = zeros(4,3);
            F_quat_bGyr = obj.Jacobian_quat_wrt_bGyro_vKike(omegaCorr, dt);
            F_quat_other = zeros(4,sizeNoINS);
             % Jacobian of the accelerometer bias
             F_bAcc = [ zeros(3, 10), eye(3), zeros(3), zeros(3, sizeNoINS) ];
             % Jacobian of the gyroscope bias
             F_bGyr = [ zeros(3, 13), eye(3),  zeros(3, sizeNoINS) ];
             % Jacobian of the other elements (phase ambiguities, clockoffsets, others...)
             F_other = [ zeros(sizeNoINS, obj.sizeState_), eye(sizeNoINS) ];
             
             %%%%%%%%%%%%%%%%
%                 F_pos_bAcc = zeros(3);
%                 F_vel_bAcc = zeros(3);
%                 F_quat_bGyr = zeros(4,3);
             %%%%%%%%%%%%%%%%
             
            Fx = [ F_pos_pos, F_pos_vel, F_pos_quat, F_pos_bAcc, F_pos_bGyr, F_pos_other; ...
                      F_vel_pos, F_vel_vel, F_vel_quat, F_vel_bAcc, F_vel_bGyr, F_vel_other; ...
                      F_quat_pos, F_quat_vel, F_quat_quat, F_quat_bAcc, F_quat_bGyr, F_quat_other;
                      F_bAcc;...
                      F_bGyr;...
                      F_other ];
            
            % Control input and its relationship to the system dynamics     
                 u  = [ gravity + centrifugalAcc + CoriolisAcc ];
                 B  = [  dt^2/2 * eye(3); ...
                    dt * eye(3);     ...
                    zeros(4,3);      ...
                    zeros(3);        ...
                    zeros(3);        ...
                    zeros(sizeNoINS,3) ];
            
            % Jacobian of the process noise                
            Fq_pos_acceNoise = 0.5*dt*dt * C_b2e;
            Fq_pos_gyroNoise = zeros(3);
            Fq_pos_biasAcceNoise = 0.5*dt*dt * C_b2e;
            Fq_pos_biasGyroNoise = zeros(3);
            Fq_pos_other = zeros(3,sizeNoINS);
            
            Fq_vel_acceNoise = dt * C_b2e;
            Fq_vel_gyroNoise = zeros(3);
            Fq_vel_biasAcceNoise = dt * dt * C_b2e;
            Fq_vel_biasGyroNoise = zeros(3);
            Fq_vel_other = zeros(3,sizeNoINS);
            
            Fq_quat_acceNoise = zeros(4,3);
            Fq_quat_gyroNoise = obj.Jacobian_quat_wrt_bGyro_vKike(omegaCorr, dt);
            Fq_quat_biasAcceNoise = zeros(4,3);
            Fq_quat_biasGyroNoise = F_quat_bGyr;
            Fq_quat_other = zeros(4,sizeNoINS);
            
            Fq_bAcc_acceNoise = zeros(3);
            Fq_bAcc_gyroNoise = zeros(3);
            Fq_bAcc_biasAcceNoise = eye(3);
            Fq_bAcc_biasGyroNoise = zeros(3);
            Fq_bAcc_other = zeros(3,sizeNoINS);
            
            Fq_bGyr_acceNoise = zeros(3);
            Fq_bGyr_gyroNoise = zeros(3);
            Fq_bGyr_biasAcceNoise = zeros(3);
            Fq_bGyr_biasGyroNoise = eye(3);
            Fq_bGyr_other = zeros(3,sizeNoINS);
            
            Fq_other_inertialNoise = zeros(sizeNoINS, 12);
            Fq_other_other = eye(sizeNoINS);
            
            %%%%%%%%
%             Fq_pos_biasAcceNoise = zeros(3);
%             Fq_vel_biasAcceNoise  = zeros(3);
%             Fq_quat_biasGyroNoise = zeros(4,3);
            %%%%%%%%
            
             Fq = [ Fq_pos_acceNoise,       Fq_pos_gyroNoise,       Fq_pos_biasAcceNoise,       Fq_pos_biasGyroNoise,       Fq_pos_other; ...
                        Fq_vel_acceNoise,           Fq_vel_gyroNoise,       Fq_vel_biasAcceNoise,       Fq_vel_biasGyroNoise,       Fq_vel_other; ...
                        Fq_quat_acceNoise,      Fq_quat_gyroNoise,      Fq_quat_biasAcceNoise,      Fq_quat_biasGyroNoise,      Fq_quat_other; ...
                        Fq_bAcc_acceNoise,      Fq_bAcc_gyroNoise,      Fq_bAcc_biasAcceNoise,      Fq_bAcc_biasGyroNoise,      Fq_bAcc_other;...
                        Fq_bGyr_acceNoise,      Fq_bGyr_gyroNoise,      Fq_bGyr_biasAcceNoise,      Fq_bGyr_biasGyroNoise,      Fq_bGyr_other;...
                        Fq_other_inertialNoise, Fq_other_other                        ];
             
              newStateMatrix = Fx*obj.state_ + B*u;   % One can compare that the update using the non-linear equations and the update using the matrix form are equivalent
                 
             % Mean of the state update
             obj.state_(1:obj.sizeState_) = [ posNew; velNew; quatNew; biasAcceNew; biasGyroNew ];
             obj.state_(7:10) = quat_normalize(obj.state_(7:10)); % Normalization of the quaternion, just to assure that it represents a proper rotation. 
             
             % Covariance Update
             obj.setSystemNoise( diag(   [ diag(obj.Qacc_); diag(obj.Qgyr_); dt*diag(obj.QbiasAcc_); dt*diag(obj.QbiasGyr_); dt*obj.Qamb_*ones(length(obj.satellitesLOS_),1) ] ) ); % Set the noise for the prediction
             newPmatrix = Fx * obj.P_ * Fx' + Fq * obj.Q_ * Fq';
             obj.P_     = newPmatrix; 
             
             % Estimation of the Euler Angles  
             C_b2e = quat_2_DCM(quatNew);
             C_e2n = ECEF_2_NED_rotation(posNew);
             C_b2n = C_e2n * C_b2e;
             obj.eulerAngles_ = DCM_2_euler(C_b2n)*180/pi;
             
            if nargout == 1
                varargout{1} = obj.state_;
            else
                varargout = [];
            end
        end
        
        
        
        %% Correction Step for non-inertial RTK
        function varargout = correctionNonINSRTK(obj,varargin)       
            %%%%%%%%%%%%%
            % TO DO:
            % 1. DONE: Dual frequency?
            % 2. DONE: We are losing the cross-covariance values of the P matrix.
            % 3. GLONASS?
            % 4. Having DD Ambiguities in the state instead of SD
            %%%%%%%%%%%%%
            
            satLOS      = varargin{1};
            satRef      = varargin{2};
            satPos      = varargin{3};
            satRefPos   = varargin{4};
            DDPhase     = varargin{5};
            DDRange     = varargin{6};
            waveLengVec = varargin{7};
            typeObs     = varargin{8};
            iGNSS       = varargin{9};
            nObs        = length(satLOS);
            
%             switch obj.model_
%                 case 'Continuous'
%                     resetAmbiguities = varargin{9};
%                     obj.cycleSlipReset(resetAmbiguities);
%                 case 'Instantaneous'
%                     obj.cycleSlipReset();
%             end
            
            
            if nObs < 2
                if obj.debugMode_ 
                    disp('Not enough observations')
                end
                if nargout == 1
                    varargout{1} = obj.state_;
                else
                    varargout = [];
                end
                return;
            end
            
            if isempty(obj.satRef_)                                         % Check whether we have a reference satellite for the Double Difference observation approach
                obj.satRef_         = satRef;
                obj.satellitesLOS_  = [];
            else
                obj.satRef_         = satRef;
            end
            [~,intNew,intOld]   = intersect(satLOS,obj.satellitesLOS_);     % The intersect of the old and new LOS satellites. tmp1 refers to the indexes of current observed satellites || tmp2 referes to the indexes of the formerly observed ones
            [~,difNew,difOld]   = setxor(satLOS,obj.satellitesLOS_);        % The intersect of the old and new LOS satellites. tmp1 refers to the indexes of current observed satellites || tmp2 referes to the indexes of the formerly observed ones
            obj.satellitesLOS_  = satLOS;                                   % Save in the class which is the PRN of the LOS satellites
            refPosIndx          = [1:1:length(satLOS)];                     % Vble to refer to the position in the vector of satellite in which the reference satellite is
            refPosIndx          = refPosIndx(ismember(satLOS,satRef));
            
            % Adapt the state to the currently observed satellites (adding/removing ambiguities from new/old satellites)
            basePosition = [4041839.1018   537121.6018  4888452.5105];
            oldState            = obj.state_;                               % Copy of the previous state
            initialAmbiguities  = (DDPhase - DDRange)./waveLengVec;         % Approximation for the initial ambiguities (for satellites that we still do not have in our state)
            newState            = [oldState(1:obj.sizeState_); initialAmbiguities.*ones(size(satLOS))];   % Create a new state whose length is equal to the number of observed satellites + 6 (position and velocity)
            newState(obj.sizeState_+intNew)  = oldState(obj.sizeState_+intOld);  % Save the ambiguities of the satellites formerly observed and present in the previous state
            obj.state_          = newState;                                 % Use the current state as the state of the object

            % Adapt the covariance to the currently observed satellites 
            oldP                = obj.P_;                                   % Copy of the previous covariance
            oldP(:,obj.sizeState_+difOld) = [];  oldP(obj.sizeState_+difOld,:) = [];                  % Eliminate from the covariance related to the no-longer observed satellites
            newP                = augmentCovarianceMatrix(oldP, obj.sizeState_+difNew); 
            for iC=1:length(difNew)                                        % Initialize the covariance for the phase ambiguities
                newP(obj.sizeState_+difNew(iC),obj.sizeState_+difNew(iC)) = obj.initialAmbCovariance_; 
            end 
            obj.P_              = newP;                                     % create the new covariance matrix, saving the relevant information from the former covariance matrix
            
            % Jacobian matrix of the solution
            [H, obj.D_]              = obj.H_( obj.state_, obj.sizeState_, satPos, satRefPos, refPosIndx, waveLengVec, typeObs, obj.leverArmIMU2GNSS_ );
            z                   = [DDPhase; DDRange];                       % Pile up phase and range measurements           
            
            % Application of the correction model using the predicted state
            z([refPosIndx,nObs+refPosIndx]) = [];                           % Eliminate the measurements corresponding to the reference satellites
            satPos([refPosIndx],:) = [];                                   
            satRefPos([refPosIndx],:) = [];                                
            waveLengVec2 = waveLengVec; waveLengVec2(refPosIndx) = [];
            h_z                 = obj.h_( obj.state_, satPos, satRefPos, obj.state_(obj.sizeState_+1:end), waveLengVec2, obj.D_, obj.leverArmIMU2GNSS_); % Observation model
            obj.y_              = z - h_z; % H*obj.state_;%h_z%;                                  % Innovation: difference between the observed measurements and the observation model which relates observations and state         
            
            % Building R matrix for the observations
            R_values            = diag(obj.R_);                             % the R values include the covariances for the phase and code measurements 
            R_phase             = R_values(1:nObs);
            R_code              = R_values(nObs+1:end);
            obj.R_              = [obj.D_*diag(R_phase)*obj.D_',   zeros(size(obj.D_*diag(R_code)*obj.D_'));        zeros(size(obj.D_*diag(R_phase)*obj.D_')),    obj.D_*diag(R_code)*obj.D_'];
             
            % Correction Step
%             obj.S_              = H*obj.P_*H' + obj.R_;                     % Innovation covariance
%             obj.K_              = obj.P_*H'/obj.S_;                         % Kalman Gain
%             innovation          = obj.K_ * obj.y_;
%             obj.state_          = obj.state_ + innovation;             % update the mean of the state
%             obj.P_              = (eye(size(obj.P_,1)) - obj.K_ * H) * obj.P_;  % update the state covariance.
            
            [obj]=KF_OD_loop(obj,z,h_z,H,basePosition,iGNSS);
            
            % Saving the single and double difference phase ambiguities
            obj.SDAmb_          = obj.state_(obj.sizeState_+1:end);         % Saving the single and double difference ambiguities of the satellites
            obj.DDAmb_          = obj.D_ * obj.SDAmb_;
            
            % Estimate the Q_RN, needed to find the fixed solution afterwards
            G_aux               = zeros(length(obj.state_)-max(typeObs),length(obj.state_));
            G_aux(1:obj.sizeState_,1:obj.sizeState_) = eye(obj.sizeState_);
            G_aux(obj.sizeState_+1:end,obj.sizeState_+1:end) = obj.D_;
            P_dd                = G_aux*obj.P_*G_aux';
            obj.QNR_            = P_dd(obj.sizeState_+1:end,1:obj.sizeState_);
            obj.QRN_            = P_dd(1:obj.sizeState_,obj.sizeState_+1:end);
            obj.QN_             = P_dd(obj.sizeState_+1:end,obj.sizeState_+1:end);
            
            %%%%%%%%%%% Trick used by Anja -> Change of the Covariance matrix
            if obj.resetDDCrossVarianceP_
                Pauxx = obj.P_;
                Pauxx(1:obj.sizeState_,:) = 0;
                Pauxx(:,1:obj.sizeState_) = 0;
                Pauxx(1:obj.sizeState_,1:obj.sizeState_) = diag(diag(obj.P_(1:obj.sizeState_,1:obj.sizeState_)));
                obj.P_ = Pauxx;
            end
            %%%%%%%%%%%
            if nargout == 1
                varargout{1} = obj.state_;
            else
                varargout = [];
            end           
        end
        
        %% Correction Step (while using INS)
        function varargout = correctionInertial(obj,varargin)       
            %%%%%%%%%%%%%
            % TO DO:
            % 1. DONE: Dual frequency?
            % 2. DONE: We are losing the cross-covariance values of the P matrix.
            % 3. GLONASS?
            % 4. Having DD Ambiguities in the state instead of SD
            %%%%%%%%%%%%%
            
            satLOS      = varargin{1};
            satRef      = varargin{2};
            satPos      = varargin{3};
            satRefPos   = varargin{4};
            DDPhase     = varargin{5};
            DDRange     = varargin{6};
            waveLengVec = varargin{7};
            typeObs     = varargin{8};
            nObs        = length(satLOS);
            
            if nObs < 2
                if obj.debugMode_ 
                    disp('Not enough observations')
                end
                if nargout == 1
                    varargout{1} = obj.state_;
                else
                    varargout = [];
                end
                return;
            end
            
            if isempty(obj.satRef_)                                         % Check whether we have a reference satellite for the Double Difference observation approach
                obj.satRef_         = satRef;
                obj.satellitesLOS_  = [];
            else
                obj.satRef_         = satRef;
            end
            [~,intNew,intOld]   = intersect(satLOS,obj.satellitesLOS_);     % The intersect of the old and new LOS satellites. tmp1 refers to the indexes of current observed satellites || tmp2 referes to the indexes of the formerly observed ones
            [~,difNew,difOld]   = setxor(satLOS,obj.satellitesLOS_);        % The intersect of the old and new LOS satellites. tmp1 refers to the indexes of current observed satellites || tmp2 referes to the indexes of the formerly observed ones
            obj.satellitesLOS_  = satLOS;                                   % Save in the class which is the PRN of the LOS satellites
            refPosIndx          = [1:1:length(satLOS)];                     % Vble to refer to the position in the vector of satellite in which the reference satellite is
            refPosIndx          = refPosIndx(ismember(satLOS,satRef));
            
            % Adapt the state to the currently observed satellites (adding/removing ambiguities from new/old satellites)
            oldState            = obj.state_;                               % Copy of the previous state
            initialAmbiguities  = (DDPhase - DDRange)./waveLengVec;         % Approximation for the initial ambiguities (for satellites that we still do not have in our state)
            newState            = [oldState(1:obj.sizeState_); initialAmbiguities.*ones(size(satLOS))];   % Create a new state whose length is equal to the number of observed satellites + 6 (position and velocity)
            newState(obj.sizeState_+intNew)  = oldState(obj.sizeState_+intOld);  % Save the ambiguities of the satellites formerly observed and present in the previous state
            obj.state_          = newState;                                 % Use the current state as the state of the object

            % Adapt the covariance to the currently observed satellites 
            oldP                = obj.P_;                                   % Copy of the previous covariance
            oldP(:,obj.sizeState_+difOld) = [];  oldP(obj.sizeState_+difOld,:) = [];                  % Eliminate from the covariance related to the no-longer observed satellites
            newP                = augmentCovarianceMatrix(oldP, obj.sizeState_+difNew); 
            for iC=1:length(difNew)                                        % Initialize the covariance for the phase ambiguities
                newP(obj.sizeState_+difNew(iC),obj.sizeState_+difNew(iC)) = obj.initialAmbCovariance_; 
            end 
            obj.P_              = newP;                                     % create the new covariance matrix, saving the relevant information from the former covariance matrix
            
            % Jacobian matrix of the solution
            [H, obj.D_]              = obj.H_( obj.state_, obj.sizeState_, satPos, satRefPos, refPosIndx, waveLengVec, typeObs, obj.leverArmIMU2GNSS_ );
            z                   = [DDPhase; DDRange];                       % Pile up phase and range measurements           
            
            % Application of the correction model using the predicted state
            z([refPosIndx,nObs+refPosIndx]) = [];                           % Eliminate the measurements corresponding to the reference satellites
            satPos([refPosIndx],:) = [];                                   
            satRefPos([refPosIndx],:) = [];                                
            waveLengVec2 = waveLengVec; waveLengVec2(refPosIndx) = [];
            h_z                 = obj.h_( obj.state_, satPos, satRefPos, obj.state_(obj.sizeState_+1:end), waveLengVec2, obj.D_, obj.leverArmIMU2GNSS_); % Observation model
            obj.y_              = z - h_z;                                  % Innovation: difference between the observed measurements and the observation model which relates observations and state         
            
            % Building R matrix for the observations
            R_values            = diag(obj.R_);                             % the R values include the covariances for the phase and code measurements 
            R_phase             = R_values(1:nObs);
            R_code              = R_values(nObs+1:end);
            obj.R_              = [obj.D_*diag(R_phase)*obj.D_',   zeros(size(obj.D_*diag(R_code)*obj.D_'));        zeros(size(obj.D_*diag(R_phase)*obj.D_')),    obj.D_*diag(R_code)*obj.D_'];
             
            % Correction Step
            obj.S_              = H*obj.P_*H' + obj.R_;                     % Innovation covariance
            obj.K_              = obj.P_*H'/obj.S_;                         % Kalman Gain
            innovation          = obj.K_ * obj.y_;
            obj.state_          = obj.state_ + innovation;             % update the mean of the state
            obj.P_              = (eye(size(obj.P_,1)) - obj.K_ * H) * obj.P_;  % update the state covariance.
            
            obj.state_(7:10)    = quat_normalize(obj.state_(7:10));
            % Estimation of the Euler Angles  
             C_b2e = quat_2_DCM(obj.state_(7:10));
             C_e2n = ECEF_2_NED_rotation(obj.state_(1:3));
             C_b2n = C_e2n * C_b2e;
             obj.eulerAngles_ = DCM_2_euler(C_b2n)*180/pi;
            

            % Saving the single and double difference phase ambiguities 
            obj.SDAmb_          = obj.state_(obj.sizeState_+1:end);         % Saving the single and double difference ambiguities of the satellites
            obj.DDAmb_          = obj.D_ * obj.SDAmb_;
            
            % Estimate the Q_RN, needed to find the fixed solution afterwards
            G_aux               = zeros(length(obj.state_)-max(typeObs),length(obj.state_));
            G_aux(1:obj.sizeState_,1:obj.sizeState_) = eye(obj.sizeState_);
            G_aux(obj.sizeState_+1:end,obj.sizeState_+1:end) = obj.D_;
            P_dd                = G_aux*obj.P_*G_aux';
            obj.QNR_            = P_dd(obj.sizeState_+1:end,1:obj.sizeState_);
            obj.QRN_            = P_dd(1:obj.sizeState_,obj.sizeState_+1:end);
            obj.QN_             = P_dd(obj.sizeState_+1:end,obj.sizeState_+1:end);
            
            
            
            
            
            
            if nargout == 1
                varargout{1} = obj.state_;
            else
                varargout = [];
            end           
        end
        
        %% Fix and Hold Mode ->  Documented in the Weekly Presentation and based on RTK Manual, pages 166-167
        function varargout = FixNHoldMode(obj, varargin)
            FixDD = varargin{1};
            nObs = length(FixDD);
            z_k = FixDD;
            H_k = [zeros(nObs,obj.sizeState_), obj.D_ ];
            
            sigmaFixHold = 0.001^2;
            R = sigmaFixHold * eye(nObs);
            
            S = H_k *obj.P_*H_k' + R;                     % Innovation covariance
            K = obj.P_*H_k'/S;                                  % Kalman Gain
            
            obj.state_          = obj.state_ + K * (z_k - H_k*obj.state_);             % update the mean of the state
            obj.P_              = obj.P_ - K * H_k * obj.P_;             % update the state covariance.
    
            if nargout == 1
                varargout{1} = obj.state_;
            else
                varargout = [];
            end           
        end
            
        %%  Loosely-Coupled Velocity Update using the Doppler measurements
        function varargout = correctionLCVelocity(obj, varargin)
            if nargin == 3
                measuredVel = varargin{1};               % we are expecting to receive directly the estimated velocity solution from maybe a SPP solver
                measurementR = varargin{2};             % we are supposed to input already a R_[3x3] matrix. You could use just a regular I_[3x3] * \Sigma_{PRR}, but also some model based on the elevation of the satellites
            end
        
            H = [zeros(3);  eye(3);  zeros( length(obj.state_)-6,3 )]';             % Jacobian matrix for the direction observation of the velocity via Doppler measurements
            KalmanGain = obj.P_ * H' / (H*obj.P_*H' +measurementR) ;
            innovation = KalmanGain * ( measuredVel - obj.state_(4:6) );
            obj.state_ = obj.state_ + innovation;
            obj.P_              = (eye(size(obj.P_,1)) - KalmanGain * H) * obj.P_;  % update the state covariance.
            
            if nargout == 1
                varargout{1} = obj.state_;
            else
                varargout = [];
            end
            
        end
        
        %% Reset of the phase ambiguities
        function cycleSlipReset( obj, resetLOS )
           if nargin == 1 % Complet reset of the ambiguities
               obj.state_  = obj.state_(1:obj.sizeState_);
               obj.P_      = obj.P_(1:obj.sizeState_, 1:obj.sizeState_);
               obj.Q_      = obj.Q_(1:obj.sizeSystemNoise_,1:obj.sizeSystemNoise_);
               obj.satellitesLOS_ = [];
               obj.satRef_ = [];
               obj.QN_ = [];
               obj.DDAmb_ = [];
           elseif nargin == 2 % Only reset those ambiguities which are detected to have a cycle slip
               [~,sat2reset] = intersect(obj.satellitesLOS_, resetLOS);
               obj.satellitesLOS_(sat2reset) = [];
               obj.state_(obj.sizeState_ + sat2reset) = [];
               obj.P_( obj.sizeState_ + sat2reset, : ) = [];  obj.P_( :, obj.sizeState_ + sat2reset ) = [];
               obj.Q_ = blkdiag(obj.Q_(1:obj.sizeSystemNoise_,1:obj.sizeSystemNoise_), obj.Qamb_ * eye(length(obj.satellitesLOS_)) );
               obj.QN_ = [];
               obj.DDAmb_ = [];
           else
              disp('EKF RTK CycleSlipReset: wrong number of inputs');
           end              
        end
        
        %% Jacobian of gravity w.r.t position
        function J_grav_pos = Jacobian_grav_wrt_pos (obj )
            % References for this function: Groves Book, chapter 12, page
            % 383 for the Jacobian. The models for the gravity in chapter 2.
            e = 0.0818191908425;
            R0 = 6378137;
            position = obj.state_(1:3);
            [ Lat, ~ , ~ ]    = cart2geod(position(1), position(2),position(3) );
            Re = R0/sqrt(1-e^2*sin(Lat)^2);
            r_eS = Re * sqrt( cos(Lat)^2 + (1-e^2)*sin(Lat)^2 );
            g0 = 9.7803253359*(1+0.001931853*sin(Lat)^2)/(sqrt(1-e^2*sin(Lat)^2));
            J_grav_pos = 2*g0/r_eS*position/norm(position)*position';
            J_grav_pos = zeros(3);
        end
        
        %% Jacobian of the rotation of a vector w.r.t. quaternion
        function J_vect_quat = Jacobian_vector_wrt_quaternion ( obj, vector )
            vector = reshape(vector,3,1);
            quaternion = obj.state_(7:10);
            qw = quaternion(1);
            qu = quaternion(2:4);
            J_vect_quat = 2 * [ qw*vector + cross(qu, vector), qu'*vector*eye(3) + qu*vector' - vector*qu' - qw*skewMatrix(vector)  ];
        end
        
        %% Jacobian of the quaternion w.r.t. bias gyroscope
        function J_quat_bGyr = Jacobian_quat_wrt_bGyro ( obj, angularRate, dt )
            quaternion = obj.state_(7:10);
            qw = quaternion(1);
            qu = quaternion(2:4);
            biasGyro = obj.state_(14:16);
            omega = angularRate;
            omegaNorm = norm(omega);
            theta = omegaNorm * dt / 2;
            J_quat_bGyr = dt/2*sin(theta);
        end
        
        function J_quat_bGyr = Jacobian_quat_wrt_bGyro_vKike( obj, angularRate, dt )
            quaternion = obj.state_(7:10);
            biasGyro = obj.state_(14:16);
            theta = norm(0.5*dt*angularRate);
            Omega = zeros(4);
            Omega = [ 0, angularRate'; -angularRate, skewMatrix(angularRate) ] ;
            omega_x = [ 0 1 0 0; -1 0 0 0; 0 0 0 -1; 0 0 1 0  ];
            omega_y = [ 0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0 ];
            omega_z = [ 0 0 0 1; 0 0 -1 0; 0 1 0 0; -1 0 0 0 ];            
            term1_x = eye(4) * sin( theta ) * biasGyro(1)/theta;
            term2_x = 1/2 * dt * omega_x * sin (theta)/theta;
            term3_x = 1/2*dt^2 * Omega * ( theta*cos(theta)-sin(theta) )/theta^3;
            term_x = -(term1_x + term2_x + term3_x);
            term1_y = eye(4) * sin( theta ) * biasGyro(2)/theta;
            term2_y = 1/2 * dt * omega_y * sin (theta)/theta;
            term3_y = 1/2*dt^2 * Omega * ( theta*cos(theta)-sin(theta) )/theta^3;
            term_y = -(term1_y + term2_y + term3_y);            
            term1_z = eye(4) * sin( theta ) * biasGyro(3)/theta;
            term2_z = 1/2 * dt * omega_z * sin (theta)/theta;
            term3_z = 1/2*dt^2 * Omega * ( theta*cos(theta)-sin(theta) )/theta^3;
            term_z = -(term1_z + term2_z + term3_z);
            J_quat_bGyr = [ term_x*quaternion, term_y*quaternion, term_z*quaternion ];
        end
        
        %% Set the Q matrix for the system's noise
        function setSystemNoise (obj, noise)
            % Function to set the matrix containing the noise of the system
            % (to be used during the prediction step)
            if size(noise,1) == size(noise,2) % the input noise is given in matrix form
                obj.Q_ = noise;
            else
                obj.Q_ = diag(noise);
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