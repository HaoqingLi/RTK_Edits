classdef ExtendedKalmanFilterIMM < AISIntegrityBasePack.AbstractKalmanFilter
    properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% process noise definition %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        eOmega_ % process noise ROT (sigma)
        eVel_ % process noise SOG (sigma)
        eAcc_ % process noise acceleration (sigma)
        ePhi_ % process noise true heading rate (sigma)
        
        Jf % time variant Jacobian of state transition function (derivative of state)
        Jc % time variant Jacobian of state transition function (derivative of error)
        
        bIMMActivated_ = true; % switch to either use both models or not
        bIMMCVModelUsed_ = false; % flag if CV Model is currently in use
        bCVModelOnly_ = false;
        
        mapPModel2VelInState_ = containers.Map({'CV','CTRV'},{4,4}); % true == CV, false == CTRV (could be adopted)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%% INTEGRITY MEASURE %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        dequeResidualsSquaredNorm % stack of residuals (squared and normalized)
        dequeResidualsNormed % stack of residuals (normalized)
        dequeResiduals % stack of residuals
        cdfResiduals % cumulative sum of residuals
        NWindow = 3 % window length of cumsum of residuals
        NMaxRes = 1001 % max size of residuals (TODO: number of snapshots, or related to time?)
        
        tracePredP_ % queue to trace P_ (after predicition)      
        traceFailure_ % queue to trace detected failures
        
        GLRInst_ % instance of GLR (Generalized Likelihood Ratio) test class
        GLRInstCV_ % instance of GLR (Generalized Likelihood Ratio) test class for CV model
        GLR_DET_THSLD = 12.1076;
        GLR_DET_THSLD_CV = 4;
        GLR_WIN_SIZE = 10;
        GLR_FTYPE = {'jump', 'step', 'ramp'};
        
        dShipL_ = 29; % length of BT2 (default)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  Fall back Process Model (CV) %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fCV_ % process model for CV (runs in parallel)
        FCV_ % Jacobian of f (w.r.t. state)
        CCV_ % Jacobian of f (w.r.t. process noise)
        QCV_ % process noise covariance
        PCV_ % error noise covariance
        JfCV % intermediate Jacobian (J) (no handles)
        JcCV % intermediate Jacobian (C) (no handles)
        stateEstCV_ % state vector for CV process model
                
        hCV_ % measurement model for CV (needs to map to CV state vector)
        HCV_ % Jacobian of h (w.r.t. state)
        DCV_ % Jacobian of h (w.r.t. process noise)
        RCV_ % measurement noise covariance
        SCV_ % covariance residual CV
        KCV_ % Kalman Gain CV
        yCV_ % innovation CV
        
        tracePredPCV_ % trace of error covariance for CV
        
        psiCVbak_ % backup of last psi value for CV
        
        stateCVbak_ % backup last state of CV (before resetting position coordinates)
        stateEstbak_ % backup last state of CTRV (before resetting position coordinates)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% chi square LUT for testing H0 that NWindow*q_m is chi square
        %%% distributed with probability 1 - alpha
        %%% q_m is expectation of normalized squared residuals of length
        %%% NWindow
        %%% degrees of freedom (key of hash maps below) is dim(z)*NWindow!
        %%% so, we look for P(N*q_m in [r1,r2]|H0) = 1 - alpha, 
        %%% where [r1,r2] is 1-alpa confidence interval
        %%% alpha = 0.05 -> we are bounded by ChiSquare0.025 and
        %%% ChiSquare0.975 (ChiSquareALPHA NOT ChiSquareP!!!)
        ChiSquareConfPerc0p975 = containers.Map((1:30),...
            [0.001 0.051 0.216 0.484 0.831 1.237 1.690 2.180 2.700 3.247 3.816 4.404 5.009 5.629 6.262 ...
            6.908 7.564 8.231 8.907 9.591 10.283 10.982 11.689 12.401 13.120 13.844 14.573 15.308 16.047 16.791]);
        ChiSquareConfPerc0p025 = containers.Map((1:30),...
            [5.024 7.378 9.348 11.143 12.833 14.449 16.013 17.535 19.023 20.483 21.920 23.337 24.736 26.119 27.488 ...
            28.845 30.191 31.526 32.852 34.170 35.479 36.781 38.076 39.364 40.646 41.923 43.195 44.461 45.722 46.979]);
        ChiSquareConfPerc0p050 = containers.Map((1:30),...
            [3.841 5.991 7.815 9.488 11.070 12.592 14.067 15.507 16.919 18.307 19.675 21.026 22.362 23.685 24.996 26.296 ...
             27.587 28.869 30.144 31.410 32.671 33.924 35.172 36.415 37.652 38.885 40.113 41.337 42.557 43.773]);
        
        dequeVelTracker_ = [];
        VEL_MIN_THSLD = .5;
        
        bDetThshldGLR_ = false; % flag to be set for estimation of detection threshold of GLR
                                % in that case, reference data is needed to
                                % compute new residuals (from EKF estimate to reference (not measurement!!))
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Interacting Multiple Model (IMM) %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        mapIMMstate2model_ = containers.Map({'CTRV','CV'},{[],[]});
        
        % error state covariance (for common elements) P^{m}_c
        mapIMMPCommon2model_ = containers.Map({'CTRV','CV'},{[],[]});
        % error state covariance (for extra elements) P^{m}_e
        mapIMMPExtra2model_ = containers.Map({'CTRV','CV'},{[],[]});
        % error state covariance (for cov common<->extra elements) P^{m}_ce
        mapIMMPce2model_ = containers.Map({'CTRV','CV'},{[],[]});
        % error state covariance (for cov extra<->common elements) P^{m}_ec
        % (should be identical to [P^{m}_ce]^T)
        mapIMMPec2model_ = containers.Map({'CTRV','CV'},{[],[]});
        
        % P^{m|n}_c
        mapIMMPcM2M_ = containers.Map({'CVCV','CVCTRV','CTRVCV','CTRVCTRV'},{[],[],[],[]});
        % P^{m|n}_e
        mapIMMPeM2M_ = containers.Map({'CVCV','CVCTRV','CTRVCV','CTRVCTRV'},{[],[],[],[]});
        % P^{m|n}_ec
        mapIMMPecM2M_ = containers.Map({'CVCV','CVCTRV','CTRVCV','CTRVCTRV'},{[],[],[],[]});
        % P^{m|n}_ce
        mapIMMPceM2M_ = containers.Map({'CVCV','CVCTRV','CTRVCV','CTRVCTRV'},{[],[],[],[]});
        
        mapIMMmodel2idx_ = containers.Map({'CTRV','CV'},{1,2});                        
        mapIMMmodelProb_ = containers.Map({'CTRV','CV'},{[],[]});
        mapIMMaposterioriState_ = containers.Map({'CTRV','CV'},{[],[]});
        mapIMMaposterioriErrCovP_ = containers.Map({'CTRV','CV'},{[],[]});
        mapIMMnormC_ = containers.Map({'CTRV','CV'},{[],[]});
        
        dIMMllhood_ = containers.Map({'CTRV','CV'},{[],[]});
        dIMMtransProb_ = [0.8 0.2; 0.2 0.8];
        dIMMmixingWeights_ = 0.5*ones(2);
        dIMMerrCovP_
        dIMMstateCombined_
                                
    end
    
    methods 
        function obj = ExtendedKalmanFilterIMM(varargin)
            
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
            
            obj@AISIntegrityBasePack.AbstractKalmanFilter(superclassargs{:});
            
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
                    elseif strcmp(varargin{1},'-dynModelFB') && ...
                            isa(varargin{2},'function_handle')
                        obj.fCV_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-measModelFB') && ...
                            isa(varargin{2},'function_handle')
                        obj.hCV_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-labelsTitles') && ...
                            iscell(varargin{2})
                        obj.cTitleLabels_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-labelsYAxes') && ...
                            iscell(varargin{2})
                        obj.cYlabelLabels_ = varargin{2};
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
                    elseif strcmp(varargin{1},'-D') && ...
                            ismatrix(varargin{2})
                        obj.D_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-FJfallback') && ...
                            ismatrix(varargin{2})
                        obj.FCV_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-Qfallback') && ...
                            ismatrix(varargin{2})
                        obj.QCV_ = varargin{2};
                        varargin{2} = [];
                   elseif strcmp(varargin{1},'-CJfallback') && ...
                            ismatrix(varargin{2})
                        obj.CCV_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-Q') && ...
                            isa(varargin{2},'function_handle')
                        obj.Q_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-HJfallback') && ...
                            ismatrix(varargin{2})
                        obj.HCV_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-DJfallback') && ...
                            ismatrix(varargin{2})
                        obj.DCV_ = varargin{2};
                        varargin{2} = [];   
                    elseif strcmp(varargin{1},'-Rfallback') && ...
                            ismatrix(varargin{2})
                        obj.RCV_ = varargin{2};
                        varargin{2} = [];                       
                    elseif strcmp(varargin{1},'-state') && ...
                            isvector(varargin{2})
                        obj.stateEst_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-stateFB') && ...
                            isvector(varargin{2})
                        obj.stateEstCV_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-bIMMActivated') && ...
                            islogical(varargin{2})
                        obj.bIMMActivated_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-bCVModelOnly') && ...
                            islogical(varargin{2})
                        obj.bCVModelOnly_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-iVelThsld') && ...
                            isnumeric(varargin{2})
                        obj.VEL_MIN_THSLD = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-ChiSWinSize') && ...
                            isnumeric(varargin{2})
                        obj.NWindow = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-GLRDetThsld') && ...
                            isa(varargin{2},'containers.Map')
                        obj.GLR_DET_THSLD = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-GLRDetThsldCV') && ...
                            isa(varargin{2},'containers.Map')
                        obj.GLR_DET_THSLD_CV = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-GLRWinSize') && ...
                            isnumeric(varargin{2})
                        obj.GLR_WIN_SIZE = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-GLRFaultType') && ...
                            iscell(varargin{2})
                        obj.GLR_FTYPE = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-GLREstimateThsldMode') && ...
                            islogical(varargin{2})
                        obj.bDetThshldGLR_ = varargin{2};
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
            
            obj.GLRInst_ = AISIntegrityBasePack.GLRClass('-thsld',obj.GLR_DET_THSLD,...
                                    '-windowlength',obj.GLR_WIN_SIZE,...
                                    '-dim',length(obj.stateEst_));
            if obj.bIMMActivated_                    
                obj.GLRInstCV_ = AISIntegrityBasePack.GLRClass('-thsld',obj.GLR_DET_THSLD_CV,...
                    '-windowlength',obj.GLR_WIN_SIZE,...
                    '-dim',length(obj.stateEstCV_));
            end
        end
        
        function init(obj)
            %%% init process model %%%
            %  define non-linear function depending on dt, state (at k - 1)
            %  and input controls (at k - 1) to propagate state to future
            if isempty(obj.f_)
                obj.f_ = @(dt,px,py,psi,omega,vel,phi) [px + dt*vel * cosd(psi+dt*omega); ...
                            py + dt * vel * sind(psi + dt*omega); ...
                            psi + dt*omega; ...
                            vel; ...
                            phi]; 
            end
            % define Jacobian matrix for f (derivative of state)
            if isempty(obj.F_)
                obj.F_ = @(dt,vel,psi,omega)[ 1 0 -dt*vel*sind(psi) dt*cosd(psi) 0; ...
                            0 1 dt*vel*cosd(psi) dt*sind(psi) 0; ...
                            0 0 1 0 0; ...
                            0 0 0 1 0; ...
                            0 0 0 0 1];
            end
            % create Jacobian of f (derivative of process noise)
            if isempty(obj.C_)
                obj.C_ = @(dt,vel,psi,omega)[dt*cosd(psi) 0 0 0; ...
                            dt*sind(psi) 0 0 0;...
                            0 dt 0 0;...
                            0 0 dt 0;...
                            0 0 0 dt];
            end
                 
            if isempty(obj.stateEst_)
                obj.stateEst_ = [0; 0; 0; 0; 0]; % zero state init
            end
            
            if isempty(obj.stateEstCV_)
                obj.stateEstCV_ = [0; 0; 0; 0]; % zero state init
            end
            
            if obj.bIMMActivated_
                obj.psiCVbak_ = obj.stateEstCV_(3);
            end
                        
            %%% init measurement model %%% 
            % create measurement matrix
            if isempty(obj.H_)
                obj.H_ = eye(length(obj.stateEst_));
            end
            
            obj.eVel_ = 1;
            obj.eOmega_ = .5/180*pi;
            %%% initializing noise covariances based on educated guess
            sigmaPosX = 3; % sigma for position estimate in x [m]
            sigmaPosY = 3; % sigma for position estimate in y [m]

            if isempty(obj.R_)
                obj.R_ = diag([sigmaPosX^2; sigmaPosY^2]); % construct to measurement covariance matrix
            end
            % sigmaProc = 0.1; % guess for sigma process noise (deg/s^2)
            if isempty(obj.Q_)
                obj.Q_ = @(~) diag([obj.eVel_^2; obj.eOmega_^2]); % init process covariance            
            end

            % init measurement noise covariance matrix
            
            obj.iR_ = .5;
            
            if isempty(obj.P_)
                obj.P_ = obj.iR_ * eye(length(obj.stateEst_));
            end
            
            if obj.bIMMActivated_ && isempty(obj.PCV_)
                obj.PCV_ = obj.iR_ * eye(length(obj.stateEstCV_));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%
            %%% init residuals %%%
            %%%%%%%%%%%%%%%%%%%%%%
            obj.dequeResidualsSquaredNorm = []; %-Inf * ones(obj.NMaxRes,1);
            obj.dequeResidualsNormed = []; %-Inf * ones(size(obj.H_,1),obj.NMaxRes);
            obj.dequeResiduals = [];% -Inf * ones(size(obj.H_,1),obj.NMaxRes);
            
            for model = obj.dIMMllhood_.keys()
                obj.dIMMllhood_(model{1}) = 1 / numel(obj.dIMMllhood_.keys());
                obj.mapIMMmodelProb_(model{1}) = 1 / numel(obj.dIMMllhood_.keys());
            end
            
            if obj.bIMMActivated_
                obj.computeIMMmixingStates(1);
                
                psiCV = obj.stateEstCV_(3);
                velCV = obj.stateEstCV_(4);
                % prediction of state covariance
                obj.JfCV = obj.FCV_(0,velCV,psiCV);
                obj.JcCV = obj.CCV_(0,velCV,psiCV);
            end

            if obj.bIMMActivated_ || ~obj.bCVModelOnly_ 
                psi = obj.stateEst_(3);
                vel = obj.stateEst_(4);
                rot = obj.stateEst_(5);
                % obj.checkConstraints();  
                % prediction of state covariance
                obj.Jf = obj.F_(0,vel,psi,rot);
                obj.Jc = obj.C_(0,vel,psi,rot);
            elseif obj.bCVModelOnly_ && ~obj.bIMMActivated_
                psi = obj.stateEst_(3);
                vel = obj.stateEst_(4);
                % obj.checkConstraints();  
                % prediction of state covariance
                obj.Jf = obj.F_(0,vel,psi);
                obj.Jc = obj.C_(0,vel,psi);                      
            end
            
        end
        
        function varargout = prediction(obj,dt)
            %%% controlIn (i.e. u == control input) is expected to be a
            %%% struct containing rot and vel values (from AIS message)
            
            % prediction of state
            % assignment of different variables strongly depends on state
            % setup! for that reason case study        
            if obj.bIMMActivated_
                %%% CV model is assumed
                psiCV = obj.stateEstCV_(3);
                % rotCV = obj.stateEstCV_(4);
                velCV = obj.stateEstCV_(4);
                %%% CTRV is assumed
                obj.stateEstCV_ = obj.fCV_(dt,obj.stateEstCV_(1),obj.stateEstCV_(2),psiCV,velCV);
                % obj.checkConstraints();  
                % prediction of state covariance
                obj.JfCV = obj.FCV_(dt,velCV,psiCV);
                obj.JcCV = obj.CCV_(dt,velCV,psiCV);
                
                obj.PCV_ = obj.JfCV * obj.PCV_ * obj.JfCV' + obj.JcCV * obj.QCV_() * obj.JcCV';
            end    
                
            if obj.bIMMActivated_ || ~obj.bCVModelOnly_ 
                %%% CTRV model is assumed
                psi = obj.stateEst_(3);
                vel = obj.stateEst_(4);
                rot = obj.stateEst_(5);
                obj.stateEst_ = obj.f_(dt,obj.stateEst_(1),obj.stateEst_(2),psi,rot,vel);
                obj.Jf = obj.F_(dt,vel,psi,rot);
                obj.Jc = obj.C_(dt,vel,psi,rot); 
            elseif obj.bCVModelOnly_ && ~obj.bIMMActivated_
                psi = obj.stateEst_(3);
                vel = obj.stateEst_(4);
                % obj.checkConstraints();  
                % prediction of state covariance
                obj.stateEst_ = obj.f_(dt,obj.stateEstCV_(1),obj.stateEstCV_(2),psi,vel);
                obj.Jf = obj.F_(dt,vel,psi);
                obj.Jc = obj.C_(dt,vel,psi);        
            end
            % obj.checkConstraints();  
            % prediction of state covariance                               
            
            obj.P_ = obj.Jf * obj.P_ * obj.Jf' + obj.Jc * obj.Q_() * obj.Jc';
            % obj.P_ = Jf * obj.P_ * Jf' + obj.Q_();
            
            obj.tracePredP_ = [obj.tracePredP_ obj.P_];
            
            if obj.bIMMActivated_
                obj.tracePredPCV_ = [obj.tracePredPCV_ obj.PCV_];
            end
            
            if obj.bIMMActivated_ && obj.bIMMCVModelUsed_
                tmpOmega = (obj.stateEstCV_(3) - psiCV) / dt;
%                 if obj.stateEstCV_(3) - psiCV > 0
%                     tmpOmega = -tmpOmega;
%                 end
                % stateTmpOut = [obj.stateEstCV_(1:3); tmpOmega; obj.stateEstCV_(4)];
                stateTmpOut = [obj.stateEstCV_; tmpOmega]; %obj.stateEst_(5)];
            else
                stateTmpOut = obj.stateEst_;
            end
            
            if nargout == 1
                varargout{1} = stateTmpOut;
            else
                varargout = [];
            end
        end
        
        function varargout = update(obj,sObs,dt,dRef)       
            persistent bSecondModel;
      
            if isempty(bSecondModel)
                bSecondModel = false;             
            end
            
            %%% compute mean velocity
            % mVel = mean(obj.dequeVelTracker_);
%             if obj.bIMMActivated_ && mVel < obj.VEL_MIN_THSLD
%                 if ~bSecondModel
%                     fprintf('switched to slow motion model ...\n')
%                     bSecondModel = true;
%                 end
%                 obj.bIMMCVModelUsed_ = true;
%             else
%                 if bSecondModel
%                     fprintf('switched to fast motion model ...\n')
%                     bSecondModel = false;
%                 end
%                 obj.bIMMCVModelUsed_ = false;
%             end                
                                             
            if obj.bIMMActivated_
                %%% fall back to second motion model (e.g. CV)
                % residual covariance
                obj.SCV_ = obj.HCV_*obj.PCV_*obj.HCV_' + ...
                    obj.DCV_ * obj.RCV_ * obj.DCV_';
            
                % Kalman Gain
                obj.KCV_ = obj.PCV_*obj.HCV_'/obj.SCV_;
                            
                % compute innovation
                obj.yCV_ = (sObs - obj.hCV_(obj.stateEstCV_));
                % obj.y_(3) = mod(abs(obj.y_(3)),2*pi);
                % update error covariance estimation.
                obj.PCV_ = (eye(size(obj.KCV_,1))-obj.KCV_*obj.HCV_)*obj.PCV_;
                % update state estimate
                obj.stateEstCV_ = obj.stateEstCV_ + obj.KCV_ * obj.yCV_;
            end
                
            % residual covariance
            obj.S_ = obj.H_*obj.P_*obj.H_' + obj.D_ * obj.R_ * obj.D_';

            % Kalman Gain
            obj.K_ = obj.P_*obj.H_'/obj.S_;

            % compute innovation
            obj.y_ = (sObs(1:size(obj.h_(obj.stateEst_),1)) - obj.h_(obj.stateEst_));
%                 if length(obj.y_) == 3
%                     obj.y_(3) = mod(abs(obj.y_(3)),2*pi);
%                 end
            % update error covariance estimation.
            obj.P_ = (eye(size(obj.K_,1))-obj.K_*obj.H_)*obj.P_; % that's the
            % same?!
                        
           
            % update the state estimate.
            obj.stateEst_ = obj.stateEst_ + obj.K_ * obj.y_;
            
            %%% perform constraint check
            obj.checkConstraints();
            
            %%% backup latest states of models
            if obj.bIMMActivated_
                obj.stateCVbak_ = obj.stateEstCV_;
            end
            obj.stateEstbak_ = obj.stateEst_;
            
            % update model likelihood
            if obj.bIMMActivated_
                %%% update model likelihoods
                obj.updateIMMlikelihoods();
                %%% combine states of models
                obj.combinationIMMstateEstimates(dt);
                %%% repeat step 1: mixing state estimates
                obj.computeIMMmixingStates(dt);
            end              
            
            if obj.bIMMActivated_ && obj.mapIMMmodelProb_('CV') >= obj.mapIMMmodelProb_('CTRV')
                if ~bSecondModel
                    fprintf('switched to %s motion model ...\n','CV')
                    bSecondModel = true;
                end
                obj.bIMMCVModelUsed_ = true;
            else
                if bSecondModel
                    fprintf('switched to %s motion model ...\n','CTRV')
                    bSecondModel = false;
                end
                obj.bIMMCVModelUsed_ = false;
            end 
            
            if obj.bIMMActivated_ && obj.bIMMCVModelUsed_
                tmpOmega = (obj.stateEstCV_(3) - obj.psiCVbak_) / dt;
                stateTmpOut = [obj.stateEstCV_(1:4); tmpOmega]; %obj.stateEst_(5)];
%                 stateTmpOut = obj.stateEstCV_;
            else
                stateTmpOut = obj.stateEst_;
            end
            
            if obj.bIMMActivated_
                obj.psiCVbak_ = obj.stateEstCV_(3);
            end
            
            if nargout == 1
                varargout{1} = stateTmpOut;
            else
                varargout = [];
            end
            
            if length(obj.dequeVelTracker_) >= obj.NWindow
                obj.dequeVelTracker_ = [obj.dequeVelTracker_(2:end), ...
                            stateTmpOut(end)];
            else
                obj.dequeVelTracker_ = [obj.dequeVelTracker_, ....
                            stateTmpOut(end)];
            end
            
            %%% update residual distribution
            obj.updateResiduals();
            
            %% perform fault detection with GLR
            %%% GLR CV 
            if obj.bIMMActivated_
                %for fault = obj.GLR_FTYPE
                if ~obj.bDetThshldGLR_
                    obj.GLRInstCV_.faultDetection('-residual',obj.yCV_,...
                                            '-errCov',obj.SCV_,...
                                            '-stateTransition',obj.JfCV,...
                                            '-measTransition',obj.HCV_,...
                                            '-kalmanGain',obj.KCV_,...
                                            '-fault',obj.GLR_FTYPE);
                else
                    obj.GLRInstCV_.faultDetection('-residual',obj.stateEstCV_(1:2) - dRef,...
                                            '-errCov',obj.SCV_,...
                                            '-stateTransition',obj.JfCV,...
                                            '-measTransition',obj.HCV_,...
                                            '-kalmanGain',obj.KCV_,...
                                            '-fault',obj.GLR_FTYPE);
                end
                % end
            end
            %%% GLR CTRV
            if ~obj.bDetThshldGLR_
                obj.GLRInst_.faultDetection('-residual',obj.y_,...
                                        '-errCov',obj.S_,...
                                        '-stateTransition',obj.Jf,...
                                        '-measTransition',obj.H_,...
                                        '-kalmanGain',obj.K_,...
                                        '-fault',obj.GLR_FTYPE);
            else
                obj.GLRInst_.faultDetection('-residual',obj.stateEst_(1:2) - dRef,...
                                        '-errCov',obj.S_,...
                                        '-stateTransition',obj.Jf,...
                                        '-measTransition',obj.H_,...
                                        '-kalmanGain',obj.K_,...
                                        '-fault',obj.GLR_FTYPE);
            end
            % end
            
            %% reset position estimates to (0,0)
            %%% since prediction-correction cycle is always performed to
            %%% latest state estimate (in ECEF)          
            if obj.bIMMActivated_
                obj.stateEstCV_(1:2) = [0;0];
            end
            obj.stateEst_(1:2) = [0;0];
            
        end
        
        function updateResiduals(obj)
            %%% updateResiduals(obj,residual) computes the squared residual
            %%% and norms it to the innovation covariance matrix, 
            %%% contained in obj
                        
            if obj.bIMMActivated_ && obj.bIMMCVModelUsed_
                residual = obj.yCV_;
                S = obj.SCV_;
            else
                residual = obj.y_;
                S = obj.S_;
            end
            
            % compute squared and normed residual 
            rSquaredNormed = residual' * (S \ residual);
            % rNormed = residual' /  sqrt(obj.S_);
            % rSquaredNormed = val^2 / obj.Rz;
            % rSquaredNormed = val;
            
            obj.dequeResidualsSquaredNorm = [obj.dequeResidualsSquaredNorm; rSquaredNormed];
            % obj.dequeResidualsNormed = [obj.dequeResidualsNormed(:,2:end); rNormed'];
            obj.dequeResiduals = [obj.dequeResiduals, residual];
            
            % append to deque
%             if ~isempty(find(isinf(obj.dequeResidualsSquaredNorm),1))
%                 idx = find(isinf(obj.dequeResidualsSquaredNorm),1);
%                 obj.dequeResidualsSquaredNorm(idx) = rSquaredNormed;
%                 %obj.dequeResidualsNormed(:,find(isinf(obj.dequeResidualsSquaredNorm),1)) = rNormed';
%                 obj.dequeResiduals(:,idx) = residual;
%             else
%                 %%% pop front
%                 obj.dequeResidualsSquaredNorm = [obj.dequeResidualsSquaredNorm(2:end); rSquaredNormed];
%                 % obj.dequeResidualsNormed = [obj.dequeResidualsNormed(:,2:end); rNormed'];
%                 obj.dequeResiduals = [obj.dequeResiduals(:,2:end), residual];
%             end
            
        end
        
    end
    
    %%%% 
    methods              
        function dRet = isFailure_ResidualHypothesisTesting(obj)
            %%% the failure test relies on three measures:
            %%%  1. chi-square testing
            %%%  2. boundary test (2*sigma) ... not used!
            %%%  3. autocorrelation test
            
            %%% chi-square table testing
            idx = length(obj.dequeResidualsSquaredNorm);
            if idx >= obj.NWindow
                q_m = mean(obj.dequeResidualsSquaredNorm(idx-obj.NWindow+1:idx));
%                 dRet = ~(obj.NWindow * q_m > ...
%                         obj.ChiSquareConfPerc0p975(length(obj.y_)*obj.NWindow) & ...
%                        obj.NWindow * q_m < ...
%                         obj.ChiSquareConfPerc0p025(length(obj.y_)*obj.NWindow));
%                 dRet = ~ (obj.NWindow * q_m < ...
%                          obj.ChiSquareConfPerc0p050((length(obj.y_)-1)*obj.NWindow));
                if obj.bIMMCVModelUsed_
                    dRet = ~ (obj.NWindow * q_m < ...
                         obj.ChiSquareConfPerc0p050((length(obj.yCV_))*obj.NWindow));
                else
                    dRet = ~ (obj.NWindow * q_m < ...
                         obj.ChiSquareConfPerc0p050((length(obj.y_))*obj.NWindow));
                end
            else
                NTemp = idx;
                q_m = mean(obj.dequeResidualsSquaredNorm(1:idx));
%                 dRet = ~(NTemp * q_m > ...
%                         obj.ChiSquareConfPerc0p975(length(obj.y_)*NTemp) & ...
%                        NTemp * q_m < ...
%                         obj.ChiSquareConfPerc0p025(length(obj.y_)*NTemp));
                if obj.bIMMCVModelUsed_
                    dRet = ~ (NTemp * q_m < ...
                        obj.ChiSquareConfPerc0p050(length(obj.yCV_)*NTemp));
                else
                    dRet = ~ (NTemp * q_m < ...
                        obj.ChiSquareConfPerc0p050(length(obj.y_)*NTemp));
                end
            end
            
            %%% chi-square table testing
%             idx = find(~isinf(obj.dequeResidualsSquaredNorm),1,'last');
%             if idx > obj.NWindow
%                 q_m = mean(obj.dequeResidualsSquaredNorm(idx-obj.NWindow+1:idx));
% %                 dRet = ~(obj.NWindow * q_m > ...
% %                         obj.ChiSquareConfPerc0p975(length(obj.y_)*obj.NWindow) & ...
% %                        obj.NWindow * q_m < ...
% %                         obj.ChiSquareConfPerc0p025(length(obj.y_)*obj.NWindow));
%                 dRet = ~ (obj.NWindow * q_m < ...
%                          obj.ChiSquareConfPerc0p050((length(obj.y_)-1)*obj.NWindow));
%             else
%                 NTemp = idx;
%                 q_m = mean(obj.dequeResidualsSquaredNorm(1:idx));
% %                 dRet = ~(NTemp * q_m > ...
% %                         obj.ChiSquareConfPerc0p975(length(obj.y_)*NTemp) & ...
% %                        NTemp * q_m < ...
% %                         obj.ChiSquareConfPerc0p025(length(obj.y_)*NTemp));
%                 dRet = ~ (NTemp * q_m < ...
%                      obj.ChiSquareConfPerc0p050((length(obj.y_)-1)*NTemp));
%             end
            
            %%% autocorrelation of residual
%             xCorr = obj.acorr(obj.dequeResiduals(:,~isinf(obj.dequeResiduals(1,:))));
%             % get length of stored residual samples
%             % N = size(obj.dequeResiduals(:,~isinf(obj.dequeResiduals(1,:))),2);
%             
%             xCorr = xCorr ./ repmat(max(xCorr,[],2),[1 size(xCorr,2)]);
%             
%             % get indexes in 2sigma boundary (property of autocorrelation function that variance is 1/N)
%             idx = abs(xCorr) < 2/sqrt(numel(xCorr));
%             
%             %%% compute final isFailure flag by OR'ing Chi-Square and
%             %%% autocorrelation method
%             dRet = dRet & length(find(idx))/(numel(xCorr)) < 0.95;
        end
        
        function bRet = isFailure_ResidualHypothesisTestingGLR(obj,model)
            
%             bOut = containers.Map({'CV','CTRV'},{[],[]});
%             if obj.bIMMActivated_ 
%                 bOut('CV') = obj.GLRInstCV_.hasFailureDetected();
%             end
%             bOut('CTRV') = obj.GLRInst_.hasFailureDetected();
%             
%             if obj.bIMMActivated_ && obj.bIMMCVModelUsed_
%                 bRet = bOut('CV');
%             else
%                 bRet = bOut('CTRV');
%             end

            if obj.bIMMActivated_
                if strcmp(model,'CV')
                    bRet = obj.GLRInstCV_.hasFailureDetected();
                elseif strcmp(model,'CTRV')
                    bRet = obj.GLRInst_.hasFailureDetected();
                else
                    mExc = MException('EKF:isFailure_ResidualHypothesisTestingGLR','unknown model requested %s',model);
                    throw(mExc)
                end
            else
                bRet = obj.GLRInst_.hasFailureDetected();
            end
            
        end            
        
        function checkConstraints(obj)
            
            % obj.stateEst_(3) = mod(obj.stateEst_(3),2*pi);
            
            obj.stateEst_(obj.mapPModel2VelInState_('CTRV')) = abs(obj.stateEst_(obj.mapPModel2VelInState_('CTRV')));
            if obj.bIMMActivated_
                obj.stateEstCV_(obj.mapPModel2VelInState_('CV')) = abs(obj.stateEstCV_(obj.mapPModel2VelInState_('CV')));
            end
        end
        
    end
    
    %%%% GETTER METHODS
    methods 
        
        function stat = getTestStatisticsGLR(obj,model)
            
%             if obj.bIMMCVModelUsed_
%                 tmpGLR = obj.GLRInstCV_;
%             else
%                 tmpGLR = obj.GLRInst_;
%             end
            if obj.bIMMActivated_
                if strcmp(model,'CV')
                    tmpGLR = obj.GLRInstCV_;
                elseif strcmp(model,'CTRV')
                    tmpGLR = obj.GLRInst_;
                else
                    mExc = MException('EKF:getTestStatisticsGLR','unknown model requested %s',model);
                    throw(mExc)
                end
            else
                tmpGLR = obj.GLRInst_;
            end
            
            tmpMap = containers.Map(tmpGLR.hashStatistics_.keys(),cell(1,length(tmpGLR.hashStatistics_.keys())));
            
            for fType = tmpGLR.hashStatistics_.keys()
                tmpMap(fType{1}) = tmpGLR.hashStatistics_(fType{1}).lt_bak;
            end
            
            stat = tmpMap;
        end
        
        function kRet = getTimeOfFailureGLR(obj,model)
            if obj.bIMMActivated_
                if strcmp(model,'CV')
                    kRet = obj.GLRInstCV_.getInstanceOfFailure();
                else
                    kRet = obj.GLRInst_.getInstanceOfFailure();
                end
            else
                kRet = obj.GLRInst_.getInstanceOfFailure();
            end
        end
        
        function trNormed = getTraceNormalized(obj)
            tmp = reshape(obj.tracePredP_,[length(obj.stateEst_) length(obj.stateEst_) size(obj.tracePredP_,2)/length(obj.stateEst_)]);
            tmpNorm = tmp ./ repmat(max(abs(tmp),[],3),[1 1 size(tmp,3)]);
            trNormed = zeros(1,size(tmpNorm,3));
            for n=1:size(tmpNorm,3)
                trNormed(n) = trace(squeeze(tmpNorm(:,:,n))) / numel(obj.stateEst_);
            end
            
            if obj.bIMMCVModelUsed_
                tmp = reshape(obj.tracePredPCV_,[length(obj.stateEstCV_) length(obj.stateEstCV_) size(obj.tracePredPCV_,2)/length(obj.stateEstCV_)]);
                tmpNorm = tmp ./ repmat(max(abs(tmp),[],3),[1 1 size(tmp,3)]);
                trNormed = zeros(1,size(tmpNorm,3));
                for n=1:size(tmpNorm,3)
                    trNormed(n) = trace(squeeze(tmpNorm(:,:,n))) / numel(obj.stateEstCV_);
                end
            end
        end
        
        function [DRet,LagsRet] = getResidualsDistribution(obj)
            [DRet,LagsRet] = hist(obj.dequeResidualsSquaredNorm(~isinf(obj.dequeResidualsSquaredNorm)),...
                length(obj.dequeResidualsSquaredNorm(~isinf(obj.dequeResidualsSquaredNorm))));
        end
        
        function dRet = getResidualsSNorm(obj)
            dRet = obj.dequeResidualsSquaredNorm(~isinf(obj.dequeResidualsSquaredNorm));
        end
        
        function dRet = getAllResiduals(obj)
            dRet = obj.dequeResiduals(:,~isinf(obj.dequeResidualsSquaredNorm));
        end
        
        function dRet = getNumelState(obj)
            dRet = length(obj.stateEst_);
        end
              
        function dRet = getNumelStateFB(obj)
            dRet = length(obj.stateEstCV_);
        end
        
        function dCov = getSCov(obj,model)
            if obj.bIMMActivated_ && nargin==2
                if strcmp(model,'CV')
                    dCov = obj.SCV_;                
                else
                    dCov = obj.S_;                
                end
            else
                dCov = obj.S_;
            end
        end
        
        function dGain = getKGain(obj,model)
            if obj.bIMMActivated_ && nargin==2
                if strcmp(model,'CV')
                    dGain = obj.KCV_;
                else
                    dGain = obj.K_;
                end
            else
                dGain = obj.K_;
            end
        end
        
        function dTrans = getStateTransition(obj,model)
            if obj.bIMMActivated_ && nargin==2
                if strcmp(model,'CV')
                    dTrans = obj.JfCV;
                else
                    dTrans = obj.Jf;
                end
            else
                dTrans = obj.Jf;
            end
        end
        
        function dTrans = getMeasTransition(obj,model)
            if obj.bIMMActivated_ && exist('model','var')
                if strcmp(model,'CV')
                    dTrans = obj.HCV_;
                else
                    dTrans = obj.H_;
                end
            else
                dTrans = obj.H_;
            end
        end
        
        function dState = getEKFStateEstimateCV(obj)
            if obj.bIMMActivated_
                dState = obj.stateCVbak_;
            else 
                fprintf('[WARNING] Currently there is no CV model running in parallel\n')
                dState = [];
            end
        end
                
        function dState = getEKFStateEstimateCTRV(obj)
            dState = obj.stateEstbak_;
        end
        
        function dState = getEKFStateEstimate(obj)
            dState = obj.stateEstbak_;
        end
        
        function dState = getIMMStateEstimate(obj)
            dState = obj.dIMMstateCombined_;
        end
        
        function dPCov = getIMMPCovEstimate(obj)
            dPCov = obj.dIMMerrCovP_;
        end
        
        function dProbs = getIMMmodelProbs(obj)
            dProbs = obj.mapIMMmodelProb_;
        end
        
        function dMllhs = getIMMModelLikelihoods(obj)
            dMllhs = obj.dIMMllhood_;
        end
        
        function bRet = isSlowMoving(obj)
            bRet = mean(obj.dequeVelTracker_) < obj.VEL_MIN_THSLD;
        end
        
        function pop_back_pcov(obj)
            if ~isempty(obj.tracePredP_)
                obj.tracePredP_(:,end-length(obj.stateEst_)+1:end) = [];
            end
            
            if obj.bIMMActivated_
                if ~isempty(obj.tracePredPCV_)
                    obj.tracePredPCV_(:,end-length(obj.stateEstCV_)+1:end) = [];
                end
            end
           
        end
        
        function dFail = evalPredCov(obj)
            
            if obj.bIMMActivated_ && obj.bIMMCVModelUsed_
                N = length(obj.stateEstCV_);
                tmpP = reshape(obj.tracePredPCV_,[N N size(obj.tracePredPCV_,2)/N]);
            else
                N = length(obj.stateEst_);
                tmpP = reshape(obj.tracePredP_,[N N size(obj.tracePredP_,2)/N]);
            end
            
            if size(tmpP,3) < obj.NWindow
                dFail = 0;
                return
            end
            
            mPCov = mean(tmpP(:,:,end-obj.NWindow +1:end),3);
            
            tPVars = diag(mPCov);
            
            %%% check IMO specs for safety region?
            % breadth = 30; % use as sigma
            
            %%% TODO: apply Fuji bounds and rotate according to current COG
            %%% to check for error bounds violation
            
            covPe = tPVars(1);
            covPn = tPVars(2);
            
            %%% counter clockwise rotation
            theta = pi / 2 - obj.stateEst_(3);
            covBody = [cos(theta) -sin(theta); sin(theta) cos(theta)] * [covPe; covPn];
                                 
%             if mean(obj.dequeVelTracker_) < obj.VEL_MIN_THSLD
% %                 errP = 5 * breadth; % 3 sigma boundary
% %                 errCOG = 25/180*pi; % COG error bound
% %                 errROT = 5/180*pi;
% %                 errSOG = 10;
%                 
%                 errP = 5 * breadth; % 3 sigma boundary
%                 errCOG = 25;%/180*pi; % COG error bound
%                 errROT = .25;%/180*pi;
%                 errSOG = 1; 
%             else
%                 errP = 3 * breadth; % 3 sigma boundary
%                 errCOG = 7;%/180*pi; % COG error bound
%                 errROT = .1;%/180*pi;
%                 errSOG = 1; 
%             end
            if mean(obj.dequeVelTracker_) < obj.VEL_MIN_THSLD
                errPx = 5 * obj.dShipL_;% errP*sin(obj.stateEst_(3));
                errPy = 5 * obj.dShipL_;% errP*cos(obj.stateEst_(3));
            else
                errPx = 4 * obj.dShipL_;% errP*sin(obj.stateEst_(3));
                errPy = 1.6 * obj.dShipL_;% errP*cos(obj.stateEst_(3));
            end
            errIDX = ( abs([covBody(2) covBody(1)]) - [errPx errPy].^2) > 0;
%             if obj.bIMMCVModelUsed_
%                 errIDX = ( [covPx covPy] - [errPx errPy].^2) > 0;
%             else
%                 errIDX = ( [covPx covPy] - [errPx errPy].^2) > 0;
%             end
            
            % obj.stateEst_(4)/2; % SOG err bound (linked to current SOG)
            % errTHDG = 10 / 180 * pi; 
%             if obj.bIMMCVModelUsed_
%                 errIDX = ( tPVars' - [errP errP errCOG errSOG].^2) > 0;
%             else
%                 errIDX = ( tPVars' - [errP errP errCOG errROT errSOG].^2) > 0;
%             end
            dFail = ~isempty(find(errIDX,1));
%             tmpPNormed = abs(tmpP)./repmat(max(abs(tmpP),[],3),[1 1 size(tmpP,3)]);
%             
%             trNormed = zeros(1,size(tmpP,3));        
%             
%             for n=1:size(tmpP,3)
%                 trNormed(n) = trace(tmpPNormed);
%             end
%             
%             trNormed = trNormed ./ max(trNormed);
%             
%             dFail = trNormed(end);
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Interacting Multiple Model %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function computeIMMmixingStates(obj,dt)
            % x1 = [x1c; x2e]; 
            obj.mapIMMstate2model_('CV') = [obj.stateEstCV_; (obj.stateEstCV_(3) - obj.psiCVbak_)/dt] ;
            % x2 = [x2c; x2e];
            obj.mapIMMstate2model_('CTRV') = [obj.stateEst_(1:4); obj.stateEst_(5)];
            % P1c = P1
            obj.mapIMMPCommon2model_('CV') = obj.PCV_;
            % P2c = P2c
            obj.mapIMMPCommon2model_('CTRV') = obj.P_(1:4,1:4);
            % P1e = P2e
            obj.mapIMMPExtra2model_('CV') = obj.P_(5,5);
            obj.mapIMMPExtra2model_('CTRV') = obj.P_(5,5);
            % P1ce = 0, P1ec = 0 (vector ...)
            obj.mapIMMPce2model_('CV') = zeros(length(obj.stateEstCV_),1);
            obj.mapIMMPec2model_('CV') = zeros(1,length(obj.stateEstCV_));
            % P2ce = cov(x2c,x2e), P2ec = cov(x2e,x2c)
            obj.mapIMMPce2model_('CTRV') = obj.P_(1:4,5);
            obj.mapIMMPec2model_('CTRV') = obj.P_(5,1:4);
            
            N = numel(obj.mapIMMstate2model_.keys());
            %% step 1: mixing state estimates
            for modelj=obj.mapIMMstate2model_.keys()
                obj.mapIMMnormC_(modelj{1}) = 0;
                %%% compute normalization constant for model j
                for modeli=obj.mapIMMstate2model_.keys()
                    obj.mapIMMnormC_(modelj{1}) = obj.mapIMMnormC_(modelj{1}) + ...
                        obj.dIMMtransProb_(obj.mapIMMmodel2idx_(modeli{1}),obj.mapIMMmodel2idx_(modelj{1})) * obj.mapIMMmodelProb_(modeli{1});
                end
                
                obj.mapIMMaposterioriState_(modelj{1}) = 0;
                obj.mapIMMaposterioriErrCovP_(modelj{1}) = 0;
                
                %%% compute mixing weights (i|j) and a priori states for 
                %%% model j, a priori/initial err Cov P for model j is also
                %%% estimated
                for modeli=obj.mapIMMstate2model_.keys()
                    obj.dIMMmixingWeights_(obj.mapIMMmodel2idx_(modeli{1}),obj.mapIMMmodel2idx_(modelj{1})) = ...
                        1 / obj.mapIMMnormC_(modelj{1}) * obj.dIMMtransProb_(obj.mapIMMmodel2idx_(modeli{1}),obj.mapIMMmodel2idx_(modelj{1})) * ...
                        obj.mapIMMmodelProb_(modeli{1});
                    
                                    
                end               
            end
            
            %% initial mixing of state (unbiased)
            %%% CV
            obj.mapIMMaposterioriState_('CV') = obj.stateEstCV_ * obj.dIMMmixingWeights_(obj.mapIMMmodel2idx_('CV'),obj.mapIMMmodel2idx_('CV')) + ...
                 obj.stateEst_(1:4) * obj.dIMMmixingWeights_(obj.mapIMMmodel2idx_('CTRV'),obj.mapIMMmodel2idx_('CV')); 
            %%% CTRV
            obj.mapIMMaposterioriState_('CTRV') = obj.stateEst_ * obj.dIMMmixingWeights_(obj.mapIMMmodel2idx_('CTRV'),obj.mapIMMmodel2idx_('CTRV')) + ...
                 obj.mapIMMstate2model_('CV') * obj.dIMMmixingWeights_(obj.mapIMMmodel2idx_('CV'),obj.mapIMMmodel2idx_('CTRV')); 
            %% initial mixing of error covariance (unbiased)
            %%% CV
            tmpP11 = (obj.stateEstCV_ - obj.mapIMMaposterioriState_('CV'))*(obj.stateEstCV_ - obj.mapIMMaposterioriState_('CV'))';
            tmpP21 = (obj.stateEst_(1:4) - obj.mapIMMaposterioriState_('CV'))*(obj.stateEst_(1:4) - obj.mapIMMaposterioriState_('CV'))';
            obj.mapIMMaposterioriErrCovP_('CV') = ...
                (obj.mapIMMPCommon2model_('CV') + tmpP11) * obj.dIMMmixingWeights_(obj.mapIMMmodel2idx_('CV'),obj.mapIMMmodel2idx_('CV')) + ...
                (obj.mapIMMPCommon2model_('CTRV') + tmpP21) * obj.dIMMmixingWeights_(obj.mapIMMmodel2idx_('CTRV'),obj.mapIMMmodel2idx_('CV'));
            %%% CTRV
            aprioriStateCV = obj.mapIMMstate2model_('CV');
            aposterioriStateCTRV = obj.mapIMMaposterioriState_('CTRV');
            % compute intermediate covariances to mitigate bias in cov
            % estimate
            tmpP12c = (obj.stateEstCV_ - aposterioriStateCTRV(1:4))*(obj.stateEstCV_ - aposterioriStateCTRV(1:4))';
            tmpP12ce = (obj.stateEstCV_ - aposterioriStateCTRV(1:4))*(aprioriStateCV(5) - aposterioriStateCTRV(5))';
            tmpP12ec = (aprioriStateCV(5) - aposterioriStateCTRV(5))*(obj.stateEstCV_ - aposterioriStateCTRV(1:4))';
            tmpP12e = (aprioriStateCV(5) - aposterioriStateCTRV(5))*(aprioriStateCV(5) - aposterioriStateCTRV(5))';
            
            tmpP22c = (obj.stateEst_(1:4) - aposterioriStateCTRV(1:4))*(obj.stateEst_(1:4) - aposterioriStateCTRV(1:4))';
            tmpP22ce = (obj.stateEst_(1:4) - aposterioriStateCTRV(1:4))*(obj.stateEst_(5) - aposterioriStateCTRV(5))';
            tmpP22ec = (obj.stateEst_(5) - aposterioriStateCTRV(5))*(obj.stateEst_(1:4) - aposterioriStateCTRV(1:4))';
            tmpP22e = (obj.stateEst_(5) - aposterioriStateCTRV(5))*(obj.stateEst_(5) - aposterioriStateCTRV(5))';
            
            obj.mapIMMaposterioriErrCovP_('CTRV') = ...
                ([obj.mapIMMPCommon2model_('CV') obj.mapIMMPce2model_('CV'); obj.mapIMMPec2model_('CV') obj.mapIMMPExtra2model_('CV')] + ...
                [tmpP12c tmpP12ce; tmpP12ec tmpP12e]) * obj.dIMMmixingWeights_(obj.mapIMMmodel2idx_('CV'),obj.mapIMMmodel2idx_('CTRV')) + ...
                ([obj.mapIMMPCommon2model_('CTRV') obj.mapIMMPce2model_('CTRV'); obj.mapIMMPec2model_('CTRV') obj.mapIMMPExtra2model_('CTRV')] + ...
                [tmpP22c tmpP22ce; tmpP22ec tmpP22e]) * obj.dIMMmixingWeights_(obj.mapIMMmodel2idx_('CTRV'),obj.mapIMMmodel2idx_('CTRV'));
                        
            %% this would conclude the original idea of IMM step 1) mixed initial conditions
            obj.stateEstCV_ = obj.mapIMMaposterioriState_('CV');
%             obj.stateEstCV_ = obj.stateEstCV_([1:3,5]);
            obj.PCV_ = obj.mapIMMaposterioriErrCovP_('CV');
            % obj.PCV_ = obj.PCV_(1:4,1:4);
            obj.stateEst_ = obj.mapIMMaposterioriState_('CTRV');
            obj.P_ = obj.mapIMMaposterioriErrCovP_('CTRV');
%             tmp = obj.mapIMMaposterioriState_('CTRV');
%             obj.stateEst_([1:3,5]) =  tmp;
%             tmp = obj.mapIMMaposterioriErrCovP_('CTRV');
%             obj.P_([1:3,5],[1:3,5]) = tmp;
            
        end
        
        function updateIMMlikelihoods(obj)
            %% step 2: compute likelihood of models (from Gaussian)
            for model = obj.dIMMllhood_.keys()
                if strcmp(model{1},'CV')
                    obj.dIMMllhood_(model{1}) = 1 / sqrt(det(2*pi*obj.SCV_)) * exp(-0.5 * obj.yCV_' * (obj.SCV_ \ obj.yCV_));
                else
                    obj.dIMMllhood_(model{1}) = 1 / sqrt(det(2*pi*obj.S_)) * exp(-0.5 * obj.y_' * (obj.S_ \ obj.y_));
                end
            end
            
            for model1 = obj.dIMMllhood_.keys()
                llh = obj.dIMMllhood_(model1{1});
                cNorm = 0;
                for model2 = obj.dIMMllhood_.keys()
                    cNorm = cNorm + obj.dIMMllhood_(model2{1}) * obj.mapIMMnormC_(model2{1});
                end
                obj.mapIMMmodelProb_(model1{1}) = llh * obj.mapIMMnormC_(model1{1}) / cNorm;
            end
        end
        
        function combinationIMMstateEstimates(obj,dt)
            obj.mapIMMstate2model_('CV') = [obj.stateEstCV_; (obj.stateEstCV_(3) - obj.psiCVbak_)/dt] ;
            % x2 = [x2c; x2e];
            obj.mapIMMstate2model_('CTRV') = [obj.stateEst_(1:4); obj.stateEst_(5)];
            % P1c = P1
            obj.mapIMMPCommon2model_('CV') = obj.PCV_;
            % P2c = P2c
            obj.mapIMMPCommon2model_('CTRV') = obj.P_(1:4,1:4);
            % P1e = P2e
            obj.mapIMMPExtra2model_('CV') = obj.P_(5,5);
            obj.mapIMMPExtra2model_('CTRV') = obj.P_(5,5);
            % P1ce = 0, P1ec = 0 (vector ...)
            obj.mapIMMPce2model_('CV') = zeros(length(obj.stateEstCV_),1);
            obj.mapIMMPec2model_('CV') = zeros(1,length(obj.stateEstCV_));
            % P2ce = cov(x2c,x2e), P2ec = cov(x2e,x2c)
            obj.mapIMMPce2model_('CTRV') = obj.P_(1:4,5);
            obj.mapIMMPec2model_('CTRV') = obj.P_(5,1:4);
            %% combination of all means (of state) estimated by each model
            obj.dIMMstateCombined_ = 0;
                       
            for modeli = obj.mapIMMstate2model_.keys()
                obj.dIMMstateCombined_ = obj.dIMMstateCombined_ + ...
                    obj.mapIMMstate2model_(modeli{1}) * obj.mapIMMmodelProb_(modeli{1});
            end
            
            %% combination of all associated covariances (of state) estimated by each model
            obj.dIMMerrCovP_ = 0;
            
            for modeli = obj.mapIMMstate2model_.keys()
                obj.dIMMerrCovP_ = obj.dIMMerrCovP_ + ...
                    obj.mapIMMmodelProb_(modeli{1}) * ...
                    ([obj.mapIMMPCommon2model_(modeli{1}) obj.mapIMMPce2model_(modeli{1});
                      obj.mapIMMPec2model_(modeli{1}) obj.mapIMMPExtra2model_(modeli{1})] + ...
                    (obj.mapIMMstate2model_(modeli{1}) - obj.dIMMstateCombined_) * (obj.mapIMMstate2model_(modeli{1}) - obj.dIMMstateCombined_)');
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% LISTENER IMPLEMENTATION %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Static)
        function eventListener__onObservation(obj,~,evt)
            %%% event is triggered on new observation
            
            if obj.status_ ~= AISIntegrityBasePack.AISBaseStates.Idle && ...
                    obj.status_ ~= AISIntegrityBasePack.AISBaseStates.Finished
                fprintf('[INFO] %s is currently busy, waiting to continue at %s \n',class(obj),obj.getNameFromStack);
                waitfor(obj.getID(),'status_',AISIntegrityBasePack.AISBaseStates.Finished)
            end
            
            %%% set internal status to active
            obj.status_ = AISIntegrityBasePack.AISBaseStates.Active;
            
            %%% update state with new measurement
            obj.update(evt.obs);
            
            %%% evaluate failure 
            obj.curFailureObs_ = obj.isFailure_ResidualHypothesisTesting();
            obj.traceFailure_ = [obj.traceFailure_ obj.curFailureObs_];
            
            %%% set internal state back to finished
            obj.status_ = AISIntegrityBasePack.AISBaseStates.Finished;
            
            if obj.plot_
                obj.plotIntermediates('update',evt.dt)
            end
        end
        
        function eventListener__onPrediction(obj,~,evt)
            %%% event is triggered on new observation
            
            fprintf('event listener activated at %s\n',...
                datestr(now,'dd-mmm-yyyy HH:MM:SS.FFF'));
            
            if obj.status_ ~= AISIntegrityBasePack.AISBaseStates.Idle && ...
                    obj.status_ ~= AISIntegrityBasePack.AISBaseStates.Finished
                fprintf('[INFO] %s is currently busy, waiting to continue at %s \n',class(obj),obj.getNameFromStack);
                waitfor(obj.getID(),'status_',AISIntegrityBasePack.AISBaseStates.Finished)
            end
            
            %%% set internal status to active
            obj.status_ = AISIntegrityBasePack.AISBaseStates.Active;
            
            %%% update state with new measurement
            obj.prediction(evt.dt,evt.controlIn);
            
            %%% evaluate predicted covariance
            tmpFailPred = obj.evalPredCov();
            obj.curFailurePred_ = tmpFailPred;
            
            %%% set internal state back to finished
            obj.status_ = AISIntegrityBasePack.AISBaseStates.Finished;
            
            if obj.plot_
                obj.plotIntermediates('prediction',evt.dt)
            end
        end
        
        function xCorr = acorr(xIn)
            
            xCorr = zeros(size(xIn,1),2*size(xIn,2)-1);
            
            N = size(xIn,2);
            
            for n=1:size(xIn,1)
                %FFT method based on zero padding
                fx = fft([xIn(n,:)'; zeros(N,1)]); % zero pad and FFT
                xCorrTmp = ifft(fx.*conj(fx)); % abs()^2 and IFFT
                % circulate to get the peak in the middle and drop one
                % excess zero to get to 2*n-1 samples
                xCorr(n,:) = [xCorrTmp(N+2:end)', xCorrTmp(1:N)'];
            end
        end
        
        
    end
end