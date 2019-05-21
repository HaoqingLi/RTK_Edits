classdef VGeneralPreparation_ExtendedKalmanFilterRTK < TutorialBasePack.AbstractKalmanFilter
    properties
        F_ % Jacobian of f (derivative w.r.t. state)
        C_ % Jacobian of f (derivative w.r.t. noise)
        H_ % Jacobian of h 
        satellitesLOS_ % vector with the satellites PRN in LOS at the moment
        satRef_ % reference satellite for the doble difference approach for RTK
        basePosition_ % position of the base or reference station, from which the position is known
        sizeState_ % number of elements which are not ambiguities in the state (in case you have position, velocity, attitude, etc etc... or only position, etc)
        WavelengthL1_ = 0.19029367;
        WavelengthL2_ = 0.24421021;
        DDAmb_ 
        SDAmb_ 
        D_
        QRN_
        QNR_
        QN_
        initialAmbCovariance_ = 100;
    end
    
    methods 
        function obj = VGeneralPreparation_ExtendedKalmanFilterRTK(varargin)
            
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
        
        function varargout = prediction(obj,dt)
            
            obj.state_ = obj.f_(dt,obj.state_);
            
            % prediction of state covariance
            Jf = obj.F_(dt,obj.state_);
            Jc = obj.C_(dt,obj.state_);
            %%%%%%%%%%%%%%
%             obj.P_ = Jf * obj.P_ * Jf' + Jc * obj.Q_ * Jc';
            Q_velocity = 0.1;
            Q  = Q_velocity * [zeros(3), eye(3)*dt^2/2, zeros(3, length(obj.state_)-6);  zeros(3)*dt^2/2, eye(3)*dt,  zeros(3, length(obj.state_)-6); zeros(length(obj.state_)-6, length(obj.state_) ) ];
            obj.P_ = Jf * obj.P_ * Jf' + Q;
            %%%%%%%%%%%%%%
                       
            if nargout == 1
                varargout{1} = obj.state_;
            else
                varargout = [];
            end
        end
        
        function varargout = correctionOriginal(obj,sObs)       
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
        
        function varargout = correction(obj,varargin)       
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
            
            if isempty(obj.satRef_)                                         % Check whether we have a reference satellite for the Double Difference observation approach
                obj.satRef_         = satRef;
                obj.satellitesLOS_  = [];
            elseif satRef ~= obj.satRef_  
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
            for iC=1:length(difNew),                                        % Initialize the covariance for the phase ambiguities
                newP(obj.sizeState_+difNew(iC),obj.sizeState_+difNew(iC)) = obj.initialAmbCovariance_; 
            end 
            obj.P_              = newP;                                     % create the new covariance matrix, saving the relevant information from the former covariance matrix
            
            % Jacobian matrix of the solution
            [H, obj.D_]              = obj.H_( obj.state_, satPos, satRefPos, refPosIndx, waveLengVec, typeObs );
            z                   = [DDPhase; DDRange];                       % Pile up phase and range measurements           
            
            % Application of the correction model using the predicted state
            z([refPosIndx,nObs+refPosIndx]) = [];                           % Eliminate the measurements corresponding to the reference satellites
            satPos([refPosIndx],:) = [];                                   
            satRefPos([refPosIndx],:) = [];                                
            waveLengVec2 = waveLengVec; waveLengVec2(refPosIndx) = [];
            h_z                 = obj.h_( obj.state_(1:3)', satPos, satRefPos, obj.state_(obj.sizeState_+1:end), waveLengVec2, obj.D_); % Observation model
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

            % Saving the single and double difference phase ambiguities 
            obj.SDAmb_          = obj.state_(obj.sizeState_+1:end);         % Saving the single and double difference ambiguities of the satellites
            obj.DDAmb_          = obj.D_ * obj.SDAmb_;
            
            % Estimate the Q_RN, needed to find the fixed solution afterwards
            G_aux               = zeros(length(obj.state_)-max(typeObs),length(obj.state_));
            G_aux(1:obj.sizeState_,1:obj.sizeState_) = eye(obj.sizeState_);
            G_aux(obj.sizeState_+1:end,obj.sizeState_+1:end) = obj.D_;
            P_dd                = G_aux*obj.P_*G_aux';
            obj.QRN_            = P_dd(obj.sizeState_+1:end,1:obj.sizeState_);
            obj.QNR_            = P_dd(1:obj.sizeState_,obj.sizeState_+1:end);
            obj.QN_             = P_dd(obj.sizeState_+1:end,obj.sizeState_+1:end);
            
            %%%%%%%%%%% Trick used by Anja -> Change of the Covariance matrix
            Pauxx = obj.P_;
            Pauxx(1:obj.sizeState_,:) = 0;
            Pauxx(:,1:obj.sizeState_) = 0;
            Pauxx(1:obj.sizeState_,1:obj.sizeState_) = diag(diag(obj.P_(1:obj.sizeState_,1:obj.sizeState_))); 
            obj.P_ = Pauxx;
            %%%%%%%%%%%
            
            if nargout == 1
                varargout{1} = obj.state_;
            else
                varargout = [];
            end           
        end
        
        
        
        
        function varargout = correctionFixNHold(obj,varargin)       
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
            DDAmbFix    = varargin{9};
            
            if isempty(obj.satRef_)                                         % Check whether we have a reference satellite for the Double Difference observation approach
                obj.satRef_         = satRef;
                obj.satellitesLOS_  = [];
            elseif satRef ~= obj.satRef_  
                obj.satRef_         = satRef;
                obj.satellitesLOS_  = [];
            end
            [~,intNew,intOld]   = intersect(satLOS,obj.satellitesLOS_);     % The intersect of the old and new LOS satellites. tmp1 refers to the indexes of current observed satellites || tmp2 referes to the indexes of the formerly observed ones
            [~,difNew,difOld]   = setxor(satLOS,obj.satellitesLOS_);        % The intersect of the old and new LOS satellites. tmp1 refers to the indexes of current observed satellites || tmp2 referes to the indexes of the formerly observed ones
            obj.satellitesLOS_  = satLOS;                                   % Save in the class which is the PRN of the LOS satellites
            refPosIndx          = [1:1:length(satLOS)];                     % Vble to refer to the position in the vector of satellite in which the reference satellite is
            refPosIndx          = refPosIndx(ismember(satLOS,satRef));
            
            oldState            = obj.state_;                               % Copy of the previous state
            oldP                = obj.P_;                                   % Copy of the previous covariance
            oldP(:,obj.sizeState_+difOld) = [];  oldP(obj.sizeState_+difOld,:) = [];                  % Eliminate from the covariance related to the no-longer observed satellites
            
            initialAmbiguities  = (DDPhase - DDRange)'./waveLengVec;         % Approximation for the initial ambiguities (for satellites that we still do not have in our state)
            newState            = [oldState(1:obj.sizeState_); initialAmbiguities.*ones(size(satLOS))];   % Create a new state whose length is equal to the number of observed satellites + 6 (position and velocity)
            newState(obj.sizeState_+intNew)  = oldState(obj.sizeState_+intOld);  % Save the ambiguities of the satellites formerly observed and present in the previous state
            obj.state_          = newState;                                 % Use the current state as the state of the object
            
            newP                = augmentCovarianceMatrix(oldP, obj.sizeState_+difNew); 
            for iC=1:length(difNew),  
                newP(obj.sizeState_+difNew(iC),obj.sizeState_+difNew(iC)) = 100; 
            end 
            obj.P_              = newP;                                     % create the new covariance matrix, saving the relevant information from the former covariance matrix

            % compute innovation
            obj.y_              = ([DDPhase'; DDRange'] - obj.h_(obj.state_(1:3), obj.basePosition_, satPos, satRefPos, obj.state_(obj.sizeState_+1:end), waveLengVec));
            obj.y_([refPosIndx,nObs+refPosIndx]) = [];          % Eliminate the measurements related to the reference satellites
                       
            % Jacobian matrix of the solution
            [H, obj.D_]              = obj.H_( obj.state_, obj.basePosition_, satPos, satRefPos, refPosIndx, waveLengVec, typeObs );
            
            aux = [DDPhase'; DDRange'];
            aux([refPosIndx,nObs+refPosIndx]) = [];
            obj.y_ = aux - H*obj.state_;
            
            R_values            = diag(obj.R_);                             % Building the R matrix
            R_phase             = obj.D_*diag(R_values(1:nObs))*obj.D_';    % As we are having combination of measurements, we have to fuse also the variance from each satellite and the reference satellite
            R_code              = obj.D_*diag(R_values(nObs+1:end))*obj.D_';
            obj.R_              = eye((nObs-max(typeObs))*2); 
            obj.R_(1:nObs-max(typeObs),1:nObs-max(typeObs)) = R_phase; 
            obj.R_(nObs-max(typeObs)+1:end,nObs-max(typeObs)+1:end) = R_code;
            
            G = zeros([nObs-max(typeObs),length(obj.state_)]);
            G(:,obj.sizeState_+1:end) = obj.D_;
            R_virtualMeas = 0.001^2*eye(nObs-max(typeObs));
            R_big = zeros(3*(nObs-max(typeObs)) );
            R_big(1:2*(nObs-max(typeObs)),1:2*(nObs-max(typeObs))) = obj.R_;
            R_big(2*(nObs-max(typeObs))+1:end,2*(nObs-max(typeObs))+1:end) = R_virtualMeas;
            
            obj.R_ = R_big;
            H = [H;G];
            obj.y_ = [obj.y_; DDAmbFix - G*obj.state_];
            
            
            
            obj.S_              = H*obj.P_*H' + obj.R_;                     % Innovation covariance
            obj.K_              = obj.P_*H'/obj.S_;                         % Kalman Gain
            
            obj.state_          = obj.state_ + obj.K_ * obj.y_;             % update the mean of the state
            obj.P_              = obj.P_ - obj.K_ * H * obj.P_;             % update the state covariance.

            obj.SDAmb_          = obj.state_(obj.sizeState_+1:end);         % Saving the single and double difference ambiguities of the satellites
            obj.DDAmb_          = obj.D_ * obj.SDAmb_;
            
            G_aux               = zeros(length(obj.state_)-max(typeObs),length(obj.state_));
            G_aux(1:obj.sizeState_,1:obj.sizeState_) = eye(obj.sizeState_);
            G_aux(obj.sizeState_+1:end,obj.sizeState_+1:end) = obj.D_;
            Q_aux               = G_aux*obj.P_*G_aux';
            obj.QRN_            = Q_aux(1:obj.sizeState_,obj.sizeState_+1:end);
            obj.QN_             = Q_aux(obj.sizeState_+1:end,obj.sizeState_+1:end);
            
            
            %%%%%%%%%%%
            auxx                = ([DDPhase'; DDRange'] - obj.h_(obj.state_(1:3), obj.basePosition_, satPos, satRefPos, obj.state_(obj.sizeState_+1:end), waveLengVec));
            %%%%%%%%%%%
            
            
            if nargout == 1
                varargout{1} = obj.state_;
            else
                varargout = [];
            end           
        end
        
        
        
        
        %% Cycle Slip Detector
        function cycleSlipDetector(obj,newAmb, oldAmb, oldLOS)  
           
            % Error if the number of ambiguities is not equal to the number of observed satellites
            if length(newAmb) ~= length(obj.satellitesLOS_)-length(obj.satRef_)
                disp('EKF Cycle Slip: The number of estimated ambiguities does not match the number of satellites in LOS')
                return;
            end
            
            % Find the list of Satellites PRN excluding the reference satellites
            oldAux = []; newAux = [];
            for i=1:length(obj.satRef_)
               newAux = [newAux; find(obj.satellitesLOS_ == obj.satRef_(i))];
               oldAux = [oldAux; find(oldLOS             == obj.satRef_(i))];
            end
            newPRN_woRef = obj.satellitesLOS_;         newPRN_woRef(newAux) = [];
            oldPRN_woRef = oldLOS;                     oldPRN_woRef(oldAux) = [];
            
            resetAmb = [];
            if ~isempty(oldAmb)
                for i=1:length(newAmb)
                   if sum(newPRN_woRef(i) == oldPRN_woRef)
                       aux2 = find(newPRN_woRef(i) == oldPRN_woRef);
                       if newAmb(i) ~= oldAmb(aux2)
                           resetAmb = [resetAmb; newPRN_woRef(i)];
                           obj.satellitesLOS_(obj.satellitesLOS_==newPRN_woRef(i)) = 0;
                       end
                   end
                end
            end
            obj.satellitesLOS_(obj.satellitesLOS_==0) = []; 
            
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