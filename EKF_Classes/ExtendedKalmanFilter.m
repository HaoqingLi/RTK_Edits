classdef ExtendedKalmanFilter < TutorialBasePack.AbstractKalmanFilter
    properties
        F_ % Jacobian of f (derivative w.r.t. state)
        C_ % Jacobian of f (derivative w.r.t. noise)
        H_ % Jacobian of h 
    end
    
    methods 
        function obj = ExtendedKalmanFilter(varargin)
            
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
                            isa(varargin{2},'function_handle')
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
            %%% controlIn (i.e. u == control input) is expected to be a
            %%% struct containing rot and vel values (from AIS message)                 
            
            obj.state_ = obj.f_(dt,obj.state_);  
            
            % prediction of state covariance
            Jf = obj.F_(dt,obj.state_);
            Jc = obj.C_(dt,obj.state_);
            obj.P_ = Jf * obj.P_ * Jf' + Jc * obj.Q_() * Jc';
                       
            if nargout == 1
                varargout{1} = obj.state_;
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