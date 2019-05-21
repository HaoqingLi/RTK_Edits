classdef AbstractKalmanFilter < TutorialBasePack.BaseObject
    properties
        %%%%%%%%%%%%%%%%%%%
        %%%   GENERAL   %%%
        %%%%%%%%%%%%%%%%%%%
        status_ = TutorialBasePack.BaseStates.Idle
       
        %%%%%%%%%%%%%%%%%%%
        %%% KF ELEMENTS %%%
        %%%%%%%%%%%%%%%%%%%
        K_ % Kalman gain
        f_ % state transition 
        S_ % residual covariance matrix
        B % control matrix: expected impact of input acceleration
        u % control input
        h_ % measurement transition        
        state_ % state vector
        Q_ % covariance matrix of measurement noise
        R_ % covariance matrix of process noise
        P_ % covariance matrix of state error
        y_ % innovation (residual) vector

        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% OBSERVER PATTERN %%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        masterObj_ = [] % master issuing events (not required, empty on default)
        
        plot_ = true    % flag to enable/disable intermediate plotting
        hPlotPred_ = [] % array of figure handles for plotting intermediate predicted state values
        hPlotUp_ = [] % array of figure handles for plotting intermediate corrected state values
        cTitleLabels_ = [] % cells of titles for plots
        cYlabelLabels_ = [] % cells of ylabels for plots
    end
    
    methods 
        function obj = AbstractKalmanFilter(varargin)
            
            while ~isempty(varargin)
                if length(varargin) > 1
                    if strcmp(varargin{1},'-masterObj') && ...
                            isobject(varargin{2})
                        obj.masterObj_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-plot') && ...
                            islogical(varargin{2})
                        obj.plot_ = varargin{2};
                        varargin{2} = [];
                    else
                        myExc = MException('AbstractKalmanFilter:Constructor', ...
                        sprintf('unknown proporty %s',varargin{1}));
                        throw(myExc);
                    end
                elseif strcmp(varargin{1},'-h') || strcmp(varargin{1},'--help')
                    fprintf('[INFO] help to %s requested:\n\n',obj.getNameFromStack)
                    obj.usage();
                    return
                elseif ~isempty(varargin{1})
                    myExc = MException('AbstractKalmanFilter:Constructor', ...
                        sprintf('unknown proporty %s',varargin{1}));
                    throw(myExc);
                end
                varargin{1} = [];
                varargin(cellfun(@isempty,varargin)) = [];
            end
            
        end
        
        function usage(obj)
            disp('******************************************************');
            eval(['help ' class(obj)])
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ABSTRACT METHODS - DECLARATION %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Abstract)
        
        init(obj)
        
        prediction(obj,dt)
        
        dRet = correction(obj,sObs)
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ABSTRACT LISTENERS - DECLARATION %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Abstract)
        eventListener__onObservation(obj,src,evt)
        
        eventListener__onPrediction(obj,src,evt)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% GETTER METHODS (not abstract) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods 
        function dRet = getPDiagElements(obj)
            dRet = diag(obj.P_);
        end
        
        function dRet = getSCov(obj)
            dRet = obj.S_;
        end
        
        function dRet = getKGain(obj)
            dRet = diag(obj.K_);
        end
        
        function dRet = getInnovation(obj)
            dRet = obj.y_;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% HELPER METHODS (not abstract) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods 
        function plotIntermediates(obj,type,dt)
            
            persistent xdataUp;
            persistent xdataPred;
            persistent xdataFailureIdx;
            persistent ydataFailureIdx;
                               
%             if ~isempty(obj.curFailureObs_) && ~isempty(obj.curFailurePred_)
%                 obj.curFailureIdx_ = obj.curFailureObs_ | obj.curFailurePred_;
%             end
            
%             if isempty(obj.hPlotFail_) && ~isempty(obj.curFailureIdx_)
%                 figure(obj.hFigOff + length(obj.state_))
%                 hold on
%                 xdataFailureIdx = dt;
%                 ydataFailureIdx = double(obj.curFailureIdx_);
%                 obj.hPlotFail_ = plot(dt,double(obj.curFailureIdx_));
%                 % set(obj.hPlotFail_,'XData',xdataFailureIdx,'YData',ydataFailureIdx);
%                 title('integrity measure (1 == Failure detected)')
%             end
            
            switch type
                case 'prediction'
                    if isempty(obj.hPlotPred_)
                        xdataPred = dt;
                        for n=1:length(obj.state_)
                            eval(['persistent ydataPred' num2str(n) ';']);
                            figure(obj.hFigOff + n-1)
                            hold on
                            
                            obj.hPlotPred_(n) = plot(dt,obj.state_(n));
                            % get(obj.hPlotPred_(n))
                            eval(['ydataPred' num2str(n) ' = obj.state_(n);']);
                            
                            if length(obj.cTitleLabels_) >= n && iscell(obj.cTitleLabels_)
                                title(obj.cTitleLabels_{n})
                            end
                            if length(obj.cYlabelLabels_) >= n && iscell(obj.cYlabelLabels_)
                                ylabel(obj.cYlabelLabels_{n})
                            end
                            xlabel('time [s]')
                            
                            % set(obj.hPlotPred_(n),'XData',xdataPred);
                            % set(obj.hPlotPred_(n),'YData',eval(['ydataPred' num2str(n)]));                            
                        end
                        
                        %%% don't continue from here
                        return 
                    end
                case 'update'
                    if isempty(obj.hPlotUp_)
                        xdataUp = dt;
                        for n=1:length(obj.state_)
                            eval(['persistent ydataUp' num2str(n) ';']);
                            figure(obj.hFigOff + n-1)
                            hold on
                            obj.hPlotUp_(n) = plot(dt,obj.state_(n));
                            eval(['ydataUp' num2str(n) ' = obj.state_(n);']);
                            % set(obj.hPlotUp_(n),'XData',xdataUp);
                            % set(obj.hPlotUp_(n),'YData',eval(['ydataUp' num2str(n)]));
                        end    
                        
                        %%% don't continue from here
                        return 
                    end
                    
            end
                        
            switch type
                case 'update'
                    xdataUp = [xdataUp dt]; % dt contains the new absolute time stamp already!!
                case 'prediction'
                    xdataPred = [xdataPred xdataPred(end)+dt]; % dt is the relative time difference from previous prediction
            end
            
            for n=1:length(obj.state_)
                figure(obj.hFigOff + n-1)
                switch type
                    case 'update'
                        eval(['persistent ydataUp' num2str(n) ';']);
                        eval(['ydataUp' num2str(n) ' = [ydataUp' num2str(n) ' obj.state_(n)];']);
                        set(obj.hPlotUp_(n),'XData',xdataUp,'YData',eval(['ydataUp' num2str(n)]));
                    case 'prediction'
                        eval(['persistent ydataPred' num2str(n) ';']);
                        eval(['ydataPred' num2str(n) ' = [ydataPred' num2str(n) ' obj.state_(n)];']);
                        set(obj.hPlotPred_(n),'XData',xdataPred,'YData',eval(['ydataPred' num2str(n)]));
                end
                refreshdata(obj.hFigOff + n-1)
            end  
            
%             if ~isempty(obj.curFailureIdx_)
%                 figure(obj.hFigOff + n)
%                 xdataFailureIdx = [xdataFailureIdx xdataFailureIdx(end)+dt];
%                 ydataFailureIdx = double([ydataFailureIdx obj.curFailureIdx_]);
%                 set(obj.hPlotFail_,'XData',xdataFailureIdx,'YData',ydataFailureIdx);
%                 refreshdata(obj.hFigOff + n)
%                 title('integrity measure (1 == Failure detected)')
%             end
            
            fprintf('[INFO] xdataUp contains %d\n',length(xdataUp))
            fprintf('[INFO] xdataPred contains %d values\n',length(xdataPred))
        end
    end
    
end