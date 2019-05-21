classdef DataPublisher < TutorialBasePack.BaseObject
    properties
        fname_ % load processed AIS data from fname_
        dtPred_ % time scale for triggering prediction of AIS states
        curVessel_ % data map for vessel to be tracked from AIS messages
        
        vesselID_ % ID of vessel of interest
        
        timeGPS_ % gps time to be followed
        
        dataIdxGlob_ % global AIS msg data index
        
        fhObs_ = []; % mapping function of measurements to needed observation data
        dataObserved_ % matrix of all observed data (subset of AIS msg)
        
        iStartIdx_ = 1 % start index in data map
        iEndIdx_ = [] % end index in data map
        
        %%%%%%% timer settings %%%%%%%%%%%%
        tmrPrediction_
        tmrUpdate_
        bConstantRatePrediction_ = true
                
        idxPModel_ = [] % indexes to match observation data to process model
        idxObsModel_ = [] % indexes to match observation data to measurement model 
        
        elUpObj_ % listener object to restart timer (needs to be deleted to 
                % remove listeners)
        elPredObj_% listener object to restart timer (needs to be deleted to 
                % remove listeners)
                
        iTmrFac_ = 1; % factor to speed up "RT" simulation
        
        dtStartPred % handle to start of time counter for prediction timer (tic)
        dtStartObs % handle to start of time counter for update timer (tic)
        
        tStampLatest_ = [];% time stamp of last activity (timeout of one of the timers)
        %%% handles to conversion functions
        hConvKts2MpS = @(x) x / 1.94384449;
        hConvMpS2Kts = @(x) x * 1.94384449;
        hConvDpM2DpS = @(x) x / 60;
        hConvDpS2DpM = @(x) x * 60;
    end
    
    events
        observation % event in case new observation (measurement was made)
        prediction % event in case new prediction shall be triggered
        restartUpdateTmr % restart correction timer event
        restartPredTmr % restart prediction timer event
    end
    
    methods
        function obj = DataPublisher(varargin)
            
            while ~isempty(varargin)
                if length(varargin) > 1
                    if strcmp(varargin{1},'-fname') && ...
                            ischar(varargin{2})
                        obj.fname_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-vesselID') && ...
                            isnumeric(varargin{2})
                        obj.vesselID_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-dtPrediction') && ...
                            isnumeric(varargin{2})
                        obj.dtPred_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-iFacRT') && ...
                            isnumeric(varargin{2})
                        obj.iTmrFac_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-hObsMapping') && ...
                        isa(varargin{2},'function_handle')
                        obj.fhObs_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-bConstantPredictionRate') && ...
                            islogical(varargin{2})
                        obj.bConstantRatePrediction_ = varargin{2};
                        varargin{2} = [];
                    else
                        myExc = MException('AISDataCaster:Constructor', ...
                        sprintf('unknown proporty "%s"',varargin{1}));
                        throw(myExc);
                    end
                elseif strcmp(varargin{1},'-h') || strcmp(varargin{1},'--help')
                    fprintf('[INFO] help to %s requested:\n\n',obj.getNameFromStack)
                    obj.usage();
                    return
                elseif ~isempty(varargin{1})
                    myExc = MException('AISDataCaster:Constructor', ...
                        sprintf('unknown proporty "%s"',varargin{1}));
                    throw(myExc);
                end
                varargin{1} = [];
                varargin(cellfun(@isempty,varargin)) = [];
            end
           
        end
        
        function loadAISdata(obj,varargin)
            
            while ~isempty(varargin)
                if length(varargin) > 1
                    if strcmp(varargin{1},'-vesselID') && ...
                            isnumeric(varargin{2})
                        obj.vesselID_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-fname') && ...
                            ischar(varargin{2})
                        obj.fname_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-istart') && ...
                            isnumeric(varargin{2})
                        obj.iStartIdx_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-iend') && ...
                            isnumeric(varargin{2})
                        obj.iEndIdx_ = varargin{2};
                        varargin{2} = [];
                    else
                        myExc = MException('AISDataCaster:loadAISdata', ...
                        sprintf('unknown proporty "%s"',varargin{1}));
                        throw(myExc);
                    end
                elseif ~isempty(varargin{1})
                    myExc = MException('AISDataCaster:loadAISdata', ...
                        sprintf('unknown proporty "%s"',varargin{1}));
                    throw(myExc);
                end
                varargin{1} = [];
                varargin(cellfun(@isempty,varargin)) = [];
            end
            
            if isempty(obj.vesselID_)
                myExc = MException('AISDataCaster:loadAISdata', ...
                        'expected particular vessel ID, but none was given');
                throw(myExc);
            end
                        
            tAISdata = load(obj.fname_,'mVesselsMap');
            
            obj.curVessel_ = tAISdata.mVesselsMap(obj.vesselID_);
            
            if isempty(obj.iEndIdx_)
                obj.iEndIdx_ = length(obj.curVessel_.gps);
            end
            
            tmp = obj.curVessel_;
            for field = fieldnames(obj.curVessel_)'
                tmp.(field{1}) = obj.curVessel_.(field{1})(obj.iStartIdx_:obj.iEndIdx_);
            end
            obj.curVessel_ = tmp;
            
            obj.timeGPS_ = obj.curVessel_.gps;
            
            obj.dataIdxGlob_ = 2;
            
            obj.assignObservationData();
        end
        
        function startDataProcessing(obj,varargin)
            while ~isempty(varargin)
                if length(varargin) > 1
                    if strcmp(varargin{1},'-dtPrediction') && ...
                            isnumeric(varargin{2})
                        obj.dtPred_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-iFacRT') && ...
                            isnumeric(varargin{2})
                        obj.iTmrFac_ = varargin{2};
                        varargin{2} = [];
                    elseif strcmp(varargin{1},'-bConstantPredictionRate') && ...
                            islogical(varargin{2})
                        obj.bConstantRatePrediction_ = varargin{2};
                        varargin{2} = [];
                    else
                        myExc = MException('AISDataCaster:startIntegrityCheck', ...
                        sprintf('unknown property "%s"',varargin{1}));
                        throw(myExc);
                    end
                elseif ~isempty(varargin{1})
                    myExc = MException('AISDataCaster:startIntegrityCheck', ...
                        sprintf('unknown property "%s"',varargin{1}));
                    throw(myExc);
                end
                varargin{1} = [];
                varargin(cellfun(@isempty,varargin)) = [];
            end
            
            if isempty(obj.dtPred_)
                myExc = MException('AISDataCaster:startIntegrityCheck', ...
                        'expected time intervall for prediction step');
                throw(myExc);
            end
            
            %% set up timers for correction and prediction steps
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%  update (correction) timer  %%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.tmrUpdate_ = timer();
            obj.tmrUpdate_.StartDelay = (obj.timeGPS_(obj.dataIdxGlob_) - obj.timeGPS_(obj.dataIdxGlob_ - 1)) / obj.iTmrFac_;
            % obj.tmrUpdate_.StartDelay = obj.dtPred_/obj.iTmrFac_;
            obj.tmrUpdate_.TimerFcn = @(~,~)obj.triggerUpdateCB();
            obj.tmrUpdate_.ExecutionMode = 'singleShot';
            obj.tmrUpdate_.StartFcn = @(~,evt)fprintf('[INFO] tmrUpdate_.%s exectued at %s\n',...
                 evt.Type, datestr(evt.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF'));
            obj.tmrUpdate_.StopFcn = @(~,evt)obj.issueUpdateTmrRestartEvent(evt);
            
            % obj.configUpdateTimer();
            
            obj.elUpObj_ = addlistener(obj,'restartUpdateTmr',@(src,evt)obj.eventListener__onUpdateTmrRestart(src,evt));
                    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% prediction timer %%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.tmrPrediction_ = timer();
            obj.tmrPrediction_.StartDelay = obj.dtPred_/obj.iTmrFac_;
            obj.tmrPrediction_.TimerFcn = @(~,~)obj.triggerPredictionCB(); %@(~,~)obj.triggerPredictionCB(obj.dtPred_);
            
            if obj.bConstantRatePrediction_
                obj.tmrPrediction_.ExecutionMode = 'fixedRate';
                obj.tmrPrediction_.Period = obj.dtPred_/obj.iTmrFac_;
                obj.tmrPrediction_.StopFcn = @(~,evt)fprintf('[INFO] tmrPrediction_.%s executed at %s\n',...
                 evt.Type, datestr(evt.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF'));
            else
                obj.tmrPrediction_.ExecutionMode = 'singleShot';
                obj.tmrPrediction_.StopFcn = @(~,evt)obj.issuePredTmrRestartEvent(evt);
                obj.elPredObj_ = addlistener(obj,'restartPredTmr',@(src,evt)obj.eventListener__onPredTmrRestart(src,evt));
            end
            
            obj.tmrPrediction_.StartFcn = @(~,evt)obj.callbackStartPredTimer(evt);       
                        
            obj.tStampLatest_ = 0;
            %% start timers 
            start(obj.tmrPrediction_);
            start(obj.tmrUpdate_);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% HELPER METHODS TO BUILD DATA EVENT %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function assignObservationData(obj)
            %%% convert geodetic to ENU (via ECEF)
            lat = obj.curVessel_.aisLat;
            lon = obj.curVessel_.aisLon;
                       
            %%% replace defaults of lat and lon with interpolated values
            %%% LAT %%%
            idxLATok = find(lat ~= TutorialBasePack.AISDefaultValues.DEF_LAT);
            lat = interp1(idxLATok,lat(idxLATok),(1:length(lat)));
            %%% LON %%%
            idxLONok = find(lon ~= TutorialBasePack.AISDefaultValues.DEF_LON);
            lon = interp1(idxLONok,lon(idxLONok),(1:length(lon)));

            %%% geodetic2ECEF
            a = 6378137.0;
            eSquared = 6.69437999014*1e-3;
            N = @(phi) a ./ sqrt(1 - eSquared*sind(phi).^2);
            h = 0;

            ECEF_X = (N(lat) + h).* cosd(lat).*cosd(lon);
            ECEF_Y = (N(lat) + h).* cosd(lat).*sind(lon);
            ECEF_Z = (N(lat)*(1 - eSquared) + h).* sind(lat);
            
            %%% find reference point close to mean of track
            dTmp = sqrt(ECEF_X.^2 + ECEF_Y.^2);
            
            CDFdist = abs(cumsum(diff(dTmp)))/abs(sum(diff(dTmp)));
            
            idx = find(CDFdist > 0.5,1);
            
            %%% ECEF2ENU
            ENU = [-sind(lon(idx)) cosd(lon(idx)) 0; ...
                -sind(lat(idx))*cosd(lon(idx)) -sind(lat(idx))*sind(lon(idx)) cosd(lat(idx)); ...
                cosd(lat(idx))*cosd(lon(idx)) cosd(lat(idx))*sind(lon(idx)) sind(lat(idx))] * ...
                [ECEF_X - ECEF_X(idx); ECEF_Y - ECEF_Y(idx); ECEF_Z - ECEF_Z(idx)];

            ENU_E = ENU(1,:);
            ENU_N = ENU(2,:);
            ENU_U = ENU(3,:);
            
            %%% in case control input should be used, it is assumed that a
            %%% different state and observation model is applied (bad,
            %%% because hard-coded)
            obj.dataObserved_ = obj.fhObs_(ENU_E, ENU_N, ...
                    obj.curVessel_.aisCOG,...
                    obj.curVessel_.aisSOG, ...
                    obj.curVessel_.aisROT);
        end
        
        function calcUpdateIntervallITU(obj,speedMpS,rotDpS)
            %%% ITU specification 5deg heading change to mean of the last 30s (ITU-R M.1371-4 p48)
            %%% this is encoded in ROT_AIS (see http://www.navcen.uscg.gov/?pageName=AISMessagesA)
            %%%    - ROT_AIS == -127: turning left at more than 5 deg per 30 s
            %%%    - ROT_AIS == +127: turning right at more than 5 deg per 30 s
            %%%
            %%% The update rate is specified in ITU-R M.1371-4 p4, Table 1
            %%%
            %%% input to this method:
            %%%    - speedMpS [m/s]
            %%%    - rotDpS as ROT_AIS [Deg/s]
            %%%
            %%%    speed will be converted to [knots], rot to [Deg/min] and 
            %%%    then to original representation (within [-127 +127])
            
            rotDpM = obj.hConvDpS2DpM(rotDpS);           
            rotOrig = 4.733*sqrt(rotDpM);
            
            isChangingCourse = @(x) abs(x) == 127;
            
            speedKts = obj.hConvMpS2Kts(speedMpS);

            obj.dtPred_ = zeros(1,length(speedKts));
            if (speedKts<14 && ~isChangingCourse(rotOrig))
                obj.dtPred_ = 10;
            elseif ((speedKts<14) && isChangingCourse(rotOrig))
                obj.dtPred_ = 10./3.;
            elseif ((speedKts>=14) && (speedKts<=23) && ~isChangingCourse(rotOrig))
                obj.dtPred_ = 6.;
            elseif ((speedKts>=14) && (speedKts<=23) && isChangingCourse(rotOrig))
                obj.dtPred_ = 2.;
            elseif ((speedKts>23) && ~isChangingCourse(rotOrig))
                obj.dtPred_ = 2.;    
            elseif ((speedKts>23) && isChangingCourse(rotOrig))
                obj.dtPred_ = 2.;    
            end
        
            fprintf('[INFO] update intervall according to ITU-R: %.3fs\n',obj.dtPred_)
            
            obj.dtPred_ = obj.dtPred_ / obj.iTmrFac_;
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CALLBACKS FROM TIMERS %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function triggerPredictionCB(obj)
            
            tElapsed = toc(obj.dtStartPred)
            if obj.bConstantRatePrediction_
                obj.dtStartPred = tic;
            end
            
            fprintf('[INFO] prediction event will be triggered at %s\n',...
                datestr(now,'dd-mmm-yyyy HH:MM:SS.FFF'));
                       
            dTime = tElapsed * obj.iTmrFac_;
            
            obj.tStampLatest_ = obj.tStampLatest_ + dTime;
                      
            tDataObj = TutorialBasePack.AISEventDataPrediction(dTime);
            
            notify(obj,'prediction',tDataObj)                          
        end
        
        function triggerUpdateCB(obj)    
            
            %if (obj.dataIdxGlob_ > 2)
            tElapsedObs = toc(obj.dtStartPred)
            % end
            
            fprintf('[INFO] new measurement arrived ... event will be triggered%s\n',...
                datestr(now,'dd-mmm-yyyy HH:MM:SS.FFF'));
                       
            if isempty(obj.idxObsModel_)
                tOMIdx = (1:size(obj.dataObserved_,1));
            else
                tOMIdx = obj.idxObsModel_;
            end
            observation = obj.dataObserved_(tOMIdx,obj.dataIdxGlob_);
%             
%             if (obj.dataIdxGlob_ > 2)
%                 % dTime = obj.tmrPrediction_.Period * obj.iTmrFac_;
%                 
%                 if tElapsedObs > 1
%                     dTime = floor(tElapsedObs) * obj.iTmrFac_ + rem(tElapsedObs,floor(tElapsedObs));
%                 else
%                     dTime = ceil(tElapsedObs) * obj.iTmrFac_ - (1 - tElapsedObs);
%                 end
%                 dTime = obj.tStampPred_ + dTime;
%             else
%                 dTime = obj.tmrUpdate_.StartDelay;
%             end
            
            dTime = ceil(tElapsedObs) * obj.iTmrFac_ + obj.tStampLatest_
            %% stop prediction timer and execute prediction
            stop(obj.tmrPrediction_);
            obj.triggerPredictionCB();
            
            %% build measurement event data
            tDataObj = TutorialBasePack.AISEventDataUpdate(dTime,observation);
            %%% notify subscribers/listeners
            notify(obj,'observation',tDataObj)
            
            obj.dataIdxGlob_ = obj.dataIdxGlob_ + 1;
            
            %% start prediction timer again
            start(obj.tmrPrediction_);            
            
            if (obj.dataIdxGlob_ > size(obj.dataObserved_,2))
                obj.killTimers();                           
                fprintf('\n[SUCCESS] simluation finished, end of batch measurement was reached!\n')
                return
            end
        end
        
        function issueUpdateTmrRestartEvent(obj,evt)
            fprintf('[INFO] tmrUpdate_.%s exectued at %s\n',...
                 evt.Type, datestr(evt.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF'))
             
             notify(obj,'restartUpdateTmr')
            
        end
        
        function issuePredTmrRestartEvent(obj,evt)
            fprintf('[INFO] tmrPrediction_.%s exectued at %s\n',...
                 evt.Type, datestr(evt.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF'))
             
             notify(obj,'restartPredTmr')
        end
        
        function configUpdateTimer(obj)
            obj.tmrUpdate_ = timer();
            obj.tmrUpdate_.StartDelay = round((obj.timeGPS_(obj.dataIdxGlob_) - obj.timeGPS_(obj.dataIdxGlob_ - 1)) / obj.iTmrFac_ * 1000)/1000;
            % obj.tmrUpdate_.StartDelay = obj.dtPred_;
            obj.tmrUpdate_.TimerFcn = @(~,~)obj.triggerUpdateCB();
            obj.tmrUpdate_.ExecutionMode = 'singleShot';
            obj.tmrUpdate_.StartFcn = @(~,evt)fprintf('[INFO] tmrUpdate_.%s exectued at %s\n',...
                 evt.Type, datestr(evt.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF'));
            obj.tmrUpdate_.StopFcn = @(~,evt)obj.issueUpdateTmrRestartEvent(evt);
        end
        
        function configPredictionTimer(obj)
            
            % calculate current update rate based on ITU specification
            % this sets obj.dtPred_ !
            obj.calcUpdateIntervallITU(obj.curVessel_.aisSOG(obj.dataIdxGlob_),obj.curVessel_.aisROT(obj.dataIdxGlob_));
            
            obj.tmrPrediction_ = timer();
            
            if obj.bConstantRatePrediction_
                obj.tmrPrediction_.StartDelay = obj.dtPred_;
            else
                obj.tmrPrediction_.StartDelay = round(obj.dtPred_ / obj.iTmrFac_*1000)/1000;
            end
            % obj.tmrUpdate_.StartDelay = obj.dtPred_;
            obj.tmrPrediction_.TimerFcn = @(~,~)obj.triggerPredictionCB();
            obj.tmrPrediction_.ExecutionMode = 'singleShot';
            obj.tmrPrediction_.StartFcn = @(~,evt)obj.callbackStartPredTimer(evt);
            obj.tmrPrediction_.StopFcn = @(~,evt)obj.issuePredTmrRestartEvent(evt);
        end
        
        function callbackStartPredTimer(obj,evt)
            obj.dtStartPred = tic;
            
            fprintf('[INFO] tmrPrediction_.%s executed at %s\n',...
                 evt.Type, datestr(evt.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF'))
        end
        
        function killTimers(obj)
            
            delete(obj.elUpObj_);
            
            if isobject(obj.elPredObj_)
                delete(obj.elPredObj_);
            end
            
            stop(obj.tmrUpdate_);
            stop(obj.tmrPrediction_);
            
            delete(obj.tmrPrediction_)
            delete(obj.tmrUpdate_)
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% GETTER/SETTER METHODS %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        function setObsIdxs(obj,idxObsModel)
            obj.idxObsModel_ = idxObsModel;
        end
        
        function dInit = getInitialState(obj,idxProcModel)

            if nargin == 1 && ~isempty(obj.idxPModel_)
                idxProcModel = obj.idxPModel_;
            else
                obj.idxPModel_ = idxProcModel;
            end

            dInit = obj.dataObserved_(idxProcModel,1);
        end
    end
    
    methods (Access=private)
        function eventListener__onUpdateTmrRestart(obj,~,~)
            
            delete(obj.tmrUpdate_)
            
            obj.configUpdateTimer();
            
            % obj.dtStartObs = tic;
            
            start(obj.tmrUpdate_);
        end
        
        function eventListener__onPredTmrRestart(obj,~,~)
            
            delete(obj.tmrPrediction_)
            
            obj.configPredictionTimer();
            
            obj.dtStartPred = tic;
            
            start(obj.tmrPrediction_);
        end
        
    end
    
    
    methods (Static)
        function killTimersAll()
            stop(timerfindall);
            delete(timerfindall);
        end
    end
    
end