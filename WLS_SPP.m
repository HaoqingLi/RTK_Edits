% -------------------------------------------------------------------------
% WLS Class for solving the Single Point Positioning (SPP), snapshot
% solution
%
% Author: 
%   Daniel Arias Medina
%   German Aerospace Center (DLR) 
%   Institute of Communication and Navigation, Nautical Systems Department
%   Neustrelitz, Germany
% Version: 14.06.16
%-------------------------------------------------------------------------
classdef WLS_SPP < handle
    properties (SetAccess = public, GetAccess = public)
        
        Position                    = [0 0 0 ]';                            % Estimated position
        CLKoffset                   = 0;                                    % Estimated clock offset
        Velocity                    = [0 0 0]';                             % Estimated velocity
        CLKoffset_rate              = 0;                                    % Estimated clock offset rate
        
        P_position                    = eye(4);                               % P covariance matrix for the estimated state: position and CLK offset
        P_velocity_matrix                    = eye(4);                               % P covariance matrix for the estimated state: position and CLK offset
        
        Sigma_PS    = [];                                   % Variance Covariance matrix of the observated pseudoranges
        Sigma_PSrate         = [];                                   % Variance Covariance matrix of the observated pseudorange rates    
        
        max_iterations              = 20;                                   % Maximum number of iterations allowed for the Least Squares method
        error_tolerance             = 10^(-3);                              % Condition to stop the iteration
        
        DOP                         = zeros(5,1);
        
        noise_model_parameter       = [0.60006, 50.6392, 33.8385];
        
        noise_model                 = 1;
        
        iteration                   = [];
        
        LS_residuals                = [];
        
        G_matrix                    = [];
        
        HPL                         = [];
        VPL                         = [];
               
        nav_requirements            = struct('mode',1,'continuity',0.997,'ContTimeInterval',60*15,'HorizontalPosAccuracy','10',...
            'VerticalPosAccuracy',15,'alertLimit',25);
        
        debugMode                   = 0;
         
        
    end
    
    methods    (Access = public)
        %% Constructor of the Class
        function obj = WLS_SPP( in_Position, in_Sigma_PS, in_Sigma_PSrate, in_max_iterations, in_error_tolerance, in_noise_model, in_noise_model_parameter, in_navigationMode )
            obj.Position                            = in_Position;
            obj.Velocity                            = [0,0,0]';
            
            obj.P_position                          = eye(4);
            obj.P_velocity_matrix                   = eye(4);
            
            obj.Sigma_PS                            = in_Sigma_PS;
            obj.Sigma_PSrate                        = in_Sigma_PSrate;
            
            obj.max_iterations                      = in_max_iterations;
            obj.error_tolerance                     = in_error_tolerance;
            
            obj.DOP                                 = zeros(5,1);
            
            obj.noise_model_parameter               = in_noise_model_parameter;
            
            obj.noise_model                         = in_noise_model;
            
            obj.nav_requirements.mode               = in_navigationMode;
            switch obj.nav_requirements.mode
                case 1 % Navigation Mode - Coastal - default navigation requirements (per operation basis) converted to algorithm requirements
                case 2 % Navigation Mode - port approach and restricted waters
                    obj.nav_requirements.continuity = 0.9997;
                    obj.nav_requirements.ContTimeInterval = 60*15; %Resolution A.1046(27)
                    obj.nav_requirements.HorizontalPosAccuracy = 10;
                case 3 % Navigation Mode - Port
                    obj.nav_requirements.continuity = 0.9997;
                    obj.nav_requirements.ContTimeInterval = 60*15; %Resolution A.1046(27)
                    obj.nav_requirements.HorizontalPosAccuracy = 1;
                    obj.nav_requirements.alertLimit = 2.5;
                case 4 % Navigation Mode - Inland Waters (IRIS Europe 2 project - category B operation)
                    obj.nav_requirements.continuity = 0.9997;
                    obj.nav_requirements.ContTimeInterval = 60*15; %Resolution A.1046(27)
                    obj.nav_requirements.HorizontalPosAccuracy = 10;
                otherwise
            end
        end
        
        
        %% Defining set methods for the class access: set group
        function obj = set_EstimatedPosition(obj, in_EstimatedPosition)     % Set the estimated position
            obj.Position = in_EstimatedPosition;
        end
        
        function obj = set_EstimatedVelocity(obj, in_v)                     % Set the estimated velocity
            obj.Position = in_v;
        end
        
        function obj = set_max_iterations (obj, in_max_iterations)
            obj.max_iterations = in_max_iterations;
        end
        
        function obj = set_error_tolerance (obj, in_error_tolerance)
            obj.error_tolerance = in_error_tolerance;
        end
        
        function obj = set_model_noise_parameter (obj, in_model_noise_parameter)
           obj.noise_model_parameter = in_model_noise_parameter;
        end
        
        
        
        %% Defining set methods for the class access: get group
        function out_EstimatedPosition = get_Position(obj)            % Get the estimated position
            out_EstimatedPosition = obj.Position(1:3);
        end
        
        function out_EstimatedClockOffset = get_CLKoffset(obj)            % Get the estimated position
            out_EstimatedClockOffset = obj.CLKoffset;
        end
        
        function out_EstimatedVelocity = get_Velocity(obj)            % Get the estimated velocity
            out_EstimatedVelocity = obj.Velocity;
        end
        
        function out_EstimatedClockOffsetRate = get_CLKoffset_rate(obj)            % Get the estimated position
            out_EstimatedClockOffsetRate = obj.CLKoffset_rate;
        end
        
        
        
        %% Snapshot for position estimation
        function calculate_position(obj, pseudorange, ionosphere, troposphere, svoffset, posX, posY, posZ, CN0, ElevationAng)
            
            % Eliminate from the vector those satellites that do not have info
            ind_2eliminate = find(isnan(pseudorange));
            ind_2use       = find(~isnan(pseudorange));
            
            obj.LS_residuals  = NaN(size(pseudorange));
            
            pseudorange(ind_2eliminate)=[];
            ionosphere(ind_2eliminate)=[];
            troposphere(ind_2eliminate)=[];
            svoffset(ind_2eliminate)=[];
            posX(ind_2eliminate)=[];
            posY(ind_2eliminate)=[];
            posZ(ind_2eliminate)=[];
            CN0(ind_2eliminate)=[];
            ElevationAng( ind_2eliminate ) = [];
            
            % Number of satellites availables
            LOS_satellites = length(pseudorange);
            if LOS_satellites < 4
                if obj.debugMode, disp('WLS Pos: Not enough satellites!'); end
                return;
            end
            
            % Adaptive pseudorange measurement noise covariance model Sigma
            %  Sigma^2 = a + b*10^((-CN0-c)/10)
            % where the approximation parameters are properties of the class and come from [On PNT Integrity in Snapshot and Recursive Positioning Algorithms for Maritime Applications].
            % By default, a = 0.60006, b = 50.6392, c = 33.8385
            switch obj.noise_model
                case 0
                    adaptive_noise = ones(size(pseudorange))*obj.Sigma_PS;
                case 1%obj.noise_model == 1
                    adaptive_noise = obj.noise_model_parameter(1) + obj.noise_model_parameter(2)*10.^(-(CN0-obj.noise_model_parameter(3))/10);
                case 2%obj.noise_model == 2
                    adaptive_noise = 1./sin(ElevationAng).^2;
                case 3%obj.noise_model == 3
                    adaptive_noise = ( 0.13 + 0.56*exp( -ElevationAng/0.1745 ) ).^2;
                otherwise
                    adaptive_noise = ones(size(pseudorange));
            end
            R = diag(adaptive_noise);
            
            % Initialization of the variables
            d      = zeros(LOS_satellites,1);
            G          = zeros(LOS_satellites,4);
            y     = zeros(LOS_satellites,1);
            deltaSolution          = [0,0,0,0]' + obj.error_tolerance + 1;   % Initialize so that it can pass the first control
            Position = [0,0,0]';
            CLKoffset = 0;
            LS_iteration        = 0;
            W = inv(R);
            
            % Iterative LS method to estimate position
            while LS_iteration < obj.max_iterations && norm(deltaSolution) > obj.error_tolerance 
                LS_iteration        = LS_iteration +1;
                for iSV = 1:LOS_satellites
                    d(iSV) = sqrt( (Position(1)-posX(iSV))^2 + (Position(2)-posY(iSV))^2 + (Position(3)-posZ(iSV))^2 );
                    y(iSV) = pseudorange(iSV) - ionosphere(iSV) - troposphere(iSV) + svoffset(iSV) - d(iSV) - CLKoffset;
                    G(iSV,:) = [ (Position(1)-posX(iSV))/d(iSV), (Position(2)-posY(iSV))/d(iSV), (Position(3)-posZ(iSV))/d(iSV), 1 ];
                end
                
                if size(W,1) == size(G,1) && size(W,1) == length(y)
                    if cond(G) < 10^3 && rank(W)>=4
                        %Least square method
                        deltaSolution          = (G' * W * G) \ G' * W * y;
                        % Gauss Newton Method to update the position and clockoffset solution
                        Position            = Position + deltaSolution(1:3);
                        CLKoffset           = CLKoffset + deltaSolution(4);
                    else
                        if obj.debugMode, disp('WLS Pos: Singular matrix'); end
                        obj.Position = nan(3,1);
                        obj.CLKoffset = NaN;
                        return;
                    end
                else
                    if obj.debugMode, disp('WLS Pos: Mismatch on the input data'); end
                    obj.Position = nan(3,1);
                    obj.CLKoffset = NaN;
                    return;
                end                
            end
            
            obj.Position = Position;
            obj.CLKoffset = CLKoffset;
            obj.P_position = inv(G'*W*G) * G' * W * R * W * G * inv(G'*W*G);
           
        end

        
        %% Snapshot for velocity calculation 
        
        function calculate_velocity(obj, posX, posY, posZ, SV_velX, SV_velY, SV_velZ, SV_doppler,CN0)
            %*********************************************************************************
            % Classical SPP for the velocity calculation using the method from Luis presentation
            %*********************************************************************************
            
            % Eliminate from the vector those satellites that do not have info
            ind_2eliminate = find(isnan(posX));
            posX(ind_2eliminate)=[];
            posY(ind_2eliminate)=[];
            posZ(ind_2eliminate)=[];
            SV_velX(ind_2eliminate)=[];
            SV_velY(ind_2eliminate)=[];
            SV_velZ(ind_2eliminate)=[];
            SV_doppler(ind_2eliminate)=[];
            CN0(ind_2eliminate)=[];
            
            % Number of satellites availables
            LOS_satellites = length(posX);
            if LOS_satellites < 4
                if obj.debugMode, disp('WLS Vel: Not enough satellites!'); end
                return;
            end
            
            % Check availability of the receiver position
            if isnan(obj.Position)
               disp('WLS Vel: User position not found') 
               obj.Velocity = nan(3,1);
               obj.CLKoffset_rate = NaN;
               return;
            end
            
            % Adaptive pseudorange measurement noise covariance model Sigma
            %  Sigma^2 = a + b*10^((-CN0-c)/10)
            % where the approximation parameters are properties of the
            % class and come from [On PNT Integrity in Snapshot and
            % Recursive Positioning Algorithms for Maritime Applications].
            % By default, a = 0.60006, b = 50.6392, c = 33.8385
            switch obj.noise_model
                case 0
                    adaptive_noise = ones(size(pseudorange))*obj.Sigma_PSrate;
                case obj.noise_model == 1
                    adaptive_noise = obj.noise_model_parameter(1) + obj.noise_model_parameter(2)*10.^(-(CN0-obj.noise_model_parameter(3))/10);
                case obj.noise_model == 2
                    adaptive_noise = obj.Sigma_PSrate./sin(ElevationAng).^2;
                case obj.noise_model == 3
                    adaptive_noise = ( 0.13 + 0.56*exp( -ElevationAng/0.1745 ) ).^2;
                otherwise
                    adaptive_noise = ones(size(SV_doppler))*obj.Sigma_PSrate;
            end
            R = diag(adaptive_noise);            
            
            % Inititiaze rho_0 with the initial position information
            GeoDistance_SV2REC      = zeros(LOS_satellites,1);
            G_Matrix                = zeros(LOS_satellites,4);
            pseudorange_rate        = zeros(LOS_satellites,1);
            c_LightSpeed            = 299792458;
            L1_freq                 = 1575.42e6;
            T_matrix                = zeros(LOS_satellites,1);
            
            
            for iSatellite = 1:LOS_satellites
                % Distance from the satellites to the iteratively-changing target position
                GeoDistance_SV2REC(iSatellite) = sqrt( (obj.Position(1)-posX(iSatellite))^2 + (obj.Position(2)-posY(iSatellite))^2 + (obj.Position(3)-posZ(iSatellite))^2 );
                
                % Construction of the G measurement matrix
                G_Matrix(iSatellite,:) = [ (obj.Position(1)-posX(iSatellite))/GeoDistance_SV2REC(iSatellite),    (obj.Position(2)-posY(iSatellite))/GeoDistance_SV2REC(iSatellite),...
                    (obj.Position(3)-posZ(iSatellite))/GeoDistance_SV2REC(iSatellite),                           1                               ];
                
                % Pseudorange rate
                pseudorange_rate(iSatellite) = [c_LightSpeed*SV_doppler(iSatellite)/L1_freq];
                
                % T Matrix (Luis presentation)
                T_matrix(iSatellite) = pseudorange_rate(iSatellite) + G_Matrix(iSatellite,:)*[SV_velX(iSatellite); SV_velY(iSatellite); SV_velZ(iSatellite); 0];
            end
            
            %Least square method
            if size(R,1) == size(G_Matrix,1) && size(R,1) == length(T_matrix)
                if cond(G_Matrix) < 10^3 && rank(R)>=4
                    Velocity_offsetrate = (G_Matrix' * inv(R) * G_Matrix) \ G_Matrix' * inv(R) * T_matrix;
                else
                   if obj.debugMode, disp('WLS Vel: Singular matrix'); end
                   obj.Velocity = nan(3,1);
                   obj.CLKoffset_rate = NaN;
                   return;
                end
            else
               if obj.debugMode, disp('WLS Vel: Mismatch on the input data'); end 
               obj.Velocity = nan(3,1);
               obj.CLKoffset_rate = NaN;
               return;
            end
                
            
            % Update Velocity & Offset change rate
            obj.Velocity  = Velocity_offsetrate(1:3);
            obj.CLKoffset_rate        = Velocity_offsetrate(4);
            obj.P_velocity_matrix = inv(G_Matrix' * inv(R) * G_Matrix);
            
        end
        
        %% Snapshot for velocity calculation (Kelly Approach)  
        function calculate_velocity_Kelly(obj, posX, posY, posZ, SV_velX, SV_velY, SV_velZ, SV_doppler)
            %**********************************************************************
            % Velocity SPP solver using Kelly method (Romney paper, VELOCITY RAIM, only velocity algorithm implementation)
            %**********************************************************************
            
            % Eliminate from the vector those satellites that do not have info
            posX(isnan(posX))=[];
            posY(isnan(posY))=[];
            posZ(isnan(posZ))=[];
            SV_velX(isnan(SV_velX))=[];
            SV_velY(isnan(SV_velY))=[];
            SV_velZ(isnan(SV_velZ))=[];
            SV_doppler(isnan(SV_doppler))=[];
            
            % Number of satellites availables
            LOS_satellites = length(SV_velX);
            if LOS_satellites < 4
                return;
            end
            
            % Inititiaze rho_0 with the initial position information
            c_LightSpeed            = 299792458;                            % Speed of Light
            L1_freq                 = 1575.42e6;                            % Frequency of GPS signals
            a_vector_SV2Rec         = zeros(LOS_satellites,3);              % Unit vector from receiver to satellite positions
            dummyRange_SV2Rec       = zeros(LOS_satellites,1);              % Norm of the vector from receiver to satellite position
            SV_velocity             = [SV_velX; SV_velY; SV_velZ]';          % Satellites' velocity
            SV_position             = [posX; posY; posZ]';                   % Satellites' position
            SV_PSrate               = c_LightSpeed*SV_doppler/L1_freq;      % Satellites' pseudorange rate
            G_Matrix                = ones(LOS_satellites, 4);              % G geometric matrix
            R_PS                    = diag(ones(LOS_satellites,1)*obj.Sigma_PS);% Covariance matrix of the pseudoranges error noise
            R_PSrate                = diag(ones(LOS_satellites,1)*obj.Sigma_PSrate); % Covariance matrix of the pseudoranges rate error noise
            B_matrix                = zeros(LOS_satellites,3);              % B matrix --> look at [ Velocity RAIM (Bruce)] paper
            Delta_d_matrix          = zeros(LOS_satellites,1);              % Difference between the SV velocity and the velocity estimated using the pseudorange rate measurement
            
            
            for iSatellite = 1:LOS_satellites
                dummyRange_SV2Rec(iSatellite) = sqrt( (obj.Position(1)-posX(iSatellite))^2 + (obj.Position(2)-posY(iSatellite))^2 + (obj.Position(3)-posZ(iSatellite))^2 );
                a_vector_SV2Rec(iSatellite,:) = (SV_position(iSatellite,:) - obj.Position')/dummyRange_SV2Rec(iSatellite);
                G_Matrix(iSatellite,1:3)      = a_vector_SV2Rec(iSatellite,:);
                B_matrix(iSatellite,:)        = SV_velocity(iSatellite,:)/dummyRange_SV2Rec(iSatellite)- [SV_position(iSatellite,:) - obj.Position']*(dot(SV_velocity(iSatellite,:),a_vector_SV2Rec(iSatellite,:)))/dummyRange_SV2Rec(iSatellite)^2;
                Delta_d_matrix(iSatellite)  = - SV_PSrate(iSatellite) + dot(SV_velocity(iSatellite,:),a_vector_SV2Rec(iSatellite,:));
            end
            
            % Covariance of the receiver's position
            J_3          = zeros(3,4);
            J_3(1:3,1:3) = eye(3);
            R_pos       = J_3*obj.P_position*J_3';
            % Covariance of the velocity residuals
            R_DeltaD    = B_matrix*R_pos*B_matrix' + R_PSrate;
%             obj.P_velocity_matrix = R_DeltaD;
            % Weighted pseudo-inverse of the measurement matrix
            H_aug       = (G_Matrix'*inv(R_DeltaD)*G_Matrix)\G_Matrix'*inv(R_DeltaD);
            % Least Squares method
            Estimated_VelocityKelly = H_aug*Delta_d_matrix;
            obj.Velocity = Estimated_VelocityKelly(1:3) ;
            obj.CLKoffset_rate   = - Estimated_VelocityKelly(4);
            
            obj.P_velocity_matrix = inv(H_aug *inv(R_DeltaD) * H_aug');
%             obj.P_velocity_matrix = inv(H_aug *inv(R_PSrate) * H_aug');
        
        
    end
    end
end


