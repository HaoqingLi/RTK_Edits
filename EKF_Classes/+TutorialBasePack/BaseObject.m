classdef BaseObject < handle
    %% base class for all AIS objects
    %  
    % This class creates a unique identifier on its construction. A
    % subclass will always call this superclass constructor during
    % creation. In this way, all child classes can be referenced by this
    % ID, useful Matlab methods like waitfor() can be called with this ID.
    % 
    % The base class itself inherits from handle, which means all internal 
    % changes (e.g. by calling a 'setter'-method of the object) are stored
    % in the object handle. If not derived from handle, you would need to
    % return the complete object to the caller instead!
    % 
    % see also handle
    properties ( GetAccess = 'public', SetAccess = 'private' )
        id
    end
    
    properties ( GetAccess = 'public', SetAccess = 'protected' )
        hFigOff = 100;
    end
    
    properties (Access = protected)
        Done = 0;
    end

    methods ( Access = 'protected' )
        function obj = BaseObject()
            obj.id = TutorialBasePack.BaseObject.increment();
        end
        
    end
    methods ( Static, Access = 'private' )
        function result = increment()
            persistent stamp;
            if isempty( stamp )
                stamp = 0;
            end
            stamp = stamp + uint32(1);
            result = stamp;
        end
    end  
    
    methods
        function ret = getID(obj)
            ret = obj.id;
        end
        
        function fStackName = getNameFromStack(obj,varargin)
        % build current stack name 
            [ST,I]=dbstack;
            fStackName = [];
            if ~isempty(varargin) && strcmpi(varargin{1},'delimiter')
                delim = varargin{2};
            else
                delim = '.';
            end
            sDot={'',delim};
            for sidx=length(ST):-1:min([length(ST), 2])
                fStackName = [fStackName sprintf(['%s' ...
                    sDot{( ((sidx == length(ST)) || (length(ST)==1)) || ...
                    (sidx==1) || (length(ST) > 2 && sidx == 3) ) + 1}],ST(sidx).name)];
            end
        end
    end
end