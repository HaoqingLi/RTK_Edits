classdef Super1 < handle
    properties (Access = private)
        mName_
    end

    methods
        function obj = Super1()
            obj.mName_ = 'Parent';
            fprintf('I am a %s\n',obj.mName_);
        end

        function setName(obj,name)
            obj.mName_ = name;
        end

        function dVal = getName(obj)
            dVal = obj.mName_;
        end
    end
end