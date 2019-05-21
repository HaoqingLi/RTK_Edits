classdef ImplBase < AbstractBase
    properties (Access = private)
        mName_
    end

    methods
        function obj = ImplBase()
            obj.mName_ = 'Implementation';
        end

        function setName(obj,name)
            obj.mName_ = name;
        end

        function dVal = getName(obj)
            dVal = obj.mName_;
        end
    end
end