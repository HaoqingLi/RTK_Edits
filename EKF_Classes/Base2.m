classdef Base2 < handle
    properties (Access = private)
        mName_
    end

    methods
        function obj = Base2()
            obj.mName_ = 'Base';
            fprintf('I am %s\n',obj.mName_);
        end

        function setName(obj,name)
            obj.mName_ = name;
        end

        function dVal = getName(obj)
            dVal = obj.mName_;
        end
    end
end