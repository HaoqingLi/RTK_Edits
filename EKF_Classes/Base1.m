classdef Base1
    properties (Access = private)
        mName_
    end

    methods
        function obj = Base1()
            obj.mName_ = 'Base';
            fprintf('I am %s\n',obj.mName_);
        end

        function obj = setName(obj,name)
            obj.mName_ = name;
        end

        function dVal = getName(obj)
            dVal = obj.mName_;
        end
    end
end