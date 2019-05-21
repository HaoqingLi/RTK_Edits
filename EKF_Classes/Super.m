classdef Super < handle
    properties
        mName_
    end

    methods
        function obj = Super()
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
    
    methods (Access=private)
        function test(obj)
            fprintf('I am private')
        end
    end
end