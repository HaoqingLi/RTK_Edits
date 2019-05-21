classdef Child1 < Super
    properties
        mAge_
    end
    methods
        function obj = Child1()
            obj@Super();
            obj.mName_ = 'Child';
            fprintf('I am a %s\n',obj.mName_);
        end
        
        function dVal = getAge(obj)
            dVal = obj.mAge_;
        end
        
        function setAge(obj,age)
            obj.mAge_ = age;
        end
    end
end