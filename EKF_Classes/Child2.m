classdef Child2 < Super1 & Super2
    properties (Access = private)
        mAge_
    end
    methods
        function obj = Child2()
            obj@Super1();
            obj@Super2();
            obj.mAge_ = 103;
        end
        
        function dVal = getAge(obj)
            dVal = obj.mAge_;
        end
        
        function setAge(obj,age)
            obj.mAge_ = age;
        end
    end
end