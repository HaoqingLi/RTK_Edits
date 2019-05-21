classdef Super2 < handle
    properties (Access = private)
        mSex_
    end

    methods
        function obj = Super2()
            obj.mSex_ = 'male';
            fprintf('I am a %s\n',obj.mSex_);
        end

        function setSex(obj,sex)
            obj.mSex_ = sex;
        end

        function dVal = getSex(obj)
            dVal = obj.mSex_;
        end
    end
end