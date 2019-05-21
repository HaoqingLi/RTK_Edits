classdef AbstractBase < handle
    methods (Abstract)
        setName(obj,name)
        
        val = getName(obj)
    end
end