classdef PolyCollector < handle
% Mutable polynomial accumulator for streaming callbacks.
properties
    coefficients
end
methods
    function obj=PolyCollector(D),obj.coefficients=cell(1,D+1);end
    function add(obj,weight,value)
        obj.coefficients{weight+1}=polar_spectrum.BigNat.add(obj.coefficients{weight+1},value);
    end
    function update(obj,poly)
        for k=find(~cellfun(@isempty,poly)),obj.add(k-1,poly{k});end
    end
end
end
