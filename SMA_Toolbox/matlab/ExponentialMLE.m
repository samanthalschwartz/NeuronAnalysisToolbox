classdef ExponentialMLE < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        T; %
        N; %
        data;
    end
    
    methods
        function obj=ExponentialMLE(data)
            if ~isvector(data)
                error('ExponenitalMLE', 'Data must be a vector')
            else
                data=data(:);
            end
            if any(data~=round(data))
                error('ExponenitalMLE', 'Data must be integer valued')
            end
            
        end
    end
    
end

