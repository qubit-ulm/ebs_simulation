classdef PoissonParams
    %POISSONPARAMS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        jump_rate; %array of rates of poisson steps in [1/s]
        delta; %array of stepsizes [nm]
        sigma_delta; %standard deviation of step size in [nm]
    end
    
    methods
        function SPobj = PoissonParams(varargin)
            if (nargin==3)
                SPobj.jump_rate = varargin{1}; 
                SPobj.delta = varargin{2}; 
                SPobj.sigma_delta = varargin{3};
            else
                SPobj.jump_rate = [10, 2]; 
                SPobj.delta = [1, -1]; 
                SPobj.sigma_delta = [0.05, 0.5];
            end
            
        end
        
        function out = getSimpleParams(obj)
            out.lambda = obj.jump_rate;
            out.delta = obj.delta;
            out.sigma_delta = obj.sigma_delta;
        end
        
    end
    
end

