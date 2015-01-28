classdef SimParams
    %SIMPARAMS class contains all parameters for simulations of SimData
    %   Detailed explanation goes here
    
    properties
        %default parameters
        N; %number of data points
        h; %time increments [s]
        kbT; %thermal energy, room temperature [pN*nm]
        gamma; %in viscosity of trapped bead  [pN*s/nm]
        kappa; %trap stiffness for trap1 and 2 [pN/nm]
        k_DNA; %stiffness of a 1.8kbp long DNA tether (L=600nm)
        L; %contour length of DNA in [nm]
        x0; %initial elongtion from trap
        k; %effective stiffness of relative coordinate x
        instrDrift; %0.1s correlation of instrument fluctuations in [1/s]
        instrNoise; %standard deviation of fluctuations in [nm^2/s]
        mix; %factor how much the beads signal is distorted by instr. Noise
        simtype; %assist or oppose
        constForce; %constant Force = true
        
        p2Pars; %Pol2Params object
        poissonPars; %PoissonParams object        
        
    end
    
    methods
        %constructor
        function SPobj = SimParams(varargin)
            
           if (nargin==14)
                SPobj.N = varargin{1}; 
                SPobj.h = varargin{2};
                SPobj.kbT = varargin{3};
                SPobj.gamma = varargin{4}; 
                SPobj.kappa = varargin{5}; 
                SPobj.k_DNA = varargin{6}; 
                SPobj.L = varargin{7};
                SPobj.x0 = varargin{8};
                SPobj.k = varargin{9};  
                SPobj.instrDrift = varargin{10}; 
                SPobj.instrNoise = varargin{11}; 
                SPobj.mix = varargin{12}; 
                SPobj.simtype = varargin{13};
                SPobj.constForce = varargout{14};
           else
                SPobj.N = 1e6; 
                SPobj.h = 2e-4;
                SPobj.gamma = 1e-5; 
                SPobj.kappa = 0.2; 
                SPobj.kbT = 4.1; 
                SPobj.k_DNA = 0.2; 
                SPobj.k = SPobj.kappa(1)+2*SPobj.k_DNA; 
                SPobj.L = 1020;
                SPobj.x0 = 0.0; 
                SPobj.instrDrift = 1; 
                SPobj.instrNoise = 5; 
                SPobj.mix = 0.05;
                SPobj.simtype = 'assist';
                SPobj.constForce = true;
           end
           
         SPobj.p2Pars = Pol2Params(); 
         SPobj.poissonPars = PoissonParams(); 
        end
        
        function obj = update_k(obj,kdna)
           obj.k = obj.kappa(1)+2*kdna; 
        end
        
        function obj = updateRates(obj)
           obj.p2Pars = obj.p2Pars.updateAllRates(); 
        end
    end
    
end

