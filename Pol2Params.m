classdef Pol2Params
    %POL2PARAMS contains all relevant parameters for Pol2 simulation.
    % The numerical values of the parameters are published in:
    % Dangkulwanich et. al.,  DOI: 10.7554/eLife.00971
    
    properties
        
        kbT;
        delta;
        F;
        k0;
        k_f;
        k_bn;
        k_b1;
        kb0
        K_D;
        c_NTP;
        k_1;
        k_rev1;
        k_10;
        k_rev10;
        k_3;
        k2_net; % approximates NTP binding + catalsis as a single step
        k1_net
        k1; %net forward elongation
        
        
   
    end
    
    methods
        %constructor
        function PolObj = Pol2Params(varargin)
            
            PolObj.kbT = 4.11; %pN*nm (298K)
            PolObj.delta = 0.17; %nm 0.5bp distance to intermediate state
            PolObj.F = 6.5; %6.5 pN external assisting force
            PolObj.k0 = 1.3; %1/s
            PolObj.kb0 = 6.9;
            PolObj.k_f = PolObj.calcForwardRate(PolObj.k0,PolObj.F);
            PolObj.k_bn = PolObj.calcBackwardRate(PolObj.k0,PolObj.F);
            PolObj.k_b1 = PolObj.calcBackwardRate(PolObj.kb0,PolObj.F);
            PolObj.K_D = 9.2; %NTP dissociation constant \mu M
            PolObj.c_NTP = 35; %NTP concentration in \mu M;
            PolObj.k_10 = 88; %1/s
            PolObj.k_rev10 = 680; %1/s
            PolObj.k_1 = PolObj.calcForwardRate(PolObj.k_10,PolObj.F);
            PolObj.k_rev1 = PolObj.calcBackwardRate(PolObj.k_rev10,PolObj.F);
            PolObj.k_3 = 35; %1/s 
            PolObj.k2_net = PolObj.c_NTP*PolObj.k_3/PolObj.K_D;
            PolObj.k1_net = PolObj.k_1*PolObj.k2_net/(PolObj.k_rev1+PolObj.k2_net);
            %elongation rate
            PolObj.k1 = 1/(1/PolObj.k1_net+1/PolObj.k2_net+1/PolObj.k_3);
            
            
            
            if(nargin == 14)
                PolObj.kbT = varargin{1}; %pN*nm (298K)
                PolObj.delta = varargin{2}; 
                PolObj.F  = varargin{3};
                PolObj.k0  = varargin{4};
                PolObj.k_f = varargin{5};
                PolObj.k_bn  = varargin{6};
                PolObj.k_b1 = varargin{7};
                PolObj.K_D  = varargin{8};
                PolObj.c_NTP = varargin{9};
                PolObj.k_1 = varargin{10};
                PolObj.k_rev1 = varargin{11};
                PolObj.k_3 = varargin{12};
                PolObj.k2_net = varargin{13};
                PolObj.simtype = varargin{14};
            end
            
        end
        
        function [k] = calcForwardRate(obj,k0,F)
            k = k0.*exp(F.*obj.delta./(obj.kbT));
        end
        
        function [k] = calcBackwardRate(obj,k0,F)
            k = k0.*exp(-F.*(0.34-obj.delta)./(obj.kbT));
        end
        
        % update rates:
        function [obj] = updateDiffusionRates(obj,varargin)
            if nargin>1
               Force = varargin{1};
               obj.F = Force;
            else
               Force = obj.F; 
            end
            obj.k_f = obj.calcForwardRate(obj.k0,Force);
            obj.k_bn = obj.calcBackwardRate(obj.k0,Force);
            obj.k_b1 = obj.calcBackwardRate(obj.kb0,Force);
            obj.k_1 = obj.calcForwardRate(obj.k_rev10,Force);
            obj.k_rev1 = obj.calcBackwardRate(obj.k_rev10,Force);
        end
        
        function obj = updateElongRates(obj,varargin)
            if nargin>1
                cntp = varargin{1};
                obj.c_NTP = cntp;
                if nargin>2
                    obj.F = varargin{2};
                end
            end
            obj.k_1 = obj.calcForwardRate(obj.k_10,obj.F);
            obj.k_rev1 = obj.calcBackwardRate(obj.k_rev10,obj.F);
            obj.k2_net = obj.c_NTP*obj.k_3/obj.K_D;
            obj.k1_net = obj.k_1*obj.k2_net/(obj.k_rev1+obj.k2_net); 
            obj.k1 = 1/(1/obj.k1_net+1/obj.k2_net+1/obj.k_3);
        end
        
        function obj = updateAllRates(obj,varargin)
            if nargin>1
                cntp = varargin{1};
                obj.c_NTP = cntp;
                if nargin>2
                    obj.F = varargin{2};
                end
            end
            obj.k_f = obj.calcForwardRate(obj.k0,obj.F);
            obj.k_bn = obj.calcBackwardRate(obj.k0,obj.F);
            obj.k_b1 = obj.calcBackwardRate(obj.kb0,obj.F);

            obj.k_1 = obj.calcForwardRate(obj.k_10,obj.F);
            obj.k_rev1 = obj.calcBackwardRate(obj.k_rev10,obj.F);
            obj.k2_net = obj.c_NTP*obj.k_3/obj.K_D;
            obj.k1_net = obj.k_1*obj.k2_net/(obj.k_rev1+obj.k2_net); 
            obj.k1 = 1/(1/obj.k1_net+1/obj.k2_net+1/obj.k_3);
        end
        
        function [k1,k_b,k_f,k_b1] = returnSimRates(obj)
            k1 = obj.k1;
            k_b = obj.k_bn;
            k_f = obj.k_f;
            k_b1 = obj.k_b1;
        end
        
        
    end
    
end

