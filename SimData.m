classdef SimData
    %SimData class contains noisy timeseries, the original piecewise 
    %constant signal and various methods to simulate such signals.
    %
    % call:
    %   sd = SimData() creates timeseries from default SimParams property
    %   sd = SimData(SimParams sp) creates timeseries with specified parameters
    %
    % example:
    %   create a poissonian noisy timeseries by:
    %   sp = SimParams();
    %   sd = SimData(sp);
    %   sd = sd.simulatePoissonJumps(); fills data and pwcs with data
    
    properties 
        time; % time axis
        data; % noisy data
        pwcs; % a priori smulated steps without noise
        p; %time points of jumps, used as seed for generating data, pwcs
        Simparams; %simulation parameters
    end
    
    properties(Access = private)
        x_ins;
        SimulationMethodSteps;
        SimulationMethodNoise;
%         pwcsNt;
    end
    
    methods 
        function SDobj = SimData(varargin)
            addpath('internals/');
            
          if (nargin==1)
           SDobj.Simparams = varargin{1};
           SDobj.Simparams = SDobj.Simparams.updateRates();
          else
             SDobj.Simparams = SimParams;
          end
          SDobj.data = zeros(SDobj.Simparams.N,1);
          SDobj.pwcs = zeros(SDobj.Simparams.N,1);
          SDobj.time = ((0:1:(SDobj.Simparams.N-1)).*SDobj.Simparams.h)';
          SDobj.SimulationMethodSteps = 'none';
          SDobj.SimulationMethodNoise = 'none';
        end
        
        

        
        %call sd = sd.simulatePoissonSteps() to generate new poissonian
        %distributed noisy steps.
        function obj = simulatePoissonSteps(obj)
           
            [obj.pwcs, obj.p] = simMultiplePoisson(obj.Simparams.N,...
                obj.Simparams.poissonPars.jump_rate.*obj.Simparams.h,...
                obj.Simparams.poissonPars.delta,obj.Simparams.x0,...
                obj.Simparams.poissonPars.sigma_delta);
            
            obj.SimulationMethodSteps = 'PoissonSteps';
        end
        

            
              
        
        
        % call: sd = sd.simulatePol2(k1,kb,kf,kb1); 
        % K_D = 9.2; %NTP dissociation constant \mu M
        % c_NTP = 15; %NTP concentration in \mu M
        % k_1 = 88; %1/s
        % k_rev1 = 680; %1/s
        % k_3 = 35;
        % k2_net = c_NTP*k_3/K_D;
        % k1_net = k_1*k2_net/(k_rev1+k2_net);
        %k1 = 1/(1/k_3+1/k2_net+1/k1_net);
        function [obj, varargout] = simulatePol2(obj,varargin)
            
            if (nargin>3)
                simparams.k1 = varargin{1};
                simparams.kb = varargin{2};
                simparams.kf = varargin{3};
                simparams.kb1 = varargin{4};
                debugFlag = false;
                [obj.pwcs, obj.p, dwells] = simSimplePol2steps(obj.Simparams.N,...
                obj.Simparams.h,obj.Simparams.x0,simparams,debugFlag);
            
            else
                [simparams.k1,simparams.kb,simparams.kf,simparams.kb1] = ...
                    obj.Simparams.p2Pars.returnSimRates();
                [obj.pwcs, obj.p, dwells] = simSimplePol2steps(obj.Simparams.N,...
                obj.Simparams.h,obj.Simparams.x0,simparams);
            end

            obj.SimulationMethodSteps = 'Pol2Steps';
            if nargout>1
                varargout{1} = dwells;
            end
            
        end
        
        
        
        function [obj, varargout] = simulatePol2variableForce(obj,varargin)
            
            if nargin>1
                debugFlag = varargin{1};
            else
                debugFlag = false;
            end
            
            [obj.pwcs, obj.p, dwells] = simSimplePol2StepsincForce(obj.Simparams.N,...
                    obj.Simparams.h,obj.Simparams.x0,obj.Simparams.p2Pars.F,obj.Simparams.kappa,...
                    obj.Simparams.L,debugFlag,obj.Simparams.p2Pars.c_NTP);
            

            obj.SimulationMethodSteps = 'Pol2StepsvariableForce';
            if nargout>1
                varargout{1} = dwells;
            end
            
        end
        
        % note that this is applied after steps are simulated
        %call sd = sd.simOptTrapNoise()
        %or  sd = sd.simOptTrapNoise(1.0,pwcs_tmp,p_tmp)
        function [obj] = simOptTrapNoise(obj, varargin)
            k = obj.Simparams.k/obj.Simparams.gamma(1);
            D = obj.Simparams.kbT/obj.Simparams.gamma(1);
            if(nargin>1)
                noisefactor = varargin{1};
                D = noisefactor*obj.Simparams.kbT/obj.Simparams.gamma(1);
                fprintf('increased noise: %d \n',D);
            end
            
            if(nargin>2)
                
               pwcspts_tmp = varargin{2};
               p_old = varargin{3};      
               
               [pwcspts, pneu] = obj.reScaleJumpInput(pwcspts_tmp, p_old);
            
                [obj.data, obj.x_ins] = simNoisyData(obj.Simparams.N,obj.Simparams.h,k,D,...
                pwcspts,pneu,obj.Simparams.instrDrift,obj.Simparams.instrNoise,...
                obj.Simparams.x0,obj.Simparams.mix);
            else
                [obj.data, obj.x_ins] = simNoisyData(obj.Simparams.N,obj.Simparams.h,k,D,...
                obj.pwcs,obj.p,obj.Simparams.instrDrift,obj.Simparams.instrNoise,...
                obj.Simparams.x0,obj.Simparams.mix);
            end
            
            obj.SimulationMethodNoise = 'simpleTrap';
        end
        
        % note that this is applied after steps are simulated
        %call [sd,k,F] = sd.simOptTweezersVariableNoise();
        %or sd = sd.simOptTweezersVariableNoise(noisefactor,false,L0,F0);
        %sd = sd.simOptTweezersVariableNoise(noisefactor,false,L0,F0,simtype,pwcspts_tmp,p);
        function [obj, varargout] = simOptTweezersVariableNoise(obj,varargin)
            %input: constForce: false= opposing force mode, true=constant
            %Lstart: contour length of DNA
            %F_0: force at beginning
            
            
            if(nargin>2)
                constForce = varargin{2};
                Lstart = varargin{3};
                F_0 = varargin{4};
                if nargin>5
                    exp_type=varargin{5};
                else
                    exp_type=obj.Simparams.simtype;
                end
            else
                constForce = false;
                Lstart = obj.Simparams.L; 
                F_0 = obj.Simparams.p2Pars.F; 
                exp_type=obj.Simparams.simtype;
            end
            
            kapp = obj.Simparams.kappa(1);
            gam = obj.Simparams.gamma(1);
            D = obj.Simparams.kbT/gam;
            if(nargin>1)
                noisefactor = varargin{1};
                D = noisefactor*obj.Simparams.kbT/gam;
            end
            
            if(nargin>6)
                
               pwcspts = varargin{6};
               pneu = varargin{7};
%                pwcspts_tmp = varargin{6};
%                p_old = varargin{7};
%                
%                [pwcspts, pneu] = obj.reScaleJumpInput(pwcspts_tmp, p_old);
            
                [obj.data, obj.time,obj.x_ins, k, F]  = simDynamicNoise(obj.Simparams.N,obj.Simparams.h,kapp,gam,D,...
                obj.Simparams.instrDrift,obj.Simparams.instrNoise,...
                obj.Simparams.x0,obj.Simparams.mix,F_0,Lstart,pwcspts,pneu,constForce,exp_type);
            else
                [obj.data, obj.time,obj.x_ins, k, F]  = simDynamicNoise(obj.Simparams.N,obj.Simparams.h,kapp,gam,D,...
                obj.Simparams.instrDrift,obj.Simparams.instrNoise,...
                obj.Simparams.x0,obj.Simparams.mix,F_0,Lstart,obj.pwcs,obj.p,constForce,exp_type);
            end
            
            
            if nargout>1
               varargout{1} = k;
               varargout{2} = F;
            end
            
            obj.SimulationMethodNoise = 'variable';
        end
        
        
        
        
        function y = getInstrNoise(obj)
            if(strcmp(obj.SimulationMethodNoise,'simpleTrap'))
                y = obj.x_ins;
                
            elseif(strcmp(obj.SimulationMethodNoise,'variable'))
                y = obj.x_ins;
            
            else
                error('SimData::whichSimulation: simulation method jump diff not used!');
                y = [];
            end
        end
        
        
        function answ = whichSimulation(obj)
            answ.steps = obj.SimulationMethodSteps;
            answ.noise = obj.SimulationMethodNoise;
        end
        
        
        function plotData(obj,varargin)
        %call sd.plotData()
        
            figure_width  = 8*2;
            figure_height = 6*2;
            FontSize = 16;
            FontName = 'Helvetica';
        if(nargin>1)
            scale = varargin{1};
            dec = floor((1/obj.Simparams.h)/scale);
        else
            dec = floor((1/obj.Simparams.h)/100); % downsample to 100 Hz
        end
            figure;
            clf;
            set(gcf, 'units', 'centimeters', 'pos', [0 0 figure_width figure_height])
            % set(gcf, 'Units', 'pixels', 'Position', [100 100 500 375]);
            set(gcf, 'PaperPositionMode', 'auto');
            set(gcf, 'Color', [1 1 1]); % Sets figure background
            set(gca, 'Color', [1 1 1]); % Sets axes background
            set(gcf, 'Renderer', 'painters')
            
            plot1  = plot(obj.time,obj.data,'.b', ...
                          'Color', [0.6,0.6,0.6], ...
                          'MarkerSize', 5);
            hold on;
            plot2  = plot(decimate(obj.time,dec,'fir'),decimate(obj.data,dec,'fir'),'g','LineWidth',3);
            plot3 = plot(obj.time,obj.pwcs,'r','MarkerSize', 1, 'LineWidth', 2.5);
            hold off;
            %title('Generation of Noisy Data and Denoising');
            %legend([plot1 plot2 plot3],'noisy full bandwidth data set','box car avg noisy data','pure step signal');
            legend([plot1 plot2 plot3], ...
                    'noisy full bandwidth data set', ...
                    strcat('factor',num2str(dec),' , box car avg'), ...
                    'pure step signal', ...
                    'FontSize', floor(0.65 * FontSize), ...
                    'FontName', FontName, ...
                    'Location', 'NorthWest');
            xlabel('time/s', 'FontSize', FontSize, 'FontName', FontName);
            ylabel('bead motion/nm', 'FontSize', FontSize, 'FontName', FontName);
            
            set(gca, ...
              'Box'         , 'off'     , ...
              'TickDir'     , 'out'     , ...
              'TickLength'  , [.02 .02] , ...
              'XMinorTick'  , 'on'      , ...
              'YMinorTick'  , 'on'      , ...
              'YGrid'       , 'on'      , ...
              'XColor'      , [.3 .3 .3], ...
              'YColor'      , [.3 .3 .3], ...
              'LineWidth'   , 1         );
        end
        
        
        % rescales pwcs and p when Simparams.N and Simparams.h are changed
        function [pwcspts, pneu] = reScaleJumpInput(obj, pwcspts_tmp, p_old)
                    sizepwcs = length(pwcspts_tmp);
                    obj.time(end);
                    old_h = obj.time(end)/sizepwcs; 
                    fprintf('update data arrays sampled with delta_t=%d to new sampling freq. %d \n',old_h,obj.Simparams.h);
                    %old_ti =  ((0:1:(sizepwcs-1)).*old_h)';
               
               if(sizepwcs ~= obj.Simparams.N || old_h ~= obj.Simparams.h)
                   pwcspts = zeros(obj.Simparams.N ,1);
                   pneu = round(p_old.*old_h./obj.Simparams.h);
                   pneu = pneu(pneu<obj.Simparams.N);
                   pneu = pneu(pneu~=0);
                    if ~isempty(p_old);
                        m_p = min(length(p_old),length(pneu));
                         for i=2:m_p
                            if (p_old(i-1)~=0 && pneu(i-1)~=0)
                              pwcspts(pneu(i-1):pneu(i)) = pwcspts_tmp(p_old(i-1));
                            end
                         end
                     pwcspts(pneu(end):obj.Simparams.N)=pwcspts_tmp(p_old(m_p));
        
                    else
                        error('SimData:simulation','no steps generated please step points empty');
                    end
                    
               else
                   pwcspts = pwcspts_tmp;
                   pneu = p_old;
               end 
        end

        
    end
    
end


