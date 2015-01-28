function [x, varargout] = simNoisyData(N,h,k,D,pwcs,p,instrDrift,instrNoise,x0,a)
%SIMNOISYDATA simulates an Ornstein-Uhlenbeck (OU) process of confined brownian motion
%which is interrupted by sudden jumps. In Addition there is a second OU
%processes which simulates instrumentNoise, modeled by instrDrift and
%instrNoise.
%   output: 
%           x:  Nx1 array of OU distributed random numbers
%           pwcs: piecewise constant signals
%       optional:
%           ti: time axis
%           y: instument noise
%   input:
%           N:  number of datapoints
%           h:  time increments delta_t
%           k:  corner frequency of harmonic potential (units inverse to h)
%           D:  Diffusion constant of brownian motion
%           instrDrift: rate which describes drift typical ms to s regime
%           instrNoise: amplitude of gaussian noise fluctuation
%           x0: initial position (elongation from minimum of harmonic
%                       potential)
%           a: mixing factor of brownian motion of beads and instrument
%           Noise
%           varargin{1}: p,pwcs: results of poisson process simulation as
%           input


if (isempty(p) || isempty(pwcs))
    error('simNoisyData:argChk', 'simulate steps at first before generating noise');
end

x = zeros(N,1);

y = simOU(N,h,instrDrift,instrNoise,0.0,0.0); %Instrument fluctuations

xstart = sqrt(D/k)*randn(1); %random start point

if ~isempty(p)
    
    x(1:p(1)) = simOU(p(1),h,k,D,xstart,x0);
    x(1:p(1)) = x(1:p(1)) + a.*y(1:p(1));
    for i=2:length(p)
        x((p(i-1)+1):p(i)) = simOU((p(i)-p(i-1)),h,k,D,x(p(i-1)),pwcs(p(i-1)));
        x((p(i-1)+1):p(i)) = x((p(i-1)+1):p(i))+a.*y((p(i-1)+1):p(i));
    end
    x((p(end)+1):N) = simOU((N-p(end)),h,k,D,x(p(end)),pwcs(p(end)));
    x((p(end)+1):N) = x(p(end)+1:N) + a.*y(p(end)+1:N);
end

 switch nargout
    
    case 2 % time axis is expected as additional output argument
        ti = 0:1:(N-1);
        ti = ti.*h;
        y = a.*y;
        varargout{1} = y;
    case 3 %time axis and intrumental noise is expected as additional output
        ti = 0:1:(N-1);
        ti = ti.*h;
        y = a.*y;
        varargout{1} = y;
        varargout{2} = ti;
        
     case 5
         ti = 0:1:(N-1);
        ti = ti.*h;
        y = a.*y;
        varargout{1} = y;
        varargout{2} = ti;
        varargout{3} = p;
        varargout{4} = pwcs;
 end



end

