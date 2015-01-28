function [x, pwcssteps, varargout] = simSimpleNoisyPol2(N,h,k,D,instrDrift,instrNoise,xstart,x0,a,varargin)
%SIMNOISYPOL2DYNAMICS Summary of this function goes here
%   Detailed explanation goes here


    
    if(nargin>9)
        simparams = varargin{1};
        k1 = simparams.k1;
        kb = simparams.kb;
        kf = simparams.kf;
        kb1 = simparams.kb1;
        [pwcssteps, p, pwcsNt] = simSimplePol2dynamics(N,h,x0,k1,kb,kf,kb1);
        
    else
        [pwcssteps, p, pwcsNt] = simSimplePol2dynamics(N,h,x0);
    end


    
    pwcssteps = 0.34.*pwcssteps;%rescale to nm;
    pwcsNt.DNA = 0.34*pwcsNt.DNA;
    pwcsNt.RNA = 0.34*pwcsNt.RNA;
    
    x = zeros(N,1);
    y = zeros(N,1);
    if(a==0)
        y = zeros(N,1);
    else
        y = simOU(N,h,instrDrift,instrNoise,0.0,0.0); %Instrument fluctuations
    end
    
    if ~isempty(p)
    
    x(1:p(1)) = simOU(p(1),h,k,D,xstart,x0);
    %x(1:p(1)) = x(1:p(1)) + a.*y(1:p(1));
        for i=2:length(p)
            x((p(i-1)+1):p(i))= simOU((p(i)-p(i-1)),h,k,D,x(p(i-1)),pwcsNt.DNA(p(i-1)));
            %x((p(i-1)+1):p(i)) = x((p(i-1)+1):p(i))+a.*y((p(i-1)+1):p(i));
        end
    x((p(end)+1):N) = simOU((N-p(end)),h,k,D,x(p(end)),pwcsNt.DNA(p(end)));
    %x((p(end)+1):N) = x(p(end)+1:N) + a.*y(p(end)+1:N);
    end
    x = x + a.*y;
    
    switch nargout
        
    case 3 %return rna position
        varargout{1} = pwcsNt;
    
    case 4 % time axis is expected as additional output argument
        varargout{1} = pwcsNt;
        ti = 0:1:(N-1);
        ti = ti.*h;
        varargout{2} = ti;
    case 6 %time axis and intrumental noise is expected as additional output
        varargout{1} = pwcsNt;
        ti = 0:1:(N-1);
        ti = ti.*h;
        varargout{2} = ti;
        y = a.*y;
        varargout{3} = y;
        varargout{4} = p;
    end
    
    

end



