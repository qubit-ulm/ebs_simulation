function [x, pwcssteps, varargout]  = simSimplePol2incForce(N,h,kappa,gamma,D,F_start,L,instrDrift,instrNoise,xstart,x0,a,varargin)
%SIMSIMPLEPOL2INCFORCE Summary of this function goes here
%   input:
%       F_start: Force between tethered beads when simulation is started
%       L: Contour length of DNA
    
    if(F_start<0)
        error('simSimplePol2incForce:negativeNumber','F_start has to be positive');
    end
    
    %compute extension length at Force F_start
    x0 = computeExtAtForce(F_start,L,100);
    fprintf('starting Force=%d, and extension x0=%d \n',F_start,x0);

    if(nargin>12)
        simparams = varargin{1};
        k1 = simparams.k1;
        kb = simparams.kb;
        kf = simparams.kf;
        kb1 = simparams.kb1;
        [pwcssteps, p, F,pwcsNt] = simPol2dynamicsincForce(N,h,F_start,L,k1,kb,kf,kb1);
        
    else
        [pwcssteps, p, F,pwcsNt] = simPol2dynamicsincForce(N,h,F_start,L,kappa);
    end


    
    pwcssteps = 0.34.*pwcssteps+x0;%rescale to nm;
    pwcsNt.DNA = 0.34*pwcsNt.DNA+x0;
    pwcsNt.RNA = 0.34*pwcsNt.RNA+x0;
    
    %compute increasing DNA stiffness
    k = computeDNAStiffness(F,L,1000);
    k = 2.*k + kappa;
    k = k./gamma;
    
    
    x = zeros(N,1);
    
    if(a==0)
        y = zeros(N,1);
    else
        y = simOU(N,h,instrDrift,instrNoise,0.0,0.0); %Instrument fluctuations
    end
    
    if ~isempty(p)
    x(1:p(1)) = simOU(p(1),h,k(1),D,x0+xstart,x0);
    x(1:p(1)) = x(1:p(1)) + a.*y(1:p(1));
        for i=2:length(p)
            x((p(i-1)+1):p(i))= simOU((p(i)-p(i-1)),h,k(i),D,x(p(i-1)),pwcsNt.DNA(p(i-1)));
            x((p(i-1)+1):p(i)) = x((p(i-1)+1):p(i))+a.*y((p(i-1)+1):p(i));
        end
    x((p(end)+1):N) = simOU((N-p(end)),h,k(end),D,x(p(end)),pwcsNt.DNA(p(end)));
    x((p(end)+1):N) = x(p(end)+1:N) + a.*y(p(end)+1:N);
    end
    
    switch nargout
        
    case 3 %return rna position
        varargout{1} = pwcsNt;
    
    case 4 % time axis is expected as additional output argument
        varargout{1} = pwcsNt;
        ti = 0:1:(N-1);
        ti = ti.*h;
        varargout{2} = ti;
    case 5 %time axis and intrumental noise is expected as additional output
        varargout{1} = pwcsNt;
        ti = 0:1:(N-1);
        ti = ti.*h;
        varargout{2} = ti;
        y = a.*y;
        varargout{3} = y;
    end



end


function x0 = computeExtAtForce(F_start,L,startval)

    fun_Force = @(x)(calcForce_WLCModel(x,L)-F_start);
    
    x0 = fzero(fun_Force,startval);
end

function f = calcForce_WLCModel(x,L)
    %input:
    %   x in nm
    %   L in nm, contour length of DNA

    P = 50; %persistance length in nm
    f = (4.1/(P*4)*(1./(1-x./L).^2 - 1 + 4.*x./L));
end

function kstiff = computeDNAStiffness(F,L,startval)

    kstiff = zeros(1,length(F));
    for i=1:length(F)
        x = computeExtAtForce(F(i),L,startval);
        kstiff(i) = stiffnessDNA(x,L);
    end
end

function kDNA = stiffnessDNA(x,L)
    kDNA = 4.1/(50*4)*(2./(L.*(1-x./L).^3) + 4/L);
end