function [x, varargout]  = simDynamicNoise(N,h,kappa,gamma,D,instrDrift,instrNoise,x0,a,F_0,Lstart,pwcs,p,const_mode,varargin)
%SIMDYNAMICNOISE Simulates poissonian steps and takes into account
%increasing noise level through shortening of tether and increasing tether
%stiffness
%   input: pwcs,p
%       F_0: initial force
%       Lstart: initial length of DNA tether
%       const_mode: true if constant force
    if nargin>15
        exp_type = varargin{1};
    else
        exp_type = 'oppose';
    end

    if ~isempty(p)
        if(strcmp(exp_type,'assist'))
            L=Lstart+pwcs(p);
            L(L<0) = 0;
            %fprintf('DNA tether length: %d \n',L);
        else
            L = Lstart-pwcs(p);
            L(L<0) = 0; %enzyme cannot transcripe more than template
        end
    
        if(const_mode)
            F = F_0*ones(1,length(p));
        else
            if(strcmp(exp_type,'assist'))
                F=-kappa*pwcs(p) + F_0;
                F(F<0) = 0;
            else
                F=kappa*pwcs(p) + F_0;
            end
        end
    else
        F=F_0;
        L = Lstart;
    end
    
    %compute increasing DNA stiffness for different Forces
    k = calculateDNAstiffnessWLC(F,L,1000,0.34);
    %fprintf('k = %d \n',k);
    %fprintf('F = %d \n',F);
    k = 2.*k + kappa;
    k = k./gamma;
    
    %simulate noise according to parameters
    x = zeros(N,1);
    if(a==0)
        y = zeros(N,1);
    else
        y = simOU(N,h,instrDrift,instrNoise,0.0,0.0); %Instrument fluctuations
    end
%     fprintf('sqrt(2*D/k(1)) = %d \n',sqrt(2*D/k(1)));

    xstart = sqrt(D(1)/k(1))*randn(1); %random start point

    if ~isempty(p)
    x(1:p(1)) = simOU(p(1),h,k(1),D,x0+xstart,x0);
    x(1:p(1)) = x(1:p(1)) + a.*y(1:p(1));
        for i=2:length(p)
%             fprintf('k(%i) = %d, sqrt(2*D/k(%i)) = %d \n',i,k(i),i,sqrt(2*D/k(i)));
            x((p(i-1)+1):p(i))= simOU((p(i)-p(i-1)),h,k(i),D,x(p(i-1)),pwcs(p(i-1)));
            x((p(i-1)+1):p(i)) = x((p(i-1)+1):p(i))+a.*y((p(i-1)+1):p(i));
        end
    x((p(end)+1):N) = simOU((N-p(end)),h,k(end),D,x(p(end)),pwcs(p(end)));
    x((p(end)+1):N) = x(p(end)+1:N) + a.*y(p(end)+1:N);
    end
    
    switch nargout
        
    
    case 2 % time axis is expected as additional output argument
        
        ti = 0:1:(N-1);
        ti = ti.*h;
        varargout{1} = ti;
    case 3 %time axis and intrumental noise is expected as additional output
        
        ti = 0:1:(N-1);
        ti = ti.*h;
        varargout{1} = ti;
        y = a.*y;
        varargout{2} = y;
        
    case 5
        
        ti = 0:1:(N-1);
        ti = ti.*h;
        varargout{1} = ti;
        y = a.*y;
        varargout{2} = y;
        varargout{3} = k.*gamma;
        varargout{4} = F;
    end

end

