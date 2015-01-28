function [pwcs, p, varargout] = simSimplePol2incForce(N,h,x0,F_start,kappa,L,varargin)

%SIMSIMPLEPOL2INCFORCE simulates a simplified Pol2 model according to 
%Dangkulwanich et.al. eLife 2013;2:e00971
%   Input:
%       N: number of data points
%       h: sampling rate (total measurement time T=h*N)
%       x0: offset (start value)
%       Fstart: initial force of constant trap distance experiment
%       varargin: rate constants: k_elong, k_b, k_f, k_b1
%   Output:
%       pwcs: piecewise constant signal of length N
%       p: stepping times as array positions in pwcs
%       varagout:
%           dwells.elongrate: elongation dwell times
%           dwells.pause_entryrate: dwell times of pause entry
%           pwcsNt: DNA and RNA state
    
%     if(F_start<0)
%         error('simSimplePol2incForce:negativeNumber','F_start has to be positive');
%     end
    
%     %compute extension length at Force F_start
%     x0 = computeExtAtForce(F_start,L,100);
%     fprintf('starting Force=%d, and extension x0=%d \n',F_start,x0);

    
elongrate = [];
pause_entryrate = [];

rng('shuffle'); %randomize seeds

switch nargin
    case 7
        debug = varargin{1};
        c_NTP = 1000;
    case 8
        debug = varargin{1};
        c_NTP = varargin{2};
    otherwise
        debug = false;
        c_NTP = 1000;
end


%read simulation parameters
[k_elong,k_b1,k_f,k_b] = getUpdatedRates(F_start,h,c_NTP);


%calculate transition probabilities 
p1 = k_elong/(k_b1+k_elong);


cstate = 'onpath';
statPos(:,1) = [0.0;0.0;0.0];
j = 1;
t = 0;

while(t<N)
    
    j = j+1;
    
    %state machine for implementing on and off pathway translocation of
    %PolII
    switch cstate
        
        case 'onpath'
            if(debug)
            fprintf('%i now(%d): %s DNA pos %d , RNA pos %d \n',j,t*h,cstate,statPos(2,j-1),statPos(3,j-1));
            end
            
            if(rand(1,1)<=p1)
                t1 = -log(rand(1,1))./k_elong;
                state = 'onpath';elongrate = [elongrate, t1];
                statPos(:,j) = [t+t1;statPos(2,j-1)+1;statPos(3,j-1)+1];
                %new forces and rates:
                Fnew = updateForces(F_start,statPos(2,j),kappa,L,debug);
                [k_elong,k_b1,k_f,k_b] = getUpdatedRates(Fnew,h,c_NTP);
            else
                tb1 = -log(rand(1,1))./k_b1;
                state = 'backtrack1';pause_entryrate = [pause_entryrate, tb1];
                statPos(:,j) = [t+tb1;statPos(2,j-1)-1;statPos(3,j-1)];
                %new forces and rates
                Fnew = updateForces(F_start,statPos(2,j),kappa,L,debug);
                [k_elong,k_b1,k_f,k_b] = getUpdatedRates(Fnew,h,c_NTP);
            end

        case 'backtrack1'
            if(debug)
            fprintf('%i now(%d): %s DNA pos %d , RNA pos %d \n',j,t*h,cstate,statPos(2,j-1),statPos(3,j-1));
            end
            t_f = -log(rand(1,1))./k_f;
            t_b = -log(rand(1,1))./k_b;
            if(t_f<t_b)
                state = 'onpath';
                statPos(:,j) = [t+t_f;statPos(2,j-1)+1;statPos(3,j-1)];
                %new forces and rates
                Fnew = updateForces(F_start,statPos(2,j),kappa,L,debug);
                [k_elong,k_b1,k_f,k_b] = getUpdatedRates(Fnew,h,c_NTP);
            else
                state = 'backtrackn';
                statPos(:,j) = [t+t_b;statPos(2,j-1)-1;statPos(3,j-1)];
                %new forces and rates
                Fnew = updateForces(F_start,statPos(2,j),kappa,L,debug);
                [k_elong,k_b1,k_f,k_b] = getUpdatedRates(Fnew,h,c_NTP);
            end
            
        case 'backtrackn'
            if(debug)
            fprintf('%i now(%d): %s DNA pos %d , RNA pos %d \n',j,t*h,cstate,statPos(2,j-1),statPos(3,j-1));
            end
            if((statPos(2,j-1)-statPos(3,j-1))>-3)
                forward = 'backtrack1';
            else             
                forward = 'backtrackn';
            end
            t_f = -log(rand(1,1))./k_f;
            t_b = -log(rand(1,1))./k_b;
            if(t_f<t_b)
                state = forward;
                statPos(:,j) = [t+t_f;statPos(2,j-1)+1;statPos(3,j-1)];
                %new forces and rates
                Fnew = updateForces(F_start,statPos(2,j),kappa,L,debug);
                [k_elong,k_b1,k_f,k_b] = getUpdatedRates(Fnew,h,c_NTP);
            else
                state = 'backtrackn';
                statPos(:,j) = [t+t_b;statPos(2,j-1)-1;statPos(3,j-1)];
                %new forces and rates
                Fnew = updateForces(F_start,statPos(2,j),kappa,L,debug);
                [k_elong,k_b1,k_f,k_b] = getUpdatedRates(Fnew,h,c_NTP);
            end
            
        otherwise
            
    end
    
    t = statPos(1,j);
    cstate = state;
       
end



tpts = round(statPos(1,:))+1;
p = tpts;
%remove identical time points
tpts = (diff(tpts)>0);
p = p(tpts);
pwcsstep = statPos(2,:);
pwcsstep = pwcsstep(tpts);
rnastep = statPos(3,:);
rnastep = rnastep(tpts);

%remove points which are out of index interval
p = p(p<N);
p = p(p~=0);

pwcs = zeros(N,1);
pwcsRNA = zeros(N,1);

if ~isempty(p);
     for i=2:length(p)
            if (p(i-1)~=0)
              pwcs(p(i-1):p(i)) = pwcsstep(i-1);
              pwcsRNA(p(i-1):p(i)) = rnastep(i-1);
            end
     end
        pwcs(p(length(p)):N)=pwcsstep(length(p));
        pwcsRNA(p(length(p)):N)=rnastep(length(p));      
else
        error('simPoisson:no steps generated please increase jump_rate and try again');  
end

pwcs = pwcs.*0.34 +x0;
pwcsRNA = pwcsRNA + x0;

pwcsNt.DNA = pwcs;
pwcsNt.RNA = pwcsRNA;
dwells.elongation = elongrate;
dwells.pause_entry = pause_entryrate;

if(nargout>2)
    varargout{1} = dwells;
    varargout{2} = pwcsNt;
end



end


function [k_elong,k_b1,k_f,k_b] = getUpdatedRates(F,h,varargin)

   if nargin>2
       c_NTP = varargin{1};
   else
       c_NTP = 1000; %NTP concentration in \mu M
   end
   kbT = 4.11; %pN*nm (298K)
   delta = 0.17; %nm 0.5bp
   k0 = 1.3; %1/s
   k0b1 = 5;%6.9;
   k_b1 = k0b1.*exp(-F.*(0.34-delta)./(kbT));
   k_b = k0.*exp(-F.*(0.34-delta)./(kbT));
   k_f = k0.*exp(F.*delta./(kbT));
   %parameters describing NTP binding rate
   K_D = 9.2; %NTP dissociation constant \mu M
   
   k_1 = 88*exp(F.*delta./(kbT)); %1/s
   k_rev1 = 680*exp(-F.*(0.34-delta)./(kbT)); %1/s
   k_3 = 35; %1/s
   %approximation under NTP concentrations used 35-2000 \mu M k_rev2,k2>>k3
   %so the ireversible rate k3 can be approximated by a net rate 
   k2_net = c_NTP*k_3/K_D;
   k1_net = k_1*k2_net/(k_rev1+k2_net); 
   k_elong = 1/(1/k1_net+1/k2_net+1/k_3);
   
    k_f = k_f*h;
    k_b = k_b*h;
    k_elong = k_elong*h;
    k_b1 = k_b1*h;
end

function Fnew = updateForces(Fstart,displacement,trap_stiff,L,debug)
    if(Fstart>0)
        Fnew = max(Fstart-trap_stiff*displacement*0.34,0);
    else
        if((L-displacement*0.34)>0)
            Fnew = Fstart-trap_stiff*displacement*0.34;
        else
            Fnew = 0.0;
            fprintf('elongation exceeded maximum length at step position=%d. Set force to zero \n',displacement);
        end
    end
    if(debug)
        fprintf('new force = %d \n',Fnew);
    end
end

