function [pwcs, p] = simPoisson(N,jump_rate,delta,x0,varargin)
%SIMPOISSON simulates a 1D poisson process with rate lambda
%   output:
%       p: temporal positions of steps p(length(p))<N
%       orig: Nx1 array of piece wise constant plateaus of poisson process
%
%   input:
%       N:      number of data points
%       lambda: rate
%       delta:  jump height
%       x0:     initial position
%       optional:
%           sigma: standard deviation of step size in [nm];
    
rng('shuffle'); %randomize seeds to current time

if nargin>4
    sigma = varargin{1};
else
    sigma = 0;
end

if jump_rate<0
    jump_rate = abs(jump_rate);
end

    pwcs = rand(N,1);
    pwcs = -log(pwcs)./jump_rate;
    pwcs = floor(cumsum(pwcs));
    p = pwcs(pwcs<N); %only keep step points which are smaller than N
    p =p(p>0);
    
    %try to draw new random numbers there are no steps within the interval
    while isempty(p)
        
        pwcs = rand(N,1);
        pwcs = -log(pwcs)./jump_rate;
        pwcs = floor(cumsum(pwcs));
        p = pwcs(pwcs<N);
    end
    
    pwcs = zeros(N,1);
    if ~isempty(p);
        for i=2:length(p)
                if (p(i-1)~=0)
                    %p((i-1):i)
                 pwcs(p(i-1):p(i)) = (i-2)*delta+sigma*randn(1);
                end
        end
        
        pwcs(p(length(p)):N)=length(p)*delta;
        
    else
        error('simPoisson:no steps generated please increase jump_rate and try again');
    
    end
    pwcs = pwcs +x0;
    
   
    
end

