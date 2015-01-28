function [pwcs varargout] = simMultiplePoisson(N,jump_rate,delta,x0,varargin)
%SIMMULTIPLEPOISSON Summary of this function goes here
%   output:
%       pwcs: Nx1 array of piece wise constant plateaus of poisson
%       processes
%       optional:
%           p: temporal positions of steps p(length(p))<N of the last generated
%               process
%
%   input:
%       N:      number of data points
%       jump_rate: array of jump probabilities
%       delta:  jump height
%       x0:     initial position
%   optional:
%       sigma: standard deviation of jump height

rng('shuffle'); %randomize seeds to current time

%determine number of jump processes
numOfProc = length(jump_rate);
% if more jump processes are expected than step sizes are given fill up the missing ones 
if (length(delta)<numOfProc)
   nDelta = length(delta);
   delta(nDelta:numOfProc)=delta(nDelta); 
end



pwcstemp = zeros(N,numOfProc);
pcell = cell(numOfProc,1);

if nargin > 4
    
    sigma = varargin{1};
    
        if (length(sigma)<numOfProc)
            nSigma = length(sigma);
            sigma(nSigma:numOfProc)=delta(nSigma); 
        end
    
    for j=1:numOfProc
        [pwcstemp(:,j),ptemp] = simPoisson(N,jump_rate(j),delta(j),x0,sigma(j));
        pcell{j} = ptemp; 
    end
    pwcs = sum(pwcstemp,2);
else
    
    
    for j=1:numOfProc
        [pwcstemp(:,j),ptemp] = simPoisson(N,jump_rate,delta,x0);
        pcell{j} = ptemp; 
    end
    pwcs = sum(pwcstemp,2);
end

p = sort(cell2mat(pcell));

if (nargout > 1)    % step points are expected as additional output argument
    varargout(1) = {p};
end

    
end

