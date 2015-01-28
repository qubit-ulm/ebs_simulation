function [x] = simOU(N,h,k,D,xinit,x0)
%SIMOU simulates an Ornstein-Uhlenbeck process
%   output: 
%           x:  Nx1 array of OU distributed random numbers
%   input:
%           N:  number of datapoints
%           h:  time increments delta_t
%           k:  corner frequency of harmonic potential (units inverse to h)
%           D:  Diffusion constant of brownian motion
%           x0: initial position (elongation from minimum of harmonic
%                       potential)

    x = zeros(N,1);
    %random diffusion term xi(t_n+1)-xi(t_n)=sqrt(2*kbT*h/gamma)*z where z is standard
    %normal distributed N(0,1).
    %pick gaussian iid random numbers:
    var = eye(1);
    var(1,1)=sqrt(2*D/k);
    xi = mvnrnd(0.0,var,N);
%    fprintf('simOU');
%     xinit
%     x0
%     x(1)
%     
%     x(1) = xinit-x0;
%     x(1)
x(1) = xinit;
    for i=1:N-1
        x(i+1) = (x(i)-x0)*exp(-k*h)+sqrt(1-exp(-2*k*h))*xi(i)+x0;
    end
    x = x;%+ x0;
end

