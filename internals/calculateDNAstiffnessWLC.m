function kstiff = calculateDNAstiffnessWLC(F,L,startval,delta)
%input: F: array of forces
%       L: length of DNA-tether between beads (shortened/prolonged during
%       transcription)
%       delta: smallest length scale (e.g. 0.34nm), to search x(F) in [0,L-delta]
%output: kstiff: array of DNA stiffness

    

    kstiff = zeros(1,length(F));
    for i=1:length(F)
        %take care that startvalue is shorter than contour length of DNA
        while startval>(L(i)*2/3)
            startval = startval-1;
        end
        %fprintf('new start value for root search= %d \n',startval);
        x = computeExtAtForce(F(i),L(i),startval,delta);
%         fprintf('x = %d \n',x);
        kstiff(i) = stiffnessDNA(x,L(i));
        %fprintf('kstiff(i) = %d \n',kstiff(i));
    end
end

function x0 = computeExtAtForce(F_start,L,startval,delta)

    fun_Force = @(x)(calcForce_WLCModel(x,L)-F_start);
    
    x0 = fzero(fun_Force,[0,L-delta]);
end

function f = calcForce_WLCModel(x,L)
    %input:
    %   x in nm
    %   L in nm, contour length of DNA

    P = 50; %persistance length in nm
    kbT = 4.1; %thermal energy
    %force according to WLC-model
    f = (kbT/(P*4)*(1./(1-x./L).^2 - 1 + 4.*x./L));
end



% force dependent stiffness according to WLC-model
function kDNA = stiffnessDNA(x,L)
    P = 50; %persistance length in nm
    kbT = 4.1; %thermal energy
    kDNA = kbT/(P*4)*(2./(L.*(1-x./L).^3) + 4/L);
end