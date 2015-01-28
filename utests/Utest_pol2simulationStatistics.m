classdef Utest_pol2simulationStatistics < matlab.unittest.TestCase
    %UTEST_POL2SIMULATIONSTATISTICS unit test suite to check statistical
    %   properties of the step generation functions of the pol2 simulation
    
    properties (TestParameter)
        numdata = struct('small',5e4,'medium',4e5,'large',2e6);
        ntp_conc = struct('very_low',10,'low',25,'medium',100,'high',500,'saturating',2000);
    end
    
    methods (Test)
        function testPausefreeVelocity(testCase, ntp_conc)
            sp = SimParams();
            sp.N = 5e5;
            sp.h = 5e-4;
            sp.p2Pars.c_NTP = ntp_conc;
            sp = sp.updateRates();
            sd = SimData(sp);
            sd = sd.simulatePol2();
            logic_backtracks = extractBacktracks(sd.time,sd.pwcs);
            transloc = sd.pwcs(logic_backtracks<1); 
            forward_rate = 1/mean(diff(sd.time(abs(diff(transloc))>0)));
            numelong = sum(diff(transloc)>0);
            relativ_tol = 0.15;%min(10/sqrt(numelong),0.1);
            fprintf('pause free velocity = %.2f , input elongation rate = %.2f, relativ tol=%d, number of elongation steps =%i \n',forward_rate,sp.p2Pars.k1,relativ_tol,numelong);
            verifyEqual(testCase, sp.p2Pars.k1, forward_rate, 'RelTol', relativ_tol, ...
                'Difference between actual and expected exceeds relative tolerance')
            %Tolerance Definition:
            %        abs(expected - actual) <= tolerance .* abs(expected)
        end
        
        function testPauseDensity(testCase, numdata, ntp_conc)
            sp = SimParams();
            sp.N = numdata;
            sp.h = 5e-4;
            sp.p2Pars.c_NTP = ntp_conc;
            sp = sp.updateRates();
            sd = SimData(sp);
            [sd,dwells] = sd.simulatePol2();
            pauseDensity = length(dwells.pause_entry)/(length(dwells.elongation)+length(dwells.pause_entry));
            [k1,kb,kf,kb1] = sp.p2Pars.returnSimRates();
            probabPauseEntry = kb1/(k1+kb1);
            
            relativ_tol = 0.15;%min(10/sqrt(numelong),0.1);
            fprintf('pause density = %.2f , input elongation rate = %.2f \n',pauseDensity,probabPauseEntry);
            verifyEqual(testCase, probabPauseEntry, pauseDensity, 'RelTol', relativ_tol, ...
                'Difference between actual and expected exceeds relative tolerance')
        end
        
    end
    
    
    
end



  
 
