classdef Utest_dataGeneration < matlab.unittest.TestCase
    %UTEST_NOISEGENERATION tests dataSimulation functions from SimData
    %class
    
    properties (TestParameter)
        poissonRates = struct('slow',0.1,'fast',10,'vector',[1,5]);
        ntp_conc = struct('slow',1,'medium',100,'fast',2000);
        noisefactor = struct('smaller',0.1,'not_modified',1.0,'bigger',10);
        force = struct('assist',6.5,'oppose',-0.1);
        
    end
    
    methods (Test)
        %test if poisson step generation works
        function test_poissonStepGeneration(testCase,poissonRates)
            sp = SimParams();
            sp.N = 1e5;
            sp.h = 5e-4;
            sp.poissonPars.jump_rate = poissonRates;
            sd = SimData(sp);
            sd = sd.simulatePoissonSteps();
            numsteps = sum(diff(sd.pwcs)~=0);
            fprintf('number of poisson steps: %i \n',numsteps);
            verifyNotEqual(testCase,numsteps,0);
        end
        
        %check if pol2 type step generation works
        function test_pol2StepGeneration(testCase, ntp_conc)
            sp = SimParams();
            sp.N = 1e5;
            sp.h = 5e-4;
            sp.p2Pars.c_NTP = ntp_conc;
            sp = sp.updateRates();
            sd = SimData(sp);
            sd = sd.simulatePol2();
            [k1,kb,kf,kb1] = sp.p2Pars.returnSimRates(); 
            numsteps = sum(diff(sd.pwcs)~=0);
            fprintf('number of steps: %i elongation rate: %d \n',numsteps,k1);
            verifyNotEqual(testCase,numsteps,0);
        end
        
        function test_pol2VarForceStepGeneration(testCase, ntp_conc,force)
            sp = SimParams();
            sp.N = 1e5;
            sp.h = 5e-4;
            sp.p2Pars.F = force;
            sp.p2Pars.c_NTP = ntp_conc;
            sd = SimData(sp);
            sd = sd.simulatePol2variableForce();
            numsteps = sum(diff(sd.pwcs)~=0);
            fprintf('number of steps: %i of variable force \n',numsteps);
            verifyNotEqual(testCase,numsteps,0);
        end
        
        %check if simple noise generation works
        function test_simpleNoiseModel(testCase, noisefactor)
            sp = SimParams();
            sp.N = 1e5;
            sp.h = 5e-4;
            sp = sp.updateRates();
            sd = SimData(sp);
            sd = sd.simulatePoissonSteps();
            sd = sd.simOptTrapNoise(noisefactor);
            stdnoise = std(sd.data);
            verifyNotEqual(testCase,stdnoise,0);
        end
        
        %check if advanced, variable noise generation works
        function test_advancedNoiseModel(testCase, noisefactor)
            sp = SimParams();
            sp.N = 1e5;
            sp.h = 5e-4;
            sp = sp.updateRates();
            sd = SimData(sp);
            sd = sd.simulatePoissonSteps();
            sd = sd.simOptTweezersVariableNoise(noisefactor);
            stdnoise = std(sd.data);
            verifyNotEqual(testCase,stdnoise,0);
        end
        
        %check simpleNoiseModel: modify std of noise
        function testNoiseMod_simpleNoiseModel(testCase)
            sp = SimParams();
            sp.N = 1e5;
            sp.h = 5e-4;
            sp = sp.updateRates();
            sd = SimData(sp);
            sd = sd.simulatePoissonSteps();
            sd = sd.simOptTrapNoise(0.1);
            stdnoise(1) = std(sd.data-sd.pwcs);
            sd = sd.simOptTrapNoise(10);
            stdnoise(2) = std(sd.data-sd.pwcs);
            verifyThat(testCase, stdnoise(1)<stdnoise(2), matlab.unittest.constraints.IsTrue());
        end
        
        %check advancedNoiseModel: modify std of noise
        function testNoiseMod_advancedNoiseModel(testCase)
            sp = SimParams();
            sp.N = 1e5;
            sp.h = 5e-4;
            sp = sp.updateRates();
            sd = SimData(sp);
            sd = sd.simulatePoissonSteps();
            sd = sd.simOptTweezersVariableNoise(0.1);
            stdnoise(1) = std(sd.data-sd.pwcs);
            sd = sd.simOptTweezersVariableNoise(10);
            stdnoise(2) = std(sd.data-sd.pwcs);
            verifyThat(testCase, stdnoise(1)<stdnoise(2), matlab.unittest.constraints.IsTrue());
        end
        
        function test_rescaleGeneratedSteps(testCase)
            sp = SimParams();
            sp.N = 1e5;
            sp.h = 5e-4;
            sp = sp.updateRates();
            sd = SimData(sp);
            sd = sd.simulatePoissonSteps();
            sd.Simparams.h = 1e-3;
            sd.Simparams.N = 5e4;
            [pwcsneu, pneu] = sd.reScaleJumpInput(sd.pwcs, sd.p);
            verifyThat(testCase, length(pwcsneu)<length(sd.pwcs), matlab.unittest.constraints.IsTrue());
        end
        
        function test_rescaleGeneratedNoiseSimple(testCase)
            sp = SimParams();
            sp.N = 1e5;
            sp.h = 5e-4;
            sp = sp.updateRates();
            sd = SimData(sp);
            sd = sd.simulatePoissonSteps();
            sd = sd.simOptTrapNoise();
            datalength(1) = length(sd.data);
            sd.Simparams.h = 1e-3;
            sd.Simparams.N = 5e4;
            [pwcsneu, pneu] = sd.reScaleJumpInput(sd.pwcs, sd.p);
            sd = sd.simOptTrapNoise(1.0,pwcsneu,pneu);
            datalength(2) = length(sd.data);
            verifyThat(testCase, datalength(2)<datalength(1), matlab.unittest.constraints.IsTrue());
        end
        
        function test_rescaleGeneratedNoiseAdvanced(testCase)
            sp = SimParams();
            sp.N = 1e5;
            sp.h = 5e-4;
            sp = sp.updateRates();
            sd = SimData(sp);
            sd = sd.simulatePoissonSteps();
            sd = sd.simOptTweezersVariableNoise();
            datalength(1) = length(sd.data);
            sd.Simparams.h = 1e-3;
            sd.Simparams.N = 5e4;
            [pwcsneu, pneu] = sd.reScaleJumpInput(sd.pwcs, sd.p);
            sd = sd.simOptTweezersVariableNoise(1.0,false,sd.Simparams.L,...
                sd.Simparams.p2Pars.F,sd.Simparams.simtype,pwcsneu,pneu);
            datalength(2) = length(sd.data);
            verifyThat(testCase, datalength(2)<datalength(1), matlab.unittest.constraints.IsTrue());
        end
    end
    
end

