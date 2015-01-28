classdef Utest_getParameters  < matlab.unittest.TestCase
    %UTEST_GETPARAMETERS checks simulation parameter classes
    %   Detailed explanation goes here
    
    properties (TestParameter) 
        
        numdata = struct('small',1e4,'medium',1e5,'large',1e6);
    end
    
    methods (Test)
        
        %check if PoissonParams() exists
        function test_existPoissonParameter(testCase)
            poisspara = PoissonParams();
            testCase.verifyClass(poisspara, 'PoissonParams');
        end
        function test_getPoissonParameter(testCase)
            poisspara = PoissonParams();
            poisson_out = poisspara.getSimpleParams();
            testCase.verifyClass(poisson_out.lambda, 'double');
        end
        
        
        %check if Pol2Params() exists
        function test_existPol2Parameter(testCase)
            p2pars = Pol2Params();
            testCase.verifyClass(p2pars, 'Pol2Params');
        end
        function test_getPol2netParameter(testCase)
            p2pars = Pol2Params();
            netparams = p2pars.returnSimRates();
            testCase.verifyEqual(length(netparams), 4);
        end
        
        %check if SimParams() exists
        function test_existSimParams(testCase)
            sp = SimParams();
            testCase.verifyClass(sp, 'SimParams');
        end
        function test_getSimParams(testCase)
            sp = SimParams();
            testCase.verifyClass(sp.k, 'double');
        end
        
%        function testSimulatePoissonSteps(testCase,numdata,rateSimplePoisson)
%            
%             N = numdata;
%             lambda = rateSimplePoisson;
%             delta = 1;
%             sigma_delta = 0.1;
%            [pwcs, p] = simPoisson(N,lambda,delta,0.0,sigma_delta);
%            fprintf('testSimulatePoissonSteps: number of elements pwcs=%i, p=%i \n',length(pwcs),length(p));
%            testCase.verifyEqual(length(pwcs),N);
%        end
    end
    
end

