%%
%clear all;
addpath('../');

suite_getParams = matlab.unittest.TestSuite.fromFile('Utest_getParameters.m');
{suite_getParams.Name}'
getParams_out = {suite_getParams.Name}';

for i=1:length(getParams_out)
    fprintf('test = %s \n',getParams_out{i});
end
suite_getParams.run