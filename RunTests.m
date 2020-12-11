% This script runs all units tests sequentially and checks the results
%   Acknowledgements: Bruce Gee

%% Clear variables and initialize code
    clearvars -global
    clear, clc, close all
    format compact
    tic;
    
    % Test VTK output
    plot2vtk = 0;
    VTKFolder ='C:\Users\b3gee\Documents\Matlab Results\';
    % suppress progress messages
    progress_on = 0;
    
    % Current directory
    curDir = pwd;
    DirFolder = 'Test Files';
    ConfigDir = fullfile(curDir, DirFolder);

%% Directories
    FuncDir = fullfile(curDir, 'Functions');
    ConfigDir = fullfile(curDir, DirFolder);
    VTKFolder = fullfile(VTKFolder, DirFolder);
    % add paths
    addpath(genpath(FuncDir));
    addpath(genpath(ConfigDir));
    
    % number of tests
    ntests = 5;
    % initialize test summary
    testpasssummary = zeros(ntests,1);
    
%% Run Tests

%% Test 1: Patch Test A - All nodal displacements prescribed
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files/RunTest1')

%% Test 2: Patch Test B - Dirichlet-Dirichlet BC
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files/RunTest2')

%% Test 3: Patch Test C - Dirichlet-Neumann BC
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files/RunTest3')

%% Test 4: Manufactured Solution - Q4 element convergence
%   Pass condition: L2-norm converges at a rate of at least h^2
%                    e-norm converges at a rate of at least h
         run('Test Files/RunTest4')
 
%% Test 5: Manufactured Solution - Q9 element convergence
%   Pass condition: L2-norm converges at a rate of at least h^3
%                    e-norm converges at a rate of at least h^2
        run('Test Files/RunTest5')
        
%% Test X: [Test Name] - Short Test Description
%   Pass Condtion:
%       run('Test Files/RunTestX')

%% Summarize test results
fprintf('\n\n%10s%10s', 'Test', 'Status')
fprintf('\n----------------------------------------------')
for test = 1:ntests
   if testpasssummary(test)
       fprintf('\n%10d%10s', test, 'PASS')
   else
      fprintf('\n%10d%10s', test, 'FAIL')
   end
end
fprintf('\n\n')