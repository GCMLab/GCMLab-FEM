% This script runs all units tests sequentially and checks the results
%   Acknowledgements: Bruce Gee, Saeed Hatefi Ardakani



%% Clear variables and initialize code
    clearvars -global
    clear, clc, close all
    format compact
    tic;
    
    % Do not output any vtk files from unit tests
    Control.vtk = 0;
    Control.vtk_dir = [];
    
    % Current time and date
    codeSubmitTime = datetime('now','Format','yyyy-MM-dd''_''HH-mm-ss');

    % Current directory
    curDir = pwd;
    DirFolder = 'Test Files';
    ConfigDir = fullfile(curDir, DirFolder);
    %% Directories
    Control.config_dir = ConfigDir;
    FuncDir = fullfile(curDir, 'Functions');
    ConfigDir = fullfile(curDir, DirFolder);
    % add paths
    addpath(genpath(FuncDir));
    addpath(genpath(ConfigDir));
    
    % number of tests
    ntests = 5;
    % initialize test summary
    testpasssummary = zeros(ntests,1);
    
%% Run Tests

%% Test 1: Patch Test A - Dirichlet-Dirichlet BC
% For Patch Test A, all nodes are restrained and nodal displacement values 
% are specfied according to the exact solution.
%
% Pass Condition: FEA solution displacements, stresses, and strains are
% exact
        run('Test Files/RunTest1')

%% Test 2: Patch Test B - Dirichlet-Neumann BC
% For Patch Test B, only nodes 1-8 (nodes in the boundaries) are restrained
% with their displacements specified according to the exact solution. 
%
% Pass Condition: FEA solution displacements, stresses, and strains are
% exact
        run('Test Files/RunTest2')

%% Test 3: Patch Test C
% Patch Test C is performed with node 1 fully restrained and nodes 4 and 8 
% restrained only in the x -direction. Nodal forces are applied to nodes 2,
% 3, and 6 in accordance with the values generated through the boundary 
% tractions by sigma(x)=2 
%
% Pass Condition: FEA solution displacements, stresses, and strains are
% exact
        run('Test Files/RunTest3')

%% Test 4: Manufactured Solution - Q4 element convergence
%   Pass condition: L2-norm converges at a rate of at least h^2
%                    e-norm converges at a rate of at least h
%
% Pass Condition: FEA solution displacements, stresses, and strains are
% exact
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