% This script runs all units tests sequentially and checks the results
%   Acknowledgements: Bruce Gee

%% Clear variables and initialize code
    clearvars -global
    clear, clc, close all
    format compact
    tic;
    
    % Test VTK output
    plot2vtk = 0;
%       VTKFolder ='C:\Users\b3gee\Documents\Matlab Results\';
%       VTKFolder = 'C:\Users\shatefia\OneDrive - University of Waterloo\Documents\Fall 2022\Github\Documents\Matlab Results\';
%       VTKFolder = 'C:\Users\shatefia\OneDrive - University of Waterloo\Documents\UWaterloo\Research\GitHub\Results';
%       VTKFolder = 'C:\Users\knbetanc\OneDrive - University of Waterloo\Documents\UWaterloo\Research\GitHub\Results\';
       VTKFolder ='C:\Users\GZingaro\Desktop\Results';
%       VTKFolder = 'C:\Users\bccampos\Downloads\Matlab Results\';
%       VTKFolder = 'C:\Users\bu_ca\Downloads\Matlab Results\';

    % suppress progress messages
    progress_on = 0;
    
    % Current directory
    curDir = pwd;
    DirFolder = 'Test Files';

%% Directories
    FuncDir = fullfile(curDir, 'Functions');
    ConfigDir = fullfile(curDir, DirFolder);
    VTKFolder = fullfile(VTKFolder, DirFolder);
    % add paths
    addpath(genpath(FuncDir));
    addpath(genpath(ConfigDir));
       
    % number of tests - Update when new tests added!
    ntests = 32; 

    nameslist = {};
    testnum = 0;

    % initialize test summary
    testpasssummary = zeros(ntests,1);
    
%% Run Tests

%% Test 1: Patch Test A - All nodal displacements prescribed Q4
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestAQ4')

%% Test 2: Patch Test B - Dirichlet-Dirichlet BC Q4
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestBQ4')

%% Test 3: Patch Test C - Dirichlet-Neumann BC Q4
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestCQ4')

%% Test 4: Manufactured Solution - Q4 element convergence
%   Pass condition: L2-norm converges at a rate of at least h^2
%                    e-norm converges at a rate of at least h
         run('Test Files\RunManSolQ4')

%% Test 5: Patch Test A - All nodal displacements prescribed Q9
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestAQ9')

%% Test 6: Patch Test B - Dirichlet-Dirichlet BC Q9
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestBQ9')

%% Test 7: Patch Test C - Dirichlet-Neumann BC Q9
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestCQ9')

%% Test 8: Manufactured Solution - Q9 element convergence
%   Pass condition: L2-norm converges at a rate of at least h^3
%                    e-norm converges at a rate of at least h^2
        run('Test Files\RunManSolQ9')

%% Test 9: Patch Test A - All nodal displacements prescribed T3
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestAT3')

%% Test 10: Patch Test B - Dirichlet-Dirichlet BC T3
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestBT3')

%% Test 11: Patch Test C - Dirichlet-Neumann BC T3
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestCT3')

%% Test 12: Manufactured Solution - T3 element convergence
%   Pass condition: L2-norm converges at a rate of at least h^2
%                    e-norm converges at a rate of at least h
         run('Test Files\RunManSolT3')

%% Test 13: Patch Test A - All nodal displacements prescribed T6
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestAT6')

%% Test 14: Patch Test B - Dirichlet-Dirichlet BC T6
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestBT6')

%% Test 15: Patch Test C - Dirichlet-Neumann BC T6
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestCT6')

%% Test 16: Manufactured Solution - T6 element convergence
%   Pass condition: L2-norm converges at a rate of at least h^2
%                    e-norm converges at a rate of at least h
         run('Test Files\RunManSolT6')

%% Test 17: Manufactured Solution - Plane Strain check
%   Pass condition: L2-norm converges at a rate of at least h^2
%                    e-norm converges at a rate of at least h
        run('Test Files\RunPstrain')


%% Test 18: Plate with Hole under tension
%  Pass Condtion: error of L2 projected stresses is less than nodal
%  averaged stresses
      run('Test Files\RunPlatewHole')

%% Test 19: Cantilever Beam
%  Pass Condition: Shear locking is prevented and error in displacement is
%  acceptable
      run('Test Files\RunCantileverBeam')

%% Test 20: Patch Test A - All nodal displacements prescribed Q8
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestAQ8')

%% Test 21: Patch Test B - Dirichlet-Dirichlet BC Q8
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestBQ8')

%% Test 22: Patch Test C - Dirichlet-Neumann BC Q8
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files\RunPatchTestCQ8')

%% Test 23: Manufactured Solution - Q8 element convergence
% Pass Condition: FEA solution displacements, stresses, and strains are exact
        run('Test Files/RunManSolQ8')

%% Test 24: One Dimensional Problem
%  Pass Condition: FEA solution displacements are exact
      run('Test Files\RunTest1D')

%% Test 25: Circular Inclusion Problem
%  Pass Condition: FEA solution stresses are close to analytical solution
      run('Test Files\RunCircularInclusion')

%% Test 26: Time Dependent Dirichlet - Time dependent dirichlet boundary conditions
%  Pass Condtion: Shear strain over time matches manufactured solution
      run('Test Files\RunDirichletTime')

%% Test 27: Plate with Hole under tension using Hypermesh input
%  Pass Condtion: error of L2 projected stresses is less than nodal
%  averaged stresses
      run('Test Files/RunPlatewHole_hypermesh')   

%% Test 28: Transient Test - 2D elastic bar with Time-Dependent Neumann B.C at free end 
%  Pass Condtion: Displacement matches analytical solution over time
      run('Test Files/RunTransientTest')

%% Test 29: Diffusion Test - Steady State Isotropic Heat Diffusion Problem
%  Pass Condtion: Temperature at reference point matches analytical  solution
      run('Test Files/RunDiffusionTest')
      

%% Test 30: Kelvin-Voigt Linear Viscoelastic Model - Relaxation
%  Pass Condtion:
     run('Test Files/RunTestVE1')     
     
 %% Test 31 Dynamic Test - Manufacture solution
%  Pass Condtion: Discrete norm of displacemets vs exact solution and
%                 stored solution from previous run
      run('Test Files\RunDynamicManufacture')
      
%% Test 32: Thermoelastic Test - Manufactured solution
%  Pass Condtion: FEA solution displacements, stresses, and strains are exact
     run('Test Files\RunThermoelasticTest')

%% Test X: [Test Name] - Short Test Description
%  Pass Condtion:
%      run('Test Files/RunTestX')

%% Summarize test results
fprintf('\n\n%-10s%-10s%-20s', 'Test', 'Status','Test Name')
fprintf('\n----------------------------------------------')
for test = 1:ntests
   if testpasssummary(test)
       fprintf('\n%-10d%-10s%-20s', test, 'PASS',nameslist{test})
   else
      fprintf('\n%-10d%-10s%-20s', test, '! - FAIL',nameslist{test})
   end
end
fprintf('\n\n')
