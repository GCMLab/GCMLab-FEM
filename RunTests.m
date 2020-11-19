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

% Template - copy this template to add on more unit tests

%% Test X: [Test Name] - Test Description
%         fprintf('\n\n Test X: [Test Name] - Test Description')
%         % Step 1 - Run Simulation
%         Control.config_name = '[Test config file name]';
%         
%         start_time = toc;
%         run('Functions/Main/main');
%         end_time = toc;
%         
%         % Step 2 - Check results
%         % run check file, script is specific to each test
%         
%         % Step 3 - Output results
%         if test_pass
%             fprintf('\nPASS')
%         else
%             fprintf('\nFAIL')
%         end
%         testpasssummary([testnumber]) = test_pass
%           
%
%         % Step 4 - Cleanup
%         clearvars -except VTKDirs ConfigFiles...
%                   curDir FuncDir  ConfigDir ...
%                   file codeSubmitTime ...
%                   Control ntests testpasssummary
       
%% Test 1: Patch Test A - Dirichlet-Dirichlet BC
% For Patch Test A, all nodes are restrained and nodal displacement values 
% are specfied according to the exact solution. 

        fprintf('\n\n Test 1: Patch Test A - Q4 elements\n')
        % Step 1 - Run Simulation
        global meshfilename
        Control.config_name = 'PatchTestA';
        meshfilename = 'Unstructured_sample.msh';
       
        run('Functions/Main/main');
        
        % Step 2 - Check results
        [disp_er, stress_er, reaction_er] = PatchTest_check(d, stress, Fext, Mesh);
        
        fprintf('\nQ4-patch test A: Displacement error is %.2f',disp_er)
        fprintf('\nQ4-patch test A: Stress error is %.2f',stress_er)
        fprintf('\nQ4-patch test A: Reaction forces error is %.2f',reaction_er)
        
        convergence_tolerance = 1e-14;
        if disp_er <= convergence_tolerance && stress_er <= convergence_tolerance && reaction_er <= convergence_tolerance 
            test_pass = 1;
        else
            test_pass = 0;
        end
        
        % Step 3 - Output results
        if test_pass
            fprintf('\nPASS Patch Test A\n')
        else
            fprintf('\n\nFAIL Patch Test A\n')
        end
        testpasssummary(1) = test_pass;

        
        % Step 4 - Cleanup
        clearvars -except VTKDirs ConfigFiles...
                      curDir FuncDir  ConfigDir ...
                      file codeSubmitTime ...
                      Control ntests testpasssummary
        

%% Test 2: Patch Test B - Dirichlet-Neumann BC
% For Patch Test B, only nodes 1-8 (nodes in the boundaries) are restrained
% with their displacements specified according to the exact solution. 

        fprintf('\n\n Test 2: Patch Test B - Q4 elements\n')
        % Step 1 - Run Simulation
        global meshfilename
        Control.config_name = 'PatchTestB';
        meshfilename = 'Unstructured_sample.msh';
       
        run('Functions/Main/main');
        
        % Step 2 - Check results
        [disp_er, stress_er, reaction_er] = PatchTest_check(d, stress, Fext, Mesh);
        
        fprintf('\nQ4-patch test B: Displacement error is %.2f',disp_er)
        fprintf('\nQ4-patch test B: Stress error is %.2f',stress_er)
        fprintf('\nQ4-patch test B: Reaction forces error is %.2f',reaction_er)
        
        convergence_tolerance = 1e-14;
        if disp_er <= convergence_tolerance && stress_er <= convergence_tolerance && reaction_er <= convergence_tolerance 
            test_pass = 1;
        else
            test_pass = 0;
        end
        
        
        % Step 3 - Output results
        if test_pass
            fprintf('\nPASS Patch Test B\n')
        else
            fprintf('\n\nFAIL Patch Test B\n')
            return
        end
        testpasssummary(2) = test_pass;

        
        % Step 4 - Cleanup
        clearvars -except VTKDirs ConfigFiles...
                      curDir FuncDir  ConfigDir ...
                      file codeSubmitTime ...
                      Control ntests testpasssummary

                  
%% Test 3: Patch Test C
% Patch Test C is performed with node 1 fully restrained and nodes 4 and 8 
% restrained only in the x -direction. Nodal forces are applied to nodes 2,
% 3, and 6 in accordance with the values generated through the boundary 
% tractions by sigma(x)=2 

        fprintf('\n\n Test 3: Patch Test C - Q4 elements\n')
        % Step 1 - Run Simulation
        global meshfilename
        Control.config_name = 'PatchTestC';
        meshfilename = 'Unstructured_sample.msh';
       
        run('Functions/Main/main');
        
        % Step 2 - Check results
        [disp_er, stress_er, reaction_er] = PatchTest_check(d, stress, Fext, Mesh);
        
        fprintf('\nQ4-patch test C: Displacement error is %.2f',disp_er)
        fprintf('\nQ4-patch test C: Stress error is %.2f',stress_er)
        fprintf('\nQ4-patch test C: Reaction forces error is %.2f',reaction_er)
        
        convergence_tolerance = 1e-14;
        if disp_er <= convergence_tolerance && stress_er <= convergence_tolerance && reaction_er <= convergence_tolerance 
            test_pass = 1;
        else
            test_pass = 0;
        end
        
        % Step 3 - Output results
        if test_pass
            fprintf('\nPASS Patch Test C\n')
        else
            fprintf('\n\nFAIL Patch Test C\n')
            return
        end
        testpasssummary(3) = test_pass;

        
        % Step 4 - Cleanup
        clearvars -except VTKDirs ConfigFiles...
                      curDir FuncDir  ConfigDir ...
                      file codeSubmitTime ...
                      Control ntests testpasssummary



%% Test 4: Manufactured Solution - Q4 element convergence
         run('Test Files/RunTest4')
    
    
%% Test 5: Manufactured Solution - Q9 element convergence
        run('Test Files/RunTest5')
        

    
%% Add new tests here

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