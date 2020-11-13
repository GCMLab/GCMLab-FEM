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
%         
%         % Step 4 - Cleanup
%         clearvars -except VTKDirs ConfigFiles...
%                   curDir FuncDir  ConfigDir ...
%                   file codeSubmitTime ...
%                   exit_when_done print_log ...
%                   plot2vtk
% 
%         clearvars -global
       
%% Test 1: Patch Test A - Dirichlet-Dirichlet BC
%         % Step 1 - Run Simulation
%         %Control.config_name = 'PatchTestA';
%         
%         %start_time = toc;
%         %run('Functions/Main/main');
%         %end_time = toc;
%         
%         % Step 2 - Check results
%         
%         % Step 3 - Output results
%         
%         % Step 4 - Cleanup
%         clearvars -except VTKDirs ConfigFiles...
%                   curDir FuncDir  ConfigDir ...
%                   file codeSubmitTime ...
%                   exit_when_done print_log ...
%                   plot2vtk
% 
%         clearvars -global
%         
        

%% Test 2: Patch Test B - Dirichlet-Neumann BC


%% Test 3: Manufactured Solution - Q4 elements
        fprintf('\n\n Test 3: Manufactured Solution - Q4 elements')
        % Step 1 - Run Simulation
        Control.config_name = 'ManufacturedSolution_coarse';
        
        start_time = toc;
        run('Functions/Main/main');
        end_time = toc;
        
        % Step 2 - Check results
        % run check file 
        
        % Step 3 - Output results
        if test_pass
            fprintf('\nPASS')
        else
            fprintf('\nFAIL')
        end
        
        % Step 4 - Cleanup
        clearvars -except VTKDirs ConfigFiles...
                  curDir FuncDir  ConfigDir ...
                  file codeSubmitTime ...
                  exit_when_done print_log ...
                  plot2vtk

        clearvars -global
    
    
