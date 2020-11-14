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
%                   Control
       
%% Test 1: Patch Test A - Dirichlet-Dirichlet BC

        

%% Test 2: Patch Test B - Dirichlet-Neumann BC


%% Test 3: Manufactured Solution - Q4 elements
        fprintf('\n\n Test 3: Manufactured Solution - Q4 elements\n')
        % Step 1 - Run Simulation
        global meshfilename quadorder 
            quadorder = 2; 
            Control.config_name = 'ManufacturedSolution';
            % Run coarse mesh
            meshfilename = 'Manufactured_coarseQ4.msh';
            run('Functions/Main/main');

            d_coarse = d;
            stress_coarse = stress;
            strain_coarse = strain;
            Mesh_coarse = Mesh;


            % Run fine mesh
            meshfilename = 'Manufactured_fineQ4.msh';
            run('Functions/Main/main');

            d_fine = d;
            stress_fine = stress;
            strain_fine = strain;
            Mesh_fine = Mesh;

            % Run finer mesh
            meshfilename = 'Manufactured_finerQ4.msh';
            run('Functions/Main/main');

            d_finer = d;
            stress_finer = stress;
            strain_finer = strain;
            Mesh_finer = Mesh;
        
        % Step 2 - Check results
            [m_L2, m_e] = ManufacturedSolution_check(d_coarse, d_fine, d_finer, ...
                stress_coarse, stress_fine, stress_finer, strain_coarse, strain_fine, ...
                strain_finer, Mesh_coarse, Mesh_fine, Mesh_finer);
            
            fprintf('\nQ4 L2-norm converges at a rate of %.2f',m_L2)
            fprintf('\nQ4  e-norm converges at a rate of %.2f',m_e)
            
            convergence_tolerance = 0.05;
            if m_L2 >= (2 - convergence_tolerance) && m_e >= (1 - convergence_tolerance)
                test_pass = 1;
            else
                test_pass = 0;
            end
        
        
        % Step 3 - Output results
            if test_pass
                fprintf('\nPASS\n')
            else
                fprintf('\n\nFAIL\n')
                return
            end
        
        % Step 4 - Cleanup
            clearvars -except VTKDirs ConfigFiles...
                      curDir FuncDir  ConfigDir ...
                      file codeSubmitTime ...
                      Control
    
    
%% Test 4: Manufactured Solution - Q9 elements
        return
        % Unit testing for convergence of Q9 elements is setup, but the
        % code is not currently capable of handling Q9 elements
        fprintf('\n\n Test 3: Manufactured Solution - Q9 elements\n')
        % Step 1 - Run Simulation
        global meshfilename quadorder 
            quadorder = 4;
            Control.config_name = 'ManufacturedSolution';
            % Run coarse mesh
            meshfilename = 'Manufactured_coarseQ9.msh';
            run('Functions/Main/main');

            d_coarse = d;
            stress_coarse = stress;
            strain_coarse = strain;
            Mesh_coarse = Mesh;


            % Run fine mesh
            meshfilename = 'Manufactured_fineQ9.msh';
            run('Functions/Main/main');

            d_fine = d;
            stress_fine = stress;
            strain_fine = strain;
            Mesh_fine = Mesh;

            % Run finer mesh
            meshfilename = 'Manufactured_finerQ9.msh';
            run('Functions/Main/main');

            d_finer = d;
            stress_finer = stress;
            strain_finer = strain;
            Mesh_finer = Mesh;
        
        % Step 2 - Check results
            [m_L2, m_e] = ManufacturedSolution_check(d_coarse, d_fine, d_finer, ...
                stress_coarse, stress_fine, stress_finer, strain_coarse, strain_fine, ...
                strain_finer, Mesh_coarse, Mesh_fine, Mesh_finer);
            
            fprintf('\nQ4 L2-norm converges at a rate of %.2f',m_L2)
            fprintf('\nQ4  e-norm converges at a rate of %.2f',m_e)
            
            convergence_tolerance = 0.05;
            if m_L2 >= (3 - convergence_tolerance) && m_e >= (2 - convergence_tolerance)
                test_pass = 1;
            else
                test_pass = 0;
            end
        
        
        % Step 3 - Output results
            if test_pass
                fprintf('\nPASS\n')
            else
                fprintf('\n\nFAIL\n')
                return
            end
        
        % Step 4 - Cleanup
            clearvars -except VTKDirs ConfigFiles...
                      curDir FuncDir  ConfigDir ...
                      file codeSubmitTime ...
                      Control
    
