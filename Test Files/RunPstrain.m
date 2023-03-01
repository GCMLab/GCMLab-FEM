% ------------------------------------------------------------------------
% Runs unit test - Q4 Plane strain manufactured solution convergence as a part of RunTests
% ------------------------------------------------------------------------
% Calculates the convergence rates of a uniform Q4 mesh using a
% manufactured solution in which 
% ux = x^5 + x*y^3 - y^6
% uy = x^5 + x*y^3 - y^6
% under plane strain conditions

        testnum = testnum + 1;
        testname = 'Plane Strain manufactured solution';
        nameslist{testnum} = testname;
       
        % Create test VTK folder
        if plot2vtk
            folname = ['\Test',num2str(testnum)];
            vtk_dir = fullfile(VTKFolder,folname);
            if ~isfolder(vtk_dir) 
                mkdir(vtk_dir)
            end
        end
        % test runs 3 meshes, only finest mesh will be saved
    
        fprintf('\n\n Test %d : %s\n', testnum, testname)
        % Step 1 - Run Simulation
        global meshfilename quadorder E nu
            E = 2230;
            nu = 0.3;
            quadorder = 2; 
            config_name = 'ManufacturedSolution_PlaneStrain';
            
            % Run coarse mesh
            meshfilename = 'Mesh Files\Manufactured_coarseQ4.msh';
%           main  % Runs calculation
            main_nonlinear % Runs calculation
            
            d_coarse = d;
            stress_coarse = stress;
            strain_coarse = strain;
            Mesh_coarse = Mesh;

            % Run fine mesh
            meshfilename = 'Mesh Files\Manufactured_fineQ4.msh';
            main

            d_fine = d;
            stress_fine = stress;
            strain_fine = strain;
            Mesh_fine = Mesh;

            % Run finer mesh
            meshfilename = 'Mesh Files\Manufactured_finerQ4.msh';
            main

            d_finer = d;
            stress_finer = stress;
            strain_finer = strain;
            Mesh_finer = Mesh;
        
        % Step 2 - Check results
            [m_L2, m_e] = ManufacturedSolution_PlaneStrain_check(d_coarse, d_fine, d_finer, ...
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
            end
        testpasssummary(testnum) = test_pass;

        
        % Step 4 - Cleanup
        clearvars -except  curDir  ConfigDir ...
                      ntests testpasssummary testnum nameslist...
                      plot2vtk VTKFolder progress_on