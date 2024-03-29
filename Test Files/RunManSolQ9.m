% ------------------------------------------------------------------------
% Runs unit test - Q9 manufactured solution convergence as a part of RunTests
% ------------------------------------------------------------------------
% Calculates the convergence rates of a uniform Q9 mesh using a
% manufactured solution in which 
% ux = x^5 + x*y^3 - y^6
% uy = x^5 + x*y^3 - y^6
% under plane stress conditions


        testnum = testnum + 1;
        testname = 'Q9 Convergence - Plane stress manufactured solution';
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
        global meshfilename quadorder
            quadorder = 4;
            config_name = 'ManufacturedSolution_PlaneStress';
            
            % Run coarse mesh
            meshfilename = 'Mesh Files\Manufactured_coarseQ9.msh';
%           main  % Runs calculation
            main_nonlinear % Runs calculation

            d_coarse = d;
            stress_coarse = stress;
            strain_coarse = strain;
            Mesh_coarse = Mesh;


            % Run fine mesh
            meshfilename = 'Mesh Files\Manufactured_fineQ9.msh';
%           main  % Runs calculation
            main_nonlinear % Runs calculation

            d_fine = d;
            stress_fine = stress;
            strain_fine = strain;
            Mesh_fine = Mesh;
            

            % Run finer mesh
            meshfilename = 'Mesh Files\Manufactured_finerQ9.msh';
%           main  % Runs calculation
            main_nonlinear % Runs calculation
            
            d_finer = d;
            stress_finer = stress;
            strain_finer = strain;
            Mesh_finer = Mesh;
        
        % Step 2 - Check results
            [m_L2, m_e] = ManufacturedSolution_PlaneStress_check(d_coarse, d_fine, d_finer, ...
                stress_coarse, stress_fine, stress_finer, strain_coarse, strain_fine, ...
                strain_finer, Mesh_coarse, Mesh_fine, Mesh_finer, Material, Control);
            
            fprintf('\nQ9 L2-norm converges at a rate of %.2f',m_L2)
            fprintf('\nQ9  e-norm converges at a rate of %.2f',m_e)
            
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
            end
        testpasssummary(testnum) = test_pass;

        
        % Step 4 - Cleanup
        clearvars -except  curDir  ConfigDir ...
                      ntests testpasssummary testnum nameslist...
                      plot2vtk VTKFolder progress_on