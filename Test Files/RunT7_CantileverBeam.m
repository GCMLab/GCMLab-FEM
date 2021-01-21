% ------------------------------------------------------------------------
% Runs unit Test 7 - Cantilever Beam as part of RunTests
% ------------------------------------------------------------------------
% Runs cantilever beam under bending. The nodes located at the left edge is
% fully restrained, and nodal forces are applied to the nodes located at
% the right edge by sigma(y) = -4e5. Then, in order to check shear locking,
% the error between the FEA and theoretical solutions is then calculated.
% The FEA approximate solution should be close to theoretical solution.
% TODO: Add description of test case to Wiki



        % Create test VTK folder
        if plot2vtk
            vtk_dir = fullfile(VTKFolder,'\Test7');
            if ~isfolder(vtk_dir) 
                mkdir(vtk_dir)
            end
        end

        fprintf('\n\n Test 7: Cantilever Beam')
        % Step 1 - Run Simulation
        global  E nu t
        t = -4e5; % applied traction
        E = 2e11;  % elastic modulus
        nu = 0.3;  % poisson's ratio
        
        config_name = 'CantileverBeam';
        main
        
        % Step 2 - Check results
        % run check file, script is specific to each test
        [error_nodal, error_L2] = PlateWithHole_check(Mesh,stress_nodal,stress_L2);
        if error_L2 < error_nodal
            test_pass = 1;
        else
            test_pass = 0;
        end
        
        
        % Step 3 - Output results
        if test_pass
            fprintf('\nPASS')
        else
            fprintf('\nFAIL')
        end
        testpasssummary(6) = test_pass;
          

        % Step 4 - Cleanup
        clearvars -except  curDir  ConfigDir ...
                      ntests testpasssummary...
                      plot2vtk VTKFolder progress_on