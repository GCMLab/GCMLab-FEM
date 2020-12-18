% ------------------------------------------------------------------------
% Runs unit Test 6 - Plate with Hole as part of RunTests
% ------------------------------------------------------------------------
% Runs plate case of Plate with hole under tension. Adapted from OpenFOAM
% tutorial problem: https://www.openfoam.com/documentation/tutorial-guide/tutorialse9.php
% TODO: Add description of test case to Wiki



        % Create test VTK folder
        if plot2vtk
            vtk_dir = fullfile(VTKFolder,'\Test6');
            if ~isfolder(vtk_dir) 
                mkdir(vtk_dir)
            end
        end

        fprintf('\n\n Test 6: Plate with Hole')
        % Step 1 - Run Simulation
        global calc_type
            config_name = 'PlateWithHole';
            % run with nodal averaging
            calc_type = 'nodal';
            main
           
            stress_nodal = stress;

            % run with L2 projection
            calc_type = 'L2projection';
            main

            stress_L2 = stress;
        
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