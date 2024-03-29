% ------------------------------------------------------------------------
% Runs unit Test - Cantilever Beam as part of RunTests
% ------------------------------------------------------------------------
% Runs cantilever beam under bending. The nodes located at the left edge is
% fully restrained, and nodal forces are applied to the nodes located at
% the right edge by sigma(y) = -4e5. Then, in order to check shear locking,
% the error between the FEA and theoretical solutions is then calculated.
% The FEA approximate solution should be close to theoretical solution.

        testnum = testnum + 1;
        testname = 'Cantilever Beam - Q4 elements, sinusoidal time-dependent force';
        nameslist{testnum} = testname;


        % Create test VTK folder
        if plot2vtk
            folname = ['\Test',num2str(testnum)];
            vtk_dir = fullfile(VTKFolder,folname);
            if ~isfolder(vtk_dir) 
                mkdir(vtk_dir)
            end
        end

        fprintf('\n\n Test %d : %s\n', testnum, testname)
        % Step 1 - Run Simulation
        config_name = 'CantileverBeam';
%         main  % Runs calculation
        main_nonlinear % Runs calculation
        
        % Step 2 - Check results
        % run check file, script is specific to each test
        [disp_er, time_er] = CantileverBeam_check(d, Material, BC, Mesh);
        tolerance_er = 1e-3;
        if disp_er < tolerance_er && time_er < tolerance_er
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
        testpasssummary(testnum) = test_pass;
          

        % Step 4 - Cleanup
        clearvars -except  curDir  ConfigDir ...
                      ntests testpasssummary testnum nameslist...
                      plot2vtk VTKFolder progress_on