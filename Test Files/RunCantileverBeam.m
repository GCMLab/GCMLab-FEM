% ------------------------------------------------------------------------
% Runs unit Test 7 - Cantilever Beam as part of RunTests
% ------------------------------------------------------------------------
% Runs cantilever beam under bending. The nodes located at the left edge is
% fully restrained, and nodal forces are applied to the nodes located at
% the right edge by sigma(y) = -4e5. Then, in order to check shear locking,
% the error between the FEA and theoretical solutions is then calculated.
% The FEA approximate solution should be close to theoretical solution.

        testnum = testnum + 1;
        testname = 'Cantilever Beam - Q4 elements';
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
        global  E nu t
        t = -4e5; % applied traction
        E = 2e11;  % elastic modulus
        nu = 0.3;  % poisson's ratio
        
        config_name = 'CantileverBeam';
        main_static
        
        % Step 2 - Check results
        % run check file, script is specific to each test
        [disp_er] = CantileverBeam_check(d, Material, BC, Mesh);
        tolerance_er = 1e-3;
        if disp_er < tolerance_er
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