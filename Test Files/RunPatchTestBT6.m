% ------------------------------------------------------------------------
% Runs unit test - Patch Test B T6 as part of RunTests
% ------------------------------------------------------------------------
% For Patch Test B, only nodes 1-8 (nodes in the boundaries) are restrained
% with their displacements specified according to the exact solution. The error between the FEA
% and exact solutions is then calculated. The FEA approximate solution
% should be exact.

        testnum = testnum + 1;
        testname = 'Patch Test B - T6 elements';
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
        global  E nu t meshfilename quadorder
        t = 3.495; % applied traction (both directions)
        E = 2540;  % elastic modulus
        nu = 0.3;  % poisson's ratio
        
        meshfilename = 'Mesh Files\PatchTestT6.msh';
        quadorder = 0;
        
        config_name = 'PatchTestB';
        main
        
        % Step 2 - Check results
        [disp_er, stress_er, reaction_er] = PatchTest_check(d, stress, Fext, Mesh);
        
        fprintf('\nQ4-patch test B: Displacement error is %.2f',disp_er)
        fprintf('\nQ4-patch test B: Stress error is %.2f',stress_er)
        fprintf('\nQ4-patch test B: Reaction forces error is %.2f',reaction_er)
        
        convergence_tolerance = 1e-10;
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
        testpasssummary(testnum) = test_pass;

        
        % Step 4 - Cleanup
        clearvars -except  curDir  ConfigDir ...
                      ntests testpasssummary testnum nameslist...
                      plot2vtk VTKFolder progress_on