% ------------------------------------------------------------------------
% Runs unit test - Patch Test C as part of RunTests
% ------------------------------------------------------------------------
% Patch Test C is performed with node 1 fully restrained and nodes 4 and 8 
% restrained only in the x -direction. Nodal forces are applied to nodes 2,
% 3, and 6 in accordance with the values generated through the boundary 
% tractions by sigma(x)=2. The error between the FEA
% and exact solutions is then calculated. The FEA approximate solution
% should be exact.

       
        testnum = testnum + 1;
        testname = 'Patch Test C - T6 elements';
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
        quadorder = 2;
        
        config_name = 'PatchTestC';
        main
        
        % Step 2 - Check results
        [disp_er, stress_er, reaction_er] = PatchTest_check(d, stress, Fext, Mesh);
        
        fprintf('\nQ4-patch test C: Displacement error is %.2f',disp_er)
        fprintf('\nQ4-patch test C: Stress error is %.2f',stress_er)
        fprintf('\nQ4-patch test C: Reaction forces error is %.2f',reaction_er)
        
        convergence_tolerance = 1e-10;
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
        testpasssummary(testnum) = test_pass;

        
        % Step 4 - Cleanup
        clearvars -except  curDir  ConfigDir ...
                      ntests testpasssummary testnum nameslist...
                      plot2vtk VTKFolder progress_on