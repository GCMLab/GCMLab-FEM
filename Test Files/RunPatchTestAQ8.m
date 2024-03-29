% ------------------------------------------------------------------------
% Runs unit test - Patch Test A Q8 as part of RunTests
% ------------------------------------------------------------------------
% For Patch Test A, all nodes are restrained and nodal displacement values 
% are specfied according to the exact solution. The error between the FEA
% and exact solutions is then calculated. The FEA approximate solution
% should be exact.

% Note: Q8 elements use Q9.msh files and central nodes are eliminated from
% the mesh
        
        testnum = testnum + 1;
        testname = 'Patch Test A - Q8 elements with reduced integration';
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
        global meshfilename quadorder
        meshfilename = 'Mesh Files\PatchTestQ9.msh';
        quadorder = 3;

        config_name = 'PatchTestA_Q8';
%       main  % Runs calculation
        main_nonlinear % Runs calculation
        
        % Step 2 - Check results
        [disp_er, stress_er, reaction_er] = PatchTest_check(d, stress, Fext, Mesh, Material, BC);
        
        fprintf('\nQ8-patch test A: Displacement error is %.2f',disp_er)
        fprintf('\nQ8-patch test A: Stress error is %.2f',stress_er)
        fprintf('\nQ8-patch test A: Reaction forces error is %.2f',reaction_er)
        
        convergence_tolerance = 1e-10;
        if disp_er <= convergence_tolerance && stress_er <= convergence_tolerance && reaction_er <= convergence_tolerance 
            test_pass = 1;
        else
            test_pass = 0;
        end
        
        % Step 3 - Output results
        if test_pass
            fprintf('\nPASS Patch Test A\n')
        else
            fprintf('\n\nFAIL Patch Test A\n')
        end
        testpasssummary(testnum) = test_pass;

        
        % Step 4 - Cleanup
        clearvars -except  curDir  ConfigDir ...
                      ntests testpasssummary testnum nameslist...
                      plot2vtk VTKFolder progress_on