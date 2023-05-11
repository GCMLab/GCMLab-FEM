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
        testname = '2D transient cantilever beam manufactured solution';
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
        config_name = 'TransientTest';
        main_nonlinear % Runs calculation

        
        % Step 2 - Check results
        error = Transient1D_check(d, Material, Mesh, Control,BC);
        
        fprintf('\nT3 Transient Cantilever Bar Test : Displacement error at free-end of cantilever bar is %.2f',error)
        
        convergence_tolerance = 1e-1;
        if error <= convergence_tolerance
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
