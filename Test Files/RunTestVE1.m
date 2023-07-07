% ------------------------------------------------------------------------
% Runs unit test - Kelvin-Voigt Linear Viscoelastic Model as part of RunTests
% ------------------------------------------------------------------------
% Test consists of a 2D cantilever bar simulation where the left is fixed.
% Boundary conditions at the fixed-end are restricted to displacements in
% the x-direction. The bottom node at the fixed-end restricts displacement
% in both the x and y directions to avoid 2d effects. An external force is
% applied to all nodes at the free end of the bar. The solution is compared
% to the 1D analytical solution whereby the displacement at node 863 
% (i.e. midpoint in y-direction at free-end) is compared over time. 
       
        testnum = testnum + 1;
        testname = '2D Kelvin-Voigt Cantilever Beam Solution';
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
        config_name = 'PatchTestVE1';
        main_nonlinear % Runs calculation

        
        % Step 2 - Check results
        error = VE1_1D_check(d, Material, Mesh, Control,BC);
        
        fprintf('\nT3 Kelvin-Voigt Cantilever Bar Test: Displacement error at the free end of the cantilever bar is %.3e', error);
        
        convergence_tolerance = 1e-1;
        if error <= convergence_tolerance
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
