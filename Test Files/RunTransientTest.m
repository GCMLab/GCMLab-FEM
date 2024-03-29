% ------------------------------------------------------------------------
% Runs unit test - Transient Solution as part of RunTests
% ------------------------------------------------------------------------
% Test consists of a 2D cantilever bar simulation where the left is fixed.
% Boundary conditions at the fixed-end are restricted to displacements in
% the x-direction. The bottom node at the fixed-end restricts displacement
% in both the x and y directions to avoid 2d effects. An external force is
% applied to all nodes at the free end of the bar. The solution is compared
% to the 1D analytical solution whereby the displacement at node 863 
% (i.e. midpoint in y-direction at free-end) is compared over time. A 
% slight variation in the numerical and analytical solution is expected due 
% to the truncation error and the comparison to an exact one-dimensional 
% solution. 

       
        testnum = testnum + 1;
        testname = '2D transient cantilever beam solution';
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
        L2d = Transient1D_check(d, Material, Mesh, Control, BC);
        
        fprintf('\nT3 Transient Cantilever Bar Test: L2-norm of the displacement error at free-end of cantilever bar is %.2e', L2d);

        convergence_tolerance = 1e-2;
        if L2d <= convergence_tolerance
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
