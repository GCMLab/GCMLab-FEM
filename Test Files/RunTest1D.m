% ------------------------------------------------------------------------
% Runs unit Test 8 - One Dimensional problem as part of RunTests
% ------------------------------------------------------------------------
% Runs One Dimensional problem under tension. The node located at the left 
% edge is fully restrained, and nodal force is applied to the node located 
% at the right edge by F = 6e5. Furthermore, a distributed body force is 
% applied to the problem by b = 2e5. Then, the error between the FEA and 
% analytical solution is calculated. The FEA approximate solution should be
% analytical one.
% TODO: Add description of test case to Wiki

        testnum = testnum + 1;
        testname = 'One Dimensional Problem - 1D elements';
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
        global  E nu t b
        t = 6e9; % applied traction [N]
        b = 2e9; % applied body force [N/m]
        E = 2e11;  % elastic modulus [Pa]
        nu = 0;  % poisson's ratio
        
        config_name = 'Test1D';
        main_static
        
        % Step 2 - Check results
        % run check file, script is specific to each test
        [disp_er] = Test1D_check(d, Material, BC, Mesh);
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