% ------------------------------------------------------------------------
% Runs unit Test 29 - Diffusion Test as part of RunTests
% ------------------------------------------------------------------------
% Runs a 2D steady-state isotropic heat transfer problem with fixed
% temperature on the edges and body loading
%

        testnum = testnum + 1;
        testname = 'Diffusion Test - 2D SS Heat Transfer'; 
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
        config_name = 'DiffusionTest'; % Update!
        main_nonlinear % Runs calculation
        
        % Step 2 - Check results
        % run check file, script is specific to each test
        TH_er = Diffusion_check(d, Mesh);            
        if TH_er < 0.01    % For the given mesh, error should be 0.0083 at the selected point
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