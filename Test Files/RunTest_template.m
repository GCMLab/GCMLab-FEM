% RunTest_template 
% Template file to add runs of new unit tests

% ------------------------------------------------------------------------
% Runs unit Test X - [Test Name] as part of RunTests
% ------------------------------------------------------------------------
% [Summary of test details]
%
% Acknowledgements: Bruce Gee, X
       
        % Create test VTK folder
        if plot2vtk
            vtk_dir = fullfile(VTKFolder,'\TestX');
            if ~isfolder(vtk_dir) 
                mkdir(vtk_dir)
            end
        end

        fprintf('\n\n Test X: [Test Name] - Test Description')
        % Step 1 - Run Simulation
        config_name = '[Test config file name]';
        main  % Runs calculation
        
        % Step 2 - Check results
        % run check file, script is specific to each test
        some_error_check = test_check();
        if some_error_check < some_test_condtion
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
        testpasssummary([testnumber]) = test_pass;
          

        % Step 4 - Cleanup
        clearvars -except  curDir  ConfigDir ...
                      ntests testpasssummary...
                      plot2vtk VTKFolder progress_on