% RunTest_template 
% Template file to add runs of new unit tests

% ------------------------------------------------------------------------
% Runs unit Test X - [Test Name] as part of RunTests
% ------------------------------------------------------------------------
% [Summary of test details]
%

        testnum = testnum + 1;
        testname = 'Test name and description'; % Update!
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
        config_name = '[Test config file name]'; % Update!
        %       main  % Runs calculation
        main_nonlinear % Runs calculation
        
        % Step 2 - Check results
        % run check file, script is specific to each test
        some_error_check = test_check(d);            % Update!
        if some_error_check < some_test_condtion    % Update!
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