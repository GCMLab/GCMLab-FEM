% ------------------------------------------------------------------------
% Runs unit Test 26 - Dirichlet Time
% ------------------------------------------------------------------------
% This test applies time dependent Dirichlet boundary conditions on a
% square domain and computes the stress as a function of time.
%

        testnum = testnum + 1;
        testname = 'Time dependent manufactured solution'; 
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
        global Omega1 Omega2 E nu
        Omega1 = 2;
        Omega2 = 3;
        E = 2.5e11;
        nu = 0.25;
        config_name = 'ManufacturedSolution_DirichletTime'; 
        main  % Runs calculation
        
        % Step 2 - Check results
        % run check file, script is specific to each test
        time_er = DirichletTime_check(sSave);          
        if time_er < 1e-10
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