% ------------------------------------------------------------------------
% Runs unit test - Thermoelastic Solution as part of RunTests
% ------------------------------------------------------------------------
% Test consists of a 2D plate where all the boundaries are fixed according
% to a prescribed manufactured solution. The error between the finite
% element solution and such manufactured solution is then calculated. The
% FE solution should be exact.

       
        testnum = testnum + 1;
        testname = '2D quasi-steady thermoelastic solution';
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
        config_name = 'ThermoElastic_Dirichlet';
        main_nonlinear % Runs calculation

        
        % Step 2 - Check results
        [d_er, reaction_er] = Thermoelasticity_check(d, Fext, Mesh);
        
        fprintf('\nT3 Quasi-steady Thermoelastic Test: Disp./Temp. error is %.2e', d_er);
        fprintf('\nT3 Quasi-steady Thermoelastic Test: Reaction forces error is %.2e', reaction_er);
        
        convergence_tolerance = 1e-10;
        if d_er <= convergence_tolerance && reaction_er <= convergence_tolerance
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
