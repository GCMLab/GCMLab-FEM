% SCRIPT
% Runs unit Test 3 - Patch Test C as part of RunTests
%
% Acknowledgements: Bruce Gee, Saeed Hatefi Ardakani
       


        % Add function folder to filepath    
        addpath('../Functions')

        fprintf('\n\n Test 3: Patch Test C - Q4 elements\n')
        % Step 1 - Run Simulation
        Control.config_name = 'PatchTestC';
        main
        
        % Step 2 - Check results
        [disp_er, stress_er, reaction_er] = PatchTest_check(d, stress, Fext, Mesh);
        
        fprintf('\nQ4-patch test C: Displacement error is %.2f',disp_er)
        fprintf('\nQ4-patch test C: Stress error is %.2f',stress_er)
        fprintf('\nQ4-patch test C: Reaction forces error is %.2f',reaction_er)
        
        convergence_tolerance = 1e-14;
        if disp_er <= convergence_tolerance && stress_er <= convergence_tolerance && reaction_er <= convergence_tolerance 
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
        testpasssummary(3) = test_pass;

        
        % Step 4 - Cleanup
        clearvars -except VTKDirs ConfigFiles...
                      curDir FuncDir  ConfigDir ...
                      file codeSubmitTime ...
                      Control ntests testpasssummary