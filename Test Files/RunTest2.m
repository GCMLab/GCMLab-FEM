% SCRIPT
% Runs unit Test 2 - Patch Test B as part of RunTests
%
% Acknowledgements: Bruce Gee, Saeed Hatefi Ardakani
       


        % Add function folder to filepath    
        addpath('../Functions')

        fprintf('\n\n Test 2: Patch Test B - Q4 elements\n')
        % Step 1 - Run Simulation
        Control.config_name = 'PatchTestB';
        main
        
        % Step 2 - Check results
        [disp_er, stress_er, reaction_er] = PatchTest_check(d, stress, Fext, Mesh);
        
        fprintf('\nQ4-patch test B: Displacement error is %.2f',disp_er)
        fprintf('\nQ4-patch test B: Stress error is %.2f',stress_er)
        fprintf('\nQ4-patch test B: Reaction forces error is %.2f',reaction_er)
        
        convergence_tolerance = 1e-14;
        if disp_er <= convergence_tolerance && stress_er <= convergence_tolerance && reaction_er <= convergence_tolerance 
            test_pass = 1;
        else
            test_pass = 0;
        end
        
        
        % Step 3 - Output results
        if test_pass
            fprintf('\nPASS Patch Test B\n')
        else
            fprintf('\n\nFAIL Patch Test B\n')
            return
        end
        testpasssummary(2) = test_pass;

        
        % Step 4 - Cleanup
        clearvars -except VTKDirs ConfigFiles...
                      curDir FuncDir  ConfigDir ...
                      file codeSubmitTime ...
                      Control ntests testpasssummary