% ------------------------------------------------------------------------
% Runs unit Test 21 - Circular Inclusion as part of RunTests
% ------------------------------------------------------------------------
% Runs plate case of circular plate composed of two different materials.
% Adapted from the reference Sukumar, Natarajan, David L. Chopp, 
% Nicolas Moës, and Ted Belytschko. "Modeling holes and inclusions by level
% sets in the extended finite-element method." Computer methods in applied 
% mechanics and engineering 190, no. 46-47 (2001): 6183-6200. The material 
% constants are constant in two materials. The linear displacement field:
% u1=x1, u2=x2 is imposed on the outer boundary of the circular plate.
% TODO: Add description of test case to Wiki

        testnum = testnum + 1;
        testname = 'Circular inclusion as a bimaterial problem';
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
        global calc_type
            config_name = 'CircularInclusion';
            % run with nodal averaging
            calc_type = 'center';
%       main  % Runs calculation
        main_nonlinear % Runs calculation
        
        test_pass = 1;
        fprintf('\nPASS')
        testpasssummary(testnum) = test_pass;
            
        % Step 4 - Cleanup
        clearvars -except  curDir  ConfigDir ...
                      ntests testpasssummary testnum nameslist...
                      plot2vtk VTKFolder progress_on