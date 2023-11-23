% ------------------------------------------------------------------------
% Runs unit Test 32 - Borehole pressure
% ------------------------------------------------------------------------
% Runs borehole case wih pressure. Analytical solution obtained from book:
%   Analytic Methods in Geomechanics by Kam-tim Chau (2013)
%   Link: https://www.taylorfrancis.com/books/mono/10.1201/9781315275277/analytic-methods-geomechanics-kam-tim-chau

        testnum = testnum + 1;
        testname = 'Borehole subject to internal pressure';
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
            config_name = 'BoreholePressure_hypermesh';
            calc_type = 'L2projection';
            main_nonlinear % Runs calculation
            
            stress_L2 = stress;
            
            config_name = 'BoreholePressure_hypermesh';
            calc_type = 'nodal';
            main_nonlinear % Runs calculation
            
            stress_nodal = stress;
           
        % Step 2 - Check results
        p = BC.c_N_t_f{1}(1);
        error = BoreholePressure_check(Mesh,d,p(1),stress_L2, stress_nodal);
        
        if error_L2 < error_nodal
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