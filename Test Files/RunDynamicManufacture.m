% ------------------------------------------------------------------------
% Runs unit test - Dynamic Q4 Plane stress manufactured solution convergence as a
% part of RunTests 
% ------------------------------------------------------------------------
% Calculates the convergence rates of a uniform Q4 mesh using a
% manufactured solution in which 
% u := (x1, x2, t) -> -sin(1/2*pi*x1)*sin(1/2*pi*x2)*sin(2*pi*t)
% v := (x1, x2, t) -> cos(1/2*pi*x1)*cos(1/2*pi*x2)*cos(2*pi*t)
% under plane stress condition 

        testnum = testnum + 1;
        testname = 'Dynamic Plane stress manufacture solution';
        nameslist{testnum} = testname;
       
        
        % test runs 3 meshes, only finest mesh will be saved
    
        fprintf('\n\n Test %d : %s\n', testnum, testname)
        
        % Step 1 - Run Simulation
        global nex quadorder E nu rho alpha tf n_steps plot_on
            alpha = -1/3 ; %           α [-1/3, 0]
            E = 2e11;
            nu = 0.3;
            quadorder = 2; 
            rho = 2400;
            tf = 3.123;
            n_steps = 100;
            plot_on = 0;

            config_name = 'ManufacturedSolution_Dynamic';
            
            
            % Run coarse mesh
                        
                % Create test VTK folder
                if plot2vtk
                    folname = ['\Test',num2str(testnum),'_mesh_1'];
                    vtk_dir = fullfile(VTKFolder,folname);
                    if ~isfolder(vtk_dir) 
                        mkdir(vtk_dir)
                    end
                end

                nex = [2;2]*1;
                main_nonlinear % Runs calculation

                d_coarse = d;
                stress_coarse = stress;
                strain_coarse = strain;
                Mesh_coarse = Mesh;

            % Run fine mesh
            
                 % Create test VTK folder
                if plot2vtk
                    folname = ['\Test',num2str(testnum),'_mesh_2'];
                    vtk_dir = fullfile(VTKFolder,folname);
                    if ~isfolder(vtk_dir) 
                        mkdir(vtk_dir)
                    end
                end
            
                nex = [2;2]*5;
                main_nonlinear % Runs calculation

                d_fine = d;
                stress_fine = stress;
                strain_fine = strain;
                Mesh_fine = Mesh;

            % Run finer mesh
            
                % Create test VTK folder
                if plot2vtk
                    folname = ['\Test',num2str(testnum),'_mesh_3'];
                    vtk_dir = fullfile(VTKFolder,folname);
                    if ~isfolder(vtk_dir) 
                        mkdir(vtk_dir)
                    end
                end
                
                nex = [2;2]*10;
                main_nonlinear % Runs calculation

                d_finer = d;
                stress_finer = stress;
                strain_finer = strain;
                Mesh_finer = Mesh;
            
        
        % Step 2 - Check results
            %Visual comparison of displacement of the center node
            [disp_err, disp_err_stored] = ManufacturedSolution_Dynamic_point_check(d_coarse, d_fine, d_finer,...
                Mesh_coarse, Mesh_fine, Mesh_finer, Control);
            
            %Convergence study
            % Not fully working
%             [m_L2, m_e] =ManufacturedSolution_Dynamic_check(d_coarse, d_fine, d_finer, ...
%                 stress_coarse, stress_fine, stress_finer, strain_coarse, strain_fine, ...
%                 strain_finer, Mesh_coarse, Mesh_fine, Mesh_finer);
%             
%             fprintf('\nQ4 L2-norm converges at a rate of %.2f',m_L2)
%             fprintf('\nQ4  e-norm converges at a rate of %.2f',m_e)
            
            convergence_exact = 0.15;
            convergence_num = 0.001;
            
            if convergence_exact >= disp_err && convergence_num >= disp_err_stored
                test_pass = 1;
            else
                test_pass = 0;
            end
            
        % Step 3 - Output results
            if test_pass
                fprintf('\nPASS\n')
            else
                fprintf('\n\nFAIL\n')
            end
        testpasssummary(testnum) = test_pass;

        
        % Step 4 - Cleanup
        clearvars -except  curDir  ConfigDir ...
                      ntests testpasssummary testnum nameslist...
                      plot2vtk VTKFolder progress_on