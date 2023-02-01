% Quasi-Static solid deformation of linear elastic material

% Acknowledgements: Matin Parchei Esfahani


%% Delete past vtk files (so don't overwrite any)
    if plot2vtk
        if progress_on
            fprintf('%.2f: Deleting old vtk files...\n', toc);
        end
        sol = fullfile(vtk_dir, '*');
        delete(sol)
    end

%% Input problem definition
    if progress_on
        fprintf('%.2f: Reading config file...\n', toc);
    end
    [Mesh, Material, BC, Control, IC] = ...
            feval(config_name, ConfigDir, progress_on);
        
%% Initialize time variables
    t = Control.StartTime;
    dt = Control.TimeStep;
    t = t + dt;
    step_count = 0;

%% Set Default Values
    [Mesh, Material, BC, Control] = setDefaults(Mesh, Material, BC, Control);

%% Check for valid inputs
    if progress_on
        fprintf('%.2f: Checking for valid inputs...\n', toc);
    end
    [Mesh, Material, BC, Control] = cleanInput(Mesh, Material, BC, Control);


%% Identify free and fixed dofs
    BC.fixed = BC.fix_disp_dof;
    BC.free = setdiff(Mesh.DOF, BC.fixed)';

%% Quadrature calculation
    % (Done at beginning of code - this assumes all elements are the 
    % same type and same order of quadrature for all integrals)
    Quad = GlobalQuad(Mesh.nsd, Mesh.type, Control.qo);

%% Compute system matrices

    % Compute stiffness matrix
    if progress_on
        disp([num2str(toc),': Assembling Stiffness Matrix...']);
    end
    K = getK(Mesh, Quad, Material);
    

%% Define initial conditions
    
    d0 = IC.d0;
    d0(BC.fix_disp_dof) = BC.fix_disp_value;
    Fext = getFext(Mesh, BC, Quad,t);
    
    % Export initial conditions
        % Strain
        [strain, stress] = getStrain(d0, Mesh, Material, Control.stress_calc, Quad);   

    % Force vectors
        Fint = K*d0;

    % Write to vtk
        if plot2vtk
            write2vtk_quasistatic(config_name, vtk_dir, Mesh, Control, BC.fixed, d0, strain, stress, ...
                            Fint, Fext, step_count);
        end
        step_count = step_count + 1;
        
 %% Solve the time-dependent problem
 while t <= Control.EndTime
    if progress_on
        fprintf(['\n', num2str(toc),': Computing Time = %.2f s, Timestep %d, Progress %.2f percent...\n'], t, step_count, t/Control.EndTime*100);
    end
    % Compute external force vector
    if progress_on
        disp([num2str(toc),': Compute Force Vector...']);
    end
    Fext = getFext(Mesh, BC, Quad, t);
    % Todo: Make the boundary conditions time dependent functions
     
     % Solve linear system of equations
        if progress_on
            disp([num2str(toc),': Solve System of Equations...']);
        end
        switch Control.LinearSolver
            case 'LinearSolver1'
                [d, f_fixed] = LinearSolver1(K, Fext, BC.fix_disp_value, ...
                                                BC.free, BC.fixed, Control.parallel);
            case 'LinearSolver2'
                [d, f_fixed] = LinearSolver2(K,Fext,BC.fix_disp_value, ...
                                                BC.free, BC.fixed, Control.parallel);
            case 'LinearSolver3'
                [d, f_fixed] = LinearSolver3(K, Fext, BC.fix_disp_value, ...
                                                BC.free, BC.fixed, Control.beta, Control.parallel);
        end
        Fext(BC.fixed) = f_fixed;

    % Strain
        if progress_on
            disp([num2str(toc),': Post-Processing...']);
        end
        [strain, stress] = getStrain(d, Mesh, Material, Control.stress_calc, Quad);   

    % Force vectors
        Fint = K*d;

    % Write to vtk
        if plot2vtk
            write2vtk_quasistatic(config_name, vtk_dir, Mesh, Control, BC.fixed, d, strain, stress, ...
                            Fint, Fext, step_count);
        end
     % Update time variables
     t = t + dt;
     step_count = step_count + 1;
     
 end
 
     if progress_on
        disp('done')
    end