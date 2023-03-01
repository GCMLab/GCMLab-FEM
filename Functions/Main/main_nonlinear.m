% Quasi-Static solid deformation of nonlinear elastic material

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
    [Mesh, Material, BC, Control] = ...
            feval(config_name, ConfigDir, progress_on);
        


%% Set Default Values
    [Mesh, Material, BC, Control] = setDefaults(Mesh, Material, BC, Control);

%% Check for valid inputs
    if progress_on
        fprintf('%.2f: Checking for valid inputs...\n', toc);
    end
    [Mesh, Material, BC, Control] = cleanInput(Mesh, Material, BC, Control);

%% Initialize time variables
    t = Control.StartTime;
    dt = Control.TimeStep;
    t = t + dt;
    step_count = 0;

%% Identify free and fixed dofs
    BC.fixed = BC.fix_disp_dof;
    BC.free = setdiff(Mesh.DOF, BC.fixed)';

%% Quadrature calculation
    % (Done at beginning of code - this assumes all elements are the 
    % same type and same order of quadrature for all integrals)
    Quad = GlobalQuad(Mesh.nsd, Mesh.type, Control.qo);

%% Compute system matrices

    % Compute linear elastic stiffness matrix
    if progress_on
        disp([num2str(toc),': Assembling Linear Elastic Stiffness Matrix...']);
    end
    Klin = getK(Mesh, Quad, Material); % Linear elastic stiffness matrix
    

%% Define initial conditions
    
    d0 = BC.IC;
    d0(BC.fix_disp_dof) = BC.fix_disp_value(t-dt);
    Fext = getFext(Mesh, BC, Quad,t-dt);
    d = d0;
    K = Klin;
    
    % Export initial conditions
        % Strain
        [strain, stress] = getStrain(d0, Mesh, Material, Control.stress_calc, Quad);   

    % Internal force vectors
        Fint = K*d0;

    % Write initial conditions to vtk
        if plot2vtk
            write2vtk_quasistatic(config_name, vtk_dir, Mesh, Control, BC.fixed, d0, strain, stress, ...
                            Fint, Fext, step_count);
        end
        step_count = step_count + 1;
        
%% Initialize variable tracking for time-dependent problems 
% Not recommended, primarily for use in testing and debugging.
% Output is saved in vtk files.

if Control.dSave
    n_timesteps = ceil((Control.EndTime - Control.StartTime)/dt);
    dSave = zeros(length(d0),n_timesteps+1);
    dSave(:,1) = d0;                                   % Save displacements
    sSave = zeros(3,Mesh.nn, n_timesteps+1);        % Save stresses
    sSave(:,:,1) = stress;
end
        
 %% Solve the time-dependent nonlinear problem
 t_tol = 1e-10; % Tolerance on the final end time
 r_tol = 1e-5; % Tolerance on residual forces
 iter_max = 50; % Maximum number of iteration in Newton-Raphson algorithm
 while true   
    if progress_on
        fprintf(['\n', num2str(toc),': Computing Time = %.2f s, Timestep %d, Progress %.2f percent...\n'], t, step_count, t/Control.EndTime*100);
    end
    
    % Initialize incremental variables in each step
    Dd = zeros(Mesh.nsd*Mesh.nn,1); % Increment of displacement vector
    Dd(BC.fix_disp_dof) = BC.fix_disp_value(t) - BC.fix_disp_value(t-dt);
    iter = 1; % Iteration counter in each step

    % Compute external force and residual vectors
    if progress_on
        disp([num2str(toc),': Compute Force Vector...']);
    end
    Fext = getFext(Mesh, BC, Quad, t);
    ResForce = Fext - Fint;
    % Todo: Make the boundary conditions time dependent functions
    
    while true
        % Solve incremental form of system of equations
        if progress_on
            disp([num2str(toc),': Solve System of Equations...']);
        end
        switch Control.LinearSolver
            case 'NonLinearSolver1'
                [Dd, ~] = NonLinearSolver1(K, ResForce, Dd(BC.fix_disp_dof), ...
                                                BC.free, BC.fixed, Control.parallel);
            case 'NonLinearSolver2'
                [Dd, ~] = NonLinearSolver2(K,ResForce,Dd, ...
                                                BC.free, BC.fixed, Control.parallel);
            case 'NonLinearSolver3'
                [Dd, ~] = NonLinearSolver3(K, ResForce, Dd, ...
                                                BC.free, BC.fixed, Control.beta, Control.parallel);
        end
        
        % Update displacement vector
        d = d + Dd;
        
        % Compute nonlinear stiffness matrix
        if progress_on
            disp([num2str(toc),': Assembling Nonlinear Stiffness Matrix...']);
        end
        Knlin = sparse(Mesh.nDOF, Mesh.nDOF); % Nonlinear stiffness matrix
        K = Klin + Knlin ; % Total stiffness matrix
        % Todo: Write a function to calculate Knlin when the material is nonlinear
        
        % Internal force vectors
        Fint = K*d;
        % Todo: Create a function for calculation of Fint when the material is nonlinear
        
        % Calculate residual vector at the end of each iteration
        Dd = zeros(Mesh.nsd*Mesh.nn,1);
        ResForce = Fext - Fint;
        ResForce(BC.fixed) = 0;
        
        % Calculate the norm of residual vector
        if norm(Fint) < 1e-2
            resScale = 1+norm(FintPre);
        else
            resScale = 1+norm(Fint);
        end
        res = norm(ResForce)/resScale;
        if res < r_tol
            break
        end
        
        iter = iter + 1;
        if iter > iter_max
            err_NR = sprintf('\t Newton-Raphson algorithm will not converge: The number of iterations in Newton-Raphson algorithm exceeds the maximum number of iteration');
            error(err_NR)
        end
        
    end
    
    if progress_on
        fprintf([num2str(toc),': Timestep %d converged with %d iterations ...\n'], step_count, iter);
    end
      
    Fext(BC.fixed) = Fint(BC.fixed);
    FintPre = Fint; 

    % Strain
        if progress_on
            disp([num2str(toc),': Post-Processing...']);
        end
        [strain, stress] = getStrain(d, Mesh, Material, Control.stress_calc, Quad);   
    % Todo: Create a new function for calculation of stress and strain in nonlinear cases


    % Write to vtk
        if plot2vtk
            write2vtk_quasistatic(config_name, vtk_dir, Mesh, Control, BC.fixed, d, strain, stress, ...
                            Fint, Fext, step_count);
        end
        
        if Control.dSave
           dSave(:,step_count+1) = d; 
           sSave(:,:,step_count + 1) = stress;
        end
        
     % Update time variables
     if abs(t-Control.EndTime) < t_tol
        break
     end
     
     t = t + dt;
     step_count = step_count + 1;
     if t > Control.EndTime
        t = t - dt;
        dt = Control.EndTime - t;
        t = Control.EndTime;
     end
     
     
 end
     if Control.dSave
         d = dSave; 
     end
 
     if progress_on
        disp('done')
    end