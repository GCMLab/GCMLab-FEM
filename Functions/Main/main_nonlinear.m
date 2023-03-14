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
    dtnm1 = dt;     % delta t bewteen timesteps n-2 to n-1
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
    
    % Compute mass matrix
    M = 0; % placeholder
    
    % Create tangent matrix function pointer
    [~,stiffnessmatrix_name] = fileparts(Material.StiffnessMatrix);
    

%% Define initial conditions
    
    d0 = BC.IC;
    d0(BC.fix_disp_dof) = BC.fix_disp_value(t-dt);
    Fext = getFext(Mesh, BC, Quad,t-dt);
    d = d0;     % d at timestep n
    dnm1 = d0;  % d at timestep n-1
    dnm2 = d0;  % d at timestep n-2
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
    
    if progress_on
        msg_len = 1;
        fprintf('Newton-Raphson - Iteration:  ')
    end
    
    while true
        if progress_on
            fprintf(repmat('\b',1,msg_len)) 
            fprintf('%d',iter');
            msg_len = numel(num2str(iter));
        end
        
        % Solve incremental form of system of equations
        switch Control.LinearSolver
            case 'LinearSolver1'
                [Dd, ~] = LinearSolver1(K, ResForce, Dd(BC.fix_disp_dof), ...
                                                BC.free, BC.fixed, Control.parallel);
            case 'LinearSolver2'
                [Dd, ~] = LinearSolver2(K,ResForce,Dd, ...
                                                BC.free, BC.fixed, Control.parallel);
            case 'LinearSolver3'
                [Dd, ~] = LinearSolver3(K, ResForce, Dd, ...
                                                BC.free, BC.fixed, Control.beta, Control.parallel);
        end
        
        % Update displacement vector
        d = d + Dd;
        
        % Compute nonlinear stiffness matrix and internal forces
        [K, Fint] = feval(stiffnessmatrix_name, Mesh, Quad, Material, Klin, M, d, dnm1, dnm2, stress, strain, dt, dtnm1) ; 

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
        
        % Calculate residual
        res = norm(ResForce)/resScale;
        
        % Check convergence
        if res < Control.r_tol
            break
        end
        
        iter = iter + 1;
        % Case of divergence in Newton Raphson algorithm
        if iter > Control.iter_max
            err_NR = sprintf('\n \t Newton-Raphson algorithm will not converge: The number of iterations in Newton-Raphson algorithm exceeds the maximum number of iteration');
            error(err_NR)
        end
        
    end
    
    if progress_on
        fprintf(['\n', num2str(toc),': Timestep %d converged with %d iterations ...\n'], step_count, iter);
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
     
     dtnm1 = dt;
     t = t + dt;
     step_count = step_count + 1;
     if t > Control.EndTime
        t = t - dt;
        dt = Control.EndTime - t;
        t = Control.EndTime;
     end
     
     % Update previous d vectors
     dnm2 = dnm1;
     dnm1 = d;
     
     
 end
     if Control.dSave
         d = dSave; 
     end
 
     if progress_on
        disp('done')
    end