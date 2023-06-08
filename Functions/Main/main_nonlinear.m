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
    switch Control.TimeCase 
        case 'dynamic'
            if progress_on
                disp([num2str(toc),': Assembling Mass Matrix...']);
            end
            M = getC(Mesh, Quad, Material); % Dynamic Case
        otherwise
            M = sparse(Mesh.nDOF, Mesh.nDOF); % Static and Transient Case
    end 
    
    % Create tangent matrix function pointer
    [~,stiffnessmatrixfile_name] = fileparts(Material.StiffnessMatrixFile);
    [~,stressstrainfile_name] = fileparts(Material.StressStrainFile);

    % Compute linear damping stiffness matrix
    switch Control.TimeCase 
        case 'static'
            C = sparse(Mesh.nDOF, Mesh.nDOF); % Static Case
        otherwise
            % for Control.TimeCase = 'transient' or 'dynamic'
            if progress_on
                disp([num2str(toc),': Assembling Mass Matrix...']);
            end
            C = getC(Mesh, Quad, Material); % Transient Case          
    end 
 
    
%% Define initial conditions
    d0 = BC.IC_d; % Initial condition for displacement
    v0 = BC.IC_v; % Initial condition for velocity
    a0 = BC.IC_a; % Initial condition for acceleration
    d0(BC.fix_disp_dof) = BC.fix_disp_value(t-dt);
    v0(BC.fix_vel_dof) = BC.fix_vel_value(t-dt);
    a0(BC.fix_acc_dof) = BC.fix_acc_value(t-dt);
    Fext = getFext(Mesh, BC, Quad,t-dt);
    Fextnm1 = Fext; % Fext at timestep n-1
    d = d0;     % d at timestep n
    dnm1 = d0;  % d at timestep n-1
    dnm2 = d0;  % d at timestep n-2
    v = v0;     % v at timestep n
    vnm1 = v0;  % v at timestep n-1
    vnm2 = v0;  % v at timestep n-2
    a = a0;     % a at timestep n
    anm1 = a0;  % a at timestep n-1
    anm2 = a0;  % a at timestep n-2
    K = Klin;
    
    % Export initial conditions
        % Strain
        [strain, stress] = getStrain(d0, Mesh, Material, Control.stress_calc, Quad);

    % Internal force vectors
        switch Control.Timecase
            case 'static'
                Fint = K*d0;
            case 'transient'
                Fint = (Control.alpha*Klin+(1/dt)*C)*d+((1-Control.alpha)*Klin-C./dt)*dnm1;
            case 'dynamic'
                gam = 1/2-Control.alpha;
                bet = (1-Control.alpha)^2/4;
                d_temp = dnm1 + dt*vnm1 + dt^2/2*(1-2*bet)*anm1;
                v_temp = vnme + dt*(1-gam)*anm1;
                Fint = (1+Control.alpha)*C*(v_temp-gam*d_temp/(dt*bet)+gam*d/(dt*bet))...
                    +(1+Control.alpha)*K*d - Control.alpha*(C*vnm1 + K*dnm1);
        end

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
    dSave(:,1) = d0;                                % Save displacements
    switch Mesh.nsd
        case 1
            dim = 1;
        case 2
            dim = 3;
        case 3
            dim = 6;
    end
    sSave = zeros(dim,Mesh.nn, n_timesteps+1);        % Save stresses
    sSave(:,:,1) = stress;
    loadSave = zeros(length(d0),n_timesteps+1);     % Save applied load
end
  
 %% Solve the time-dependent nonlinear problem
 % Set tolerance on the final end time
    t_tol = 1e-10; 
 while true   % Timestep loop
    % Output progress
        if progress_on
            fprintf(['\n', num2str(toc),': Computing Time = %.2f s, Timestep %d, Progress %.2f percent...\n'], t, step_count, t/Control.EndTime*100);
        end
    
    % Initialize incremental variables in each step
        Dd = zeros(Mesh.nsd*Mesh.nn,1);                                         % Increment of displacement vector
        d(BC.fix_disp_dof) = BC.fix_disp_value(t);                              % Apply fixed DoF Boundary conditions
        iter = 1;                                                               % Iteration counter in each step
        converged = 0;                                                          % convergence status
        FintPrev = 0;                                                            % Internal force vector from previous iteration

    % Compute external force and residual vectors
        if progress_on
            disp([num2str(toc),': Compute Force Vector...']);
        end
        Fext = getFext(Mesh, BC, Quad, t);
    
    % Output progress to command window
        if progress_on
            msg_len = 1;
            fprintf('Newton-Raphson - Iteration:  ')
        end
    
    % Begin Newton-Raphson Loop
    while ~converged
        % Update iteration counter in command window
            if progress_on
                fprintf(repmat('\b',1,msg_len)) 
                fprintf('%d',iter');
                msg_len = numel(num2str(iter));
            end
      
        % Compute nonlinear stiffness matrix and internal forces
            [K, ResForce, Fint] = feval(stiffnessmatrixfile_name, Mesh, Quad, Material, Fext, Fextnm1, Klin, M, d, dnm1, dnm2, dt, dtnm1, C, Control.alpha); 
            ResForce(BC.fixed) = 0;
        
        % Calculate the norm of residual vector
            if norm(Fint) < 1e-2
                resScale = 1+norm(FintPrev);
            else
                resScale = 1+norm(Fint);
            end
        
        % Calculate residual
            res = norm(ResForce)/resScale;
        
        % Check convergence
            if res < Control.r_tol
                converged = 1;
            else
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

                % Update iteration number
                    iter = iter + 1;
                % Update internal force vector for normalization
                    FintPrev = Fint; 


            end
              
        
        % Case of divergence in Newton Raphson algorithm
            if iter > Control.iter_max 
                err_NR = sprintf('\n \t Newton-Raphson algorithm will not converge: The number of iterations in Newton-Raphson algorithm exceeds the maximum number of iteration');
                error(err_NR)
            end
        
    end
    
    % Output progress
        if progress_on
            fprintf(['\n', num2str(toc),': Timestep %d converged with %d iterations ...\n'], step_count, iter);
        end

    % Strain calculation
        if progress_on
            disp([num2str(toc),': Post-Processing...']);
        end
        [strain, stress] = feval(stressstrainfile_name, d, Mesh, Material, Control.stress_calc, Quad, dnm1);   


    % Write to vtk
        Fext(BC.fixed) = Fint(BC.fixed);   % Set external forces as equal to reaction forces at fixed dof for output
        
        if plot2vtk
            write2vtk_quasistatic(config_name, vtk_dir, Mesh, Control, BC.fixed, d, strain, stress, ...
                            Fint, Fext, step_count);
        end
        
        if Control.dSave
           dSave(:,step_count+1) = d; 
           sSave(:,:,step_count + 1) = stress;
           loadSave(:,step_count+1) = Fext;
        end
        
     
     % Break out of loop at end time
         if abs(t-Control.EndTime) < t_tol
            break
         end
         
     % Update time variables  
         dtnm1 = dt;                    % Timestep from previous step
         t = t + dt;                    % New time n 
         step_count = step_count + 1;   
         if t > Control.EndTime         % Catch end time overshoot
            t = t - dt;
            dt = Control.EndTime - t;
            t = Control.EndTime;
         end
     
     % Update vectors from previous timesteps
     dnm2 = dnm1;                       % d vector from timestep n-2
     dnm1 = d;                          % d vector from timestep n-1
     Fextnm1 = Fext;                    % Fext from timestep n-1
 
 end
     if Control.dSave
         d = dSave; 
     end
     
     if Control.plotLoadDispl
         plotLoadVsDispl(loadSave, dSave, Control);
     end
 
     if progress_on
        disp('done')
     end
