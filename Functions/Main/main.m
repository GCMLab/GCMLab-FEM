% Static solid deformation of linear elastic material

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
    
    % Compute external force vector
    if progress_on
        disp([num2str(toc),': Compute Force Vector...']);
    end
    Fext = getFext(Mesh, BC, Quad);

%% Define initial conditions
    % For static solution, always assume that free nodes are initially at
    % zero (undisplaced)
    d0 = zeros(Mesh.nDOF, 1);
    d0(BC.fix_disp_dof) = BC.fix_disp_value;

%% Solve linear system of equations
    switch Control.LinearSolver
        case 'LinearSolver1'
            [d, f_fixed] = LinearSolver1(K, Fext, BC.fix_disp_value, ...
                                            BC.free, BC.fixed);
        case 'LinearSolver2'
            [d, f_fixed] = LinearSolver2(K,Fext,BC.fix_disp_value, ...
                                            BC.free, BC.fixed);
        case 'LinearSolver3'
            [d, f_fixed] = LinearSolver3(K, Fext, BC.fix_disp_value, ...
                                            BC.free, BC.fixed, Control.beta);
    end
    Fext(BC.fixed) = f_fixed;

% Strain
    [strain, stress] = getStrain(d, Mesh, Material, Control.stress_calc, Quad);   

% Force vectors
    Fint = K*d;

%% Write to vtk
    if plot2vtk
        write2vtk_static(config_name, vtk_dir, Mesh, Control, BC.fixed, d, strain, stress, ...
                        Fint, Fext);
    end
    
    if progress_on
        disp('done')
    end