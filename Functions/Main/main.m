% Static solid deformation of linear elastic material

% Acknowledgements: Matin Parchei Esfahani, Endrina Rivas

%% Delete past vtk files (so don't overwrite any)
    if Control.vtk
        fprintf('%.2f: Deleting old vtk files...\n', toc);
        sol = fullfile(Control.vtk_dir, '*');
        delete(sol)
    end

%% Input problem definition
    fprintf('%.2f: Reading config file...\n', toc);
    [Mesh, Material, BC, Control] = ...
            feval(Control.config_name, Control.config_dir);

%% Set Default Values
    [Mesh, Material, BC, Control] = setDefaults(Mesh, Material, BC, Control);

%% Identify free and fixed dofs
    BC.fixed = BC.fix_disp_dof;
    BC.free = setdiff(Mesh.DOF, BC.fixed)';

%% Quadrature calculation
    % (Done at beginning of code - this assumes all elements are the 
    % same type and same order of quadrature for all integrals)
    Quad = GlobalQuad(Mesh.nsd, Mesh.type, Control.qo);

%% Compute system matrices

    % Compute stiffness matrix
    disp([num2str(toc),': Assembling Stiffness Matrix...']);
    K = getK(Mesh, Quad, Material);
    
    
    % Compute external force vector
    disp([num2str(toc),': Compute Force Vector...']);
    Fext = getFext(Mesh, BC, Material, Quad, Control);

%% Define initial conditions
    % For static solution, always assume that free nodes are initially at
    % zero (undisplaced)
    d0 = zeros(Mesh.nDOF, 1);
    d0(BC.fix_disp_dof) = BC.fix_disp_value;

%% Solve linear system of equations
    switch Control.LinearSolver
        case 'LinearSolver1'
            [d, f_fixed] = LinearSolver1(K, Fext, BC.fix_disp_value, ...
                                            Mesh.nDOF, BC.free, BC.fixed);
        case 'LinearSolver2'
            [d, f_fixed] = LinearSolver2(K,Fext,BC.fix_disp_value, ...
                                            BC.free, BC.fixed);
        case 'LinearSolver3'
            [d, f_fixed] = LinearSolver3(K, Fext, BC.fix_disp_value, ...
                                            BC.free, BC.fixed, Control.beta);
    end
    Fext(BC.fixed) = f_fixed;

% Strain
    [strain, stress] = getStrain(d, Mesh, Control, Material);   

% Force vectors
    Fint = K*d;

%% Write to vtk
    if Control.vtk
        write2vtk_static(Mesh, Control, BC.fixed, d, strain, stress, ...
                        Fint, Fext);
    end
    
disp('done')