function [Mesh, Material, BC, Control] = BoreholePressure_hypermesh(config_dir, progress_on)


%% Material Properties (Solid)

    % NOTES-------------------------------------------------------------
                                
        % NOTE: anonymous functions are defined with respect to the variable x,
        % which is a vector [x(1) x(2) x(3)] = [x y z]

        % NOTE: Material properties must be continuous along an element, 
        % otherwise, quadrature order must be increased significantly
        
        % NOTE: Number of material properties can be more than one. Properties
        % for different materials are saved in Material.Prop.
        % For example, Young's modulus and Poisson's ratio of ith material will be saved in
        % Material.Prop(i).E and Material.Prop(i).nu, respectively.

    % Specify Material Model
        % LE1 - Linear elasticity
        % LET1 - Linear elastic with mass based damping
        % LED1 - Dynamic linear elasticity
        % ST1 - Stiffening model with 1st invariant of strain
        % ST2 - Softening model with 1st invariant of strain
        % TR2 - Stiffening model with mass based damping with 1st invariant of strain
        % VE1 - Viscoelaticity with stiffness based damping
        % TH1 - Thermal Diffusion (Steady-State)
        % TH2 - Thermal Diffusion (Transient)
    Material.Model = 'LE1';
    
    % number of material properties
    Material.nmp = 1;
        
    % Properties material 1
    Material.Prop(1).E0 = 2e11; % Young's modulus [Pa]
    Material.Prop(1).nu = 0.3; % Poisson's ratio

    % Constitutive law: 'PlaneStrain' or 'PlaneStress' 
    Material.Dtype = 'PlaneStress'; 

    % Thickness (set as default to 1)
    Material.t = @(x) 1;

    % Alternatively, import a material file
    % Material = Material_shale();
	
    [Material, ~, ~] = setMaterialModel(Material);


%% Mesh Properties
    if progress_on
        disp([num2str(toc),': Building Mesh...']);
    end

    meshFileName = 'Mesh Files\borehole_pressure.fem';
%     meshFileName = 'Mesh Files\borehole_pressure_finer.fem';
    % number of space dimensions 
    nsd = 2;

    Mesh = BuildMesh_imported(meshFileName, nsd, config_dir, progress_on, 0, Material.ProblemType);

%% Assign Materials to Mesh
    % type of material per element
    Mesh.MatList = zeros(Mesh.ne, 1, 'int8');
    
    % assign material type to elements
    Mesh.MatList(:) = 1;


%% Boundary Conditions

    % Dirichlet boundary conditions (essential)
    % Note on test case:
    %       As an example, the present case solves the problem for two sets
    %       of both Dirichlet and Neumann BC.
    % -----------------------------------------------------------------
        % column vector of prescribed displacement dof  
        temp_DOF = [Mesh.DOF(Mesh.c_BC_nE{1},:); Mesh.DOF(Mesh.c_BC_nE{2},:)];
        temp_BC_flag = [Mesh.c_BC_E{1}; Mesh.c_BC_E{2}];
        temp = temp_DOF.*temp_BC_flag;
        BC.fix_disp_dof = nonzeros(reshape(temp, length(temp)*Mesh.nsd,1));

        % prescribed displacement for each dof [u1; u2; ...] [m]
        BC.fix_disp_value = @(t) zeros(length(BC.fix_disp_dof),1);  

    %% Neumann BC
    % -----------------------------------------------------------------
        % column vector of prescribed traction dofs
        BC.traction_force_dof = [];

        % magnitude of prescribed tractions [N]
        BC.traction_force_dof_value = [];
        
        %--------------------------------------------------------------
        % Edge elements and normals are obtained from Mesh structure.
        % Tractions boundary conditions can be applied as:
        %    → function for vector on x and y 
        %    → function for vector on n and s (normal
        %           and tangential directions of each point of the boundary)
        % These functions are stored in BC.c_N_t_f
        % To recognize if the function is case 1 or 2, BC.c_N_t_flag must include
        % a vector with:
        %   0 → function for vector on x and y
        %   1 → function for vector on n and s
        % -------------------------------------------------------------
        % Note: The total number of traction functions must be compatible
        % with the number of sets defined for edge elements
        % -------------------------------------------------------------
        BC.c_N_t_f = cell(1,1);
        BC.c_N_t_f{1} = @(x,t) [-1;0];
  
        BC.c_N_t_flag = [1,1];
        %--------------------------------------------------------------
        % prescribed traction [t1x t1y;t2x t2y;...] [N]
        BC.traction_force_node = [];
        
        BC.traction_force_value = [];

        %--------------------------------------------------------------
        % body force     
        BC.b = @(x)[];    

%% Computation controls

        % quadrature order
        Control.qo = 2;

        % Calculation of values for discontinuous variables 
        % (i.e. stress/strain)
        % 'none': calculated at each node for each element separately; 
        %           no output in vtk
        % 'nodal': averaged at each node for all elements attached to 
        %           the node; output as nodal values in vtk
        % 'center': calculated at the center of each element; output as 
        %           single value for each element in vtk
        % 'L2projection': Least squares projection of stress and strain,
        %           output as nodal values
        global calc_type
        Control.stress_calc = calc_type;

        % penalty parameter for solution of static problem with 
        % LinearSolver3
        Control.beta = 10^10;

        % method used for solving linear problem:
        % 'LinearSolver1': Partitioning
        % 'LinearSolver2': Zeroing DOFs in stiffness matrix 
        %                   corresponding to essential boundaries
        % 'LinearSolver3': Penalty method
        Control.LinearSolver = 'LinearSolver1'; 

        % parallel inversion
        % Use parallel processing to invert the matrix.
        % Usually more efficient at 2e5 dofs
        Control.parallel = 2;
        
        % time integration parameter
        % for 1st order problem (transient diffusion, viscoelastic)
        % 1 = Backward Euler, 0.5 = Crank-Nicolson
        % for 2nd order problem (dynamic)
        % range = [-1/3, 0], use 0 by default
        Control.alpha = 0.5; 
        
        % Newton Raphson controls
        Control.r_tol = 1e-5; % Tolerance on residual forces
        Control.iter_max = 50; % Maximum number of iteration in Newton Raphson algorithm
        
end