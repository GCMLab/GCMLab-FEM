function [Mesh, Material, BC, Control] = PatchTestA(config_dir, progress_on)
    global E nu traction quadorder meshfilename
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
    Material.Prop(1).E0 = E; % Young's modulus [Pa]
    Material.Prop(1).nu = nu; % Poisson's ratio

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
    
    % Mesh formats: 
    %   'MANUAL'- In-house structured meshing
    % 	'IMPORTED'  - Import .msh file from GMSH, or .fem from HYPERMESH structured or unstructured
    %   'EXCEL' - Import .xlsx file, structured or unstructured
    MeshType = 'IMPORTED';        

    
    switch MeshType
        case 'MANUAL'
            % location of initial node [m] [x0;y0;z0] 
            x1 = [0;0;0];
            % number of space dimensions 
            nsd = 2;
            % size of domain [m] [Lx;Ly;Lz] 
            L = [1;1];
            % number of elements in each direction [nex; ney; nez] 
            nex = [2;2]*20;
            % element type ('Q4')
            type = 'Q4';
            
            Mesh = BuildMesh_structured(nsd, x1, L, nex, type, progress_on, Material.ProblemType);
        case 'IMPORTED'
            % Allows input of files from GMSH
            % Note: the only currently supported .msh file formatting is
            % Version 2 ASCII
            % Ctrl + e to export the mesh, specify extension .msh, specify
            % format Version 2 ASCII
            meshFileName = meshfilename;
            % number of space dimensions 
            nsd = 2;
            
            Mesh = BuildMesh_imported(meshFileName, nsd, config_dir, progress_on, 0, Material.ProblemType);
        case 'EXCEL'
            meshFileName = 'CircularInclusion.xlsx';
            % number of space dimensions
            nsd = 2;
            
            Mesh = BuildMesh_EXCEL(meshFileName, nsd, config_dir, progress_on, Material.ProblemType);
    end
    



%% Assign Materials to Mesh
    % type of material per element
    Mesh.MatList = zeros(Mesh.ne, 1, 'int8');
    
    % assign material type to elements
    Mesh.MatList(:) = 1;


%% Boundary Conditions
    % {TIPS}------------------------------------------------------------
        % TIP selecting edges:
        % bottom_nodes = find(Mesh.x(:,2)==0); 
        % top_nodes = find(Mesh.x(:,2)==2);
        % left_nodes = find(Mesh.x(:,1)==0);
        % right_nodes = find(Mesh.x(:,1)==4);
        % bottom_dof = [bottom_nodes*2 - 1; bottom_nodes*2];
        % top_dof = [top_nodes*2 - 1;top_nodes*2];

    % Dirichlet boundary conditions (essential) according to exact solution
    % ux = (1-nu)*traction/E*x
    % uy = (1-nu)*traction/E*y
    % -----------------------------------------------------------------
        
        BC.UU = @(x) (1-nu)*traction/E*x(:,1);
        BC.VV = @(x) (1-nu)*traction/E*x(:,2);
        
        % column vector of prescribed displacement dof  
        BC.fix_disp_dof = 1:Mesh.nDOF;

        % prescribed displacement for each dof [u1; u2; ...] [m]
        BC.fix_disp_value = zeros(length(BC.fix_disp_dof),1);  
        BC.fix_disp_value(1:2:end) = BC.UU(Mesh.x);
        BC.fix_disp_value(2:2:end) = BC.VV(Mesh.x);  
        BC.fix_disp_value = @(t) BC.fix_disp_value;
        
    %% Neumann BC
    % -----------------------------------------------------------------
        % column vector of prescribed traction dofs
        BC.traction_force_dof = [];

        % magnitude of prescribed tractions [N]
        BC.traction_force_dof_value = [];

        % NOTE: this is slower than prescribing tractions at dofs
        % column vector of prescribed traction nodes 
        BC.traction_force_node = Mesh.right_nodes;  

        % prescribed traction [t1x t1y;t2x t2y;...] [N]
        Fnode = 1/(length(BC.traction_force_node) - 1);
        BC.traction_force_value = Fnode*[zeros(size(BC.traction_force_node)), zeros(size(BC.traction_force_node))];
        
        % find the nodes in the top right and bottom right corners
        toprightnode = find(Mesh.x(BC.traction_force_node,2) == max(Mesh.x(:,2)));
        botrightnode = find(Mesh.x(BC.traction_force_node,2) == min(Mesh.x(:,2)));
        
        BC.traction_force_value(toprightnode,1) = BC.traction_force_value(toprightnode,1)/2;
        BC.traction_force_value(botrightnode,1) = BC.traction_force_value(botrightnode,1)/2;
        
        % Empty function for application of tractions using edge elements
        BC.c_N_t_f = @(x,t)[];
          
        % NOTE: point loads at any of the element nodes can also be 
        % added as a traction.

        % magnitude of distributed body force [N/m] [bx;by]
            % 1D: [N/m], 2D: [N/m2]
        	% NOTE: if no body force, use '@(x)[]'
         	% NOTE: anonymous functions is defined with respect to the 
            %      variable x,  which is a vector [x(1) x(2)] = [x y]
        BC.b = @(x)[];    

%% Computation controls

        % quadrature order
        Control.qo = quadorder;

        % Nodal averaging for discontinuous variables (stress/strain)
        % 'none', 'nodal', 'center'
        Control.stress_calc = 'nodal';

        % penalty parameter for solution of static problem with 
        % LinearSolver3
        Control.beta = 10^10;

        % method used for solving linear problem:
        % 'LinearSolver1': Partitioning
        % 'LinearSolver2': Zeroing DOFs in stiffness matrix 
        %                   corresponding to essential boundaries
        % 'LinearSolver3': Penalty method
        Control.LinearSolver = 'LinearSolver1'; 
 
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
