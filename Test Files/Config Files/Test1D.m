function [Mesh, Material, BC, Control] = Test1D(config_dir, progress_on)
%TEST1D Mesh, material parameters, boundary conditions, 
%and control parameters
%   Mesh = MASTERCONFIGFILE() is a structure array with the
%   following fields: 
%       .type:          the topological class of finite element; it is in 
%                       the general form 'topology-#of nodes' ie a three 
%                       node triangle is T3 a four node quadralateral is 
%                       Q4 a 4 node tetrahedra is H4 a 27 node brick is 
%                       B27 etc. Presently defined are L2, Q4, and Q9. 
%       .nsd:           Number of spatial dimensions
%       .ne:            Total number of elements in the mesh
%       .nne:           Vector of number of nodes per element (size nsd x 1)
%       .nn:            Total number of nodes 
%       .nDOFe:         Number of DOFs per element
%       .nDOF:          Total number of DOFs
%       .x:             Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)
%       .conn:          Array of element connectivity (size ne x nne)
%       .eneighbours:   Array of neighbouring elements (size ne x nneighbours
%                       in which nneighbours is 1 for 1D elements and 4
%                       for 2D elements)
%       .DOF:           Array of DOF indices (size nn x nsd)
%       .nodeconn:      Array of nodal connectivity (size nn x 8)
%                       containing the indices of elements connected to 
%                       each node
%       .left_nodes     Nodes on the left edge of the domain
%       .left_dof       DOFs on the left edge of the domain
%       .right_nodes    Nodes on the right edge of the domain
%       .right_dof      DOFs on the right edge of the domain
%       .xdofs          DOFs in the x-direction
%       .ydofs          DOFs in the y-direction
%       .zdofs          DOFs in the z-direction
%   Two-dimensional meshes also contain the fields,
%       .top_nodes      Nodes on the top edge of the domain
%       .top_dof        DOFs on the top edge of the domain
%       .top_dofx       DOFs on the top boundary in the x-direction
%       .top_dofy       DOFs on the top boundary in the y-direction
%       .bottom_nodes   Nodes on the bottom edge of the domain
%       .bottom_dof     DOFs on the bottom edge of the domain
%       .bottom_dofx    DOFs on the bottom boundary in the x-direction
%       .bottom_dofy    DOFs on the bottom boundary in the y-direction
%       .left_dofx      DOFs on the left boundary in the x-direction
%       .left_dofy      DOFs on the left boundary in the y-direction
%       .right_dofx     DOFs on the right boundary in the x-direction
%       .right_dofy     DOFs on the right boundary in the y-direction
%   Three-dimensional meshes also contain the fields, 
%       .near_nodes     nodes on the nearest face of the domain
%       .near_dof       DOFs on the nearest face of the domain
%       .near_dofx      DOFs on the near face in the x-direction
%       .near_dofy      DOFs on the near face in the y-direction
%       .near_dofz      DOFs on the near face in the z-direction
%       .far_nodes      Nodes on the farthest face of the domain
%       .far_dof        DOFs on the farthest face of the domain
%       .far_dofx       DOFs on the far face in the x-direction
%       .far_dofy       DOFs on the far face in the y-direction
%       .far_dofz       DOFs on the far face in the z-direction
%       .left_dofz      DOFs on the left face in the z-direction
%       .right_dofz     DOFs on the right face in the z-direction
%       .top_dofz       DOFs on the top face in the z-direction
%       .bottom_dofz    DOFs on the bottom face in the z-direction
%       
%   Mesh = TEST1D(config_dir) defines the mesh using GMSH file 
%   import located in the directory config_dir
%
%   [Mesh, Material] = TEST1D() also returns a
%   structure array with the following fields: 
%       .nmp:           number of material properties
%       .Prop:          Material properties
%       .Prop.E:        Modulus of elasticity
%       .Prop.nu:       Poisson's ratio
%       .Prop.Dtype:    2D approximation ('PlaneStrain' or 'PlainStress')
%       .Prop.t:        Material thickness
% 
%   [Mesh, Material, BC] = TEST1D() also returns a structure
%   array with the following fields: 
%       .fix_disp_dof:              Column vector of degrees of freedom 
%                                   with prescribed displacements
%                                   (size nfixed x 1)
%       .fix_disp_value             Column vector of prescribed 
%                                   displacements (size nfixed x 1)
%       .traction_force_dof         Column vector of degrees of freedom
%                                   with prescribed tractions
%       .traction_force_dof_value   Column vector of prescribed tractions
%                                   on DOF
%       .traction_force_node        Column vector of nodes with 
%                                   prescribed tractions
%       .traction_force_value       Column vector of prescribed tractions
%                                   on nodes
%       .b                          Anonymous function of distributed
%                                   body force (size 1 x nsd)
% 
%   [Mesh, Material, BC, Control] = TEST1D() also returns a 
%   structure array with the following fields: 
%       .qo:            Quadrature order
%       .MagCoef:       Displacement magnification coefficient
%                       for visualization
%       .contour:       Nodal averaging for discontinous variables
%                       ('none', 'nodal')
%       .beta:          Penalty parameter  
%       .LinearSolver   Method used for solving linear problem:
%                       'LinearSolver1': Partitioning
%                       'LinearSolver2': Zeroing DOFs in stiffness matrix 
%                                        corresponding to essential boundaries
%                       'LinearSolver3': Penalty method
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   config_dir:     (OPTIONAL) File path for the directory where 
%                   unstructured mesh is stored

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
        % NLTH1 - Nonlinear thermal transient
    Material.Model = 'LE1';

    % number of material properties
    Material.nmp = 1;
        
    % Properties material 1
    Material.Prop(1).E0 = 2e11; % Young's modulus [Pa]
    Material.Prop(1).nu = 0; % Poisson's ratio

    % Constitutive law: 'PlaneStrain' or 'PlaneStress' 
    Material.Dtype = ''; 

    % Thickness (set as default to 1)
    % 1D: [m2], 2D: [m]
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
    MeshType = 'MANUAL';        
    
    switch MeshType
        case 'MANUAL'
            % location of initial node [m] [x0;y0;z0] 
            x1 = 0;
            % number of space dimensions 
            nsd = 1;
            % size of domain [m] [Lx;Ly;Lz] 
            L = 3.4;
            % number of elements in each direction [nex; ney; nez] 
            nex = 52;
            % element type ('L2')
            type = 'L2';
            
            Mesh = BuildMesh_structured(nsd, x1, L, nex, type, progress_on, Material.ProblemType);
        case 'IMPORTED'
            % Allows input of files from GMSH
            % Note: the only currently supported .msh file formatting is
            % Version 2 ASCII
            % Ctrl + e to export the mesh, specify extension .msh, specify
            % format Version 2 ASCII
            meshFileName = 'Mesh Files\Test1D.msh';
            % number of space dimensions
            nsd = 1;
            
            Mesh = BuildMesh_imported(meshFileName, nsd, config_dir, progress_on, 0, Material.ProblemType);
        case 'EXCEL'
            meshFileName = 'CircularInclusion.xlsx';
            % number of space dimensions
            nsd = 2;
            
            Mesh = BuildMesh_EXCEL(meshFileName, nsd, config_dir, progress_on, Material.ProblemType);
    end

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

    % Dirichlet boundary conditions (essential)
    % -----------------------------------------------------------------
        % column vector of prescribed displacement dof  
        BC.fix_disp_dof = Mesh.left_dof;

        % prescribed displacement for each dof [u1; u2; ...] [m]
        BC.fix_disp_value = @(t) zeros(length(BC.fix_disp_dof),1);

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
        % tensile force applied to right edge
        BC.traction = 6e9; % applied traction [N]
        BC.traction_force_value = @(t) BC.traction;
            
        % NOTE: point loads at any of the element nodes can also be 
        % added as a traction.

        % magnitude of distributed body force [N/m] [bx;by]
            % 1D: [N/m], 2D: [N/m2]
        	% NOTE: if no body force, use '@(x)[]'
         	% NOTE: anonymous functions is defined with respect to the 
            %      variable x,  which is a vector [x(1) x(2)] = [x y]
        BC.b_force = 2e9; % applied body force [N/m]
        BC.b = @(x,t) BC.b_force;    

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
