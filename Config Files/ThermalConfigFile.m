function [Mesh, Material, BC, Control] = ThermalConfigFile(config_dir, progress_on)
%THERMALCONFIGFILE Mesh, material parameters, boundary conditions, 
%and control parameters
%   Mesh = THERMALCONFIGFILE() is a structure array with the
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
%   Mesh = MASTERCONFIGFILE(config_dir) defines the mesh using GMSH file 
%   import located in the directory config_dir
%
%   [Mesh, Material] = MASTERCONFIGFILE() also returns a
%   structure array with the following fields: 
%       .nmp:           number of material properties
%       .Prop:          Material properties
%       .Prop.E:        Modulus of elasticity
%       .Prop.nu:       Poisson's ratio
%       .Prop.Dtype:    2D approximation ('PlaneStrain' or 'PlainStress')
%       .Prop.t:        Material thickness
% 
%   [Mesh, Material, BC] = MASTERCONFIGFILE() also returns a structure
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
%   [Mesh, Material, BC, Control] = MASTERCONFIGFILE() also returns a 
%   structure array with the following fields: 
%       .qo:            Quadrature order
%       .stress_calc    Calculation of values for discontinous variables
%                       ('none', 'nodal', 'center', 'L2projection')
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
    Material.Model = 'TH1';
    
    % number of material properties
    Material.nmp = 1;

    % Properties material 1
    Material.Prop(1).k1 = 50; % Conductivity in the x-direction [W/mK]
    Material.Prop(1).k2 = 25; % Conductivity in the y-direction [W/mK]
    Material.Prop(1).C = 5000;   % Heat Capacity (specific heat * density) = [J/K kg] * [kg/m^3] = [J/K m^3]
   
    % Constitutive law: 'ISO' or 'ORTHO'
    Material.Dtype = 'ORTHO'; 
    
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
            x1 = [0;0;0];
            % number of space dimensions 
            nsd = 2;
            % size of domain [m] [Lx;Ly;Lz] 
            L = [1;1];
            % number of elements in each direction [nex; ney; nez] 
            nex = [2;2]*10;
            % element type ('Q4')
            type = 'Q4';
            
            Mesh = BuildMesh_structured(nsd, x1, L, nex, type, progress_on, Material.ProblemType);
        case 'IMPORTED'
            % Allows input of files from GMSH
            % Note: the only currently supported .msh file formatting is
            % Version 2 ASCII
            % Ctrl + e to export the mesh, specify extension .msh, specify
            % format Version 2 ASCII
            meshFileName = 'Unstructured_sample.msh';
            % number of space dimensions 
            nsd = 2;
            % Optional 5th input in case Q8 with reduced integration is desired
            Q8_reduced = 'Q8'; %Do not consider this input if a case different than Q8 with reduced integration is desired
            
            Mesh = BuildMesh_imported(meshFileName, nsd, config_dir, progress_on, 0, Material.ProblemType);            
%             Mesh = BuildMesh_imported(meshFileName, nsd, config_dir, progress_on,Q8_reduced);  
        case 'EXCEL'
            meshFileName = 'CricularInclusion.xlsx';
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
        % TIP selecting edges - GMESH file or MANUAL mesh:
        % bottom_nodes = find(Mesh.x(:,2)==0); 
        % top_nodes = find(Mesh.x(:,2)==2);
        % left_nodes = find(Mesh.x(:,1)==0);
        % right_nodes = find(Mesh.x(:,1)==4);
        % bottom_dof = [bottom_nodes*2 - 1; bottom_nodes*2];
        % top_dof = [top_nodes*2 - 1;top_nodes*2];
        
        % TIP for .fem file - HYPERMESH
        %   Fixed BC
        %   temp = Mesh.DOF(Mesh.BC_nE,:).*Mesh.BC_E;
        %   BC.fix_disp_dof = nonzeros(reshape(temp, length(temp)*Mesh.nsd,1));
        %
        %   temp = Mesh.BC_nN_n;
        %   BC.traction_force_node = temp;

    % Dirichlet boundary conditions (essential)
    % -----------------------------------------------------------------
        % column vector of prescribed displacement dof  
        BC.fix_disp_dof = [Mesh.left_dof; Mesh.bottom_dof];

        % prescribed displacement for each dof [u1; u2; ...] [m]
        BC.fix_disp_value = @(t) [zeros(length(Mesh.left_dof),1); 25*ones(length(Mesh.bottom_dof),1)];  

    %% Neumann BC
    % -----------------------------------------------------------------
        %--------------------------------------------------------------
        % Method 1:
        %--------------------------------------------------------------
        % column vector of prescribed traction dofs
        BC.traction_force_dof = [];

        % magnitude of prescribed tractions [N]
        BC.traction_force_dof_value = [];
        
        %--------------------------------------------------------------
        % Method 2:
        %--------------------------------------------------------------
        % Edge elements and normals are obtained from Mesh structure.
        %
        % Tractions boundary conditions can be applied as:
        %   1 → function for vector on x and y (2 by 1 vector)
        %   2 → function for vector on n and s (normal
        %           and tangential directions of each point of the boundary)
        %           (2 by 1 vector)
        % These functions are stored in BC.c_N_t_f{i} (Neuman _ traction _ function)
        %
        % To recognize if the function is case 1 or 2, BC.c_N_t_flag must include
        % a vector with:
        %   0 → function for vector on x and y
        %   1 → function for vector on n and s
        %
        % To rectognize if the function is applied on the mechanical or diffusive DOF
        % BC.c_N_t_flag_case must include a vector with
        %   1 → if the set {i} is applied to the diffusive part
        %   2 → if the set {i} is applied to the mechanical part
        %
        % -------------------------------------------------------------
        % Note: The total number of traction functions must be compatible
        % with the number of sets defined for edge elements
        %
        % Note: BC.c_N_t_f, BC.c_N_t_flag, and BC.c_N_t_flag_case must have the same length 
        %
        % Note: If set 'i' is applied in the diffusive part, in BC.c_N_t_f{i} bust be
        % defined in normal components, and must be a single term, as
        % tangential components are not defined in a diffusive problem.
        % Moreover, BC.c_N_t_flag{i} must be set as 1, and
        % BC.c_N_t_flag_case{i} must be set as 1.
        % -------------------------------------------------------------
        BC.c_N_t_f = cell(1,1); % change to cell(1,number of sets) or
        % set as @(x,t) []  if no tractions of this type are applied
        BC.c_N_t_f{1} = @(x,t) []; % call BC.c_N_t_f{1.. number of stes} and define all sets
        BC.c_N_t_flag = []; %change to row vector if required. Example → [0,0,1]
        BC.c_N_t_flag_case = []; %change to row vector if required. Example → [1,1,2]

        %--------------------------------------------------------------
        % Method 3:
        %--------------------------------------------------------------
        % NOTE: this is slower than prescribing tractions at dofs
        % column vector of prescribed traction nodes 
        BC.traction_force_node = Mesh.right_nodes;  

        % prescribed traction [t1x t1y;t2x t2y;...] [N]
        Fnode = 0/(length(BC.traction_force_node) - 1);
        BC.traction_force_value = Fnode*[ones(size(BC.traction_force_node)), zeros(size(BC.traction_force_node))];
        
        % find the nodes in the top right and bottom right corners
        toprightnode = find(Mesh.x(BC.traction_force_node,2) == max(Mesh.x(:,2)));
        botrightnode = find(Mesh.x(BC.traction_force_node,2) == min(Mesh.x(:,2)));
        
        BC.traction_force_value(toprightnode,1) = BC.traction_force_value(toprightnode,1)/2;
        BC.traction_force_value(botrightnode,1) = BC.traction_force_value(botrightnode,1)/2;
        
        % Make the vector into an anonymous function in time
        BC.traction_force_value = @(t) BC.traction_force_value; 
    
        % NOTE: point loads at any of the element nodes can also be 
        % added as a traction.

        % magnitude of distributed body force [N/m] [bx;by]
            % 1D: [N/m], 2D: [N/m2]
        	% NOTE: if no body force, use '@(x)[]'
         	% NOTE: anonymous functions is defined with respect to the 
            %      variable x,  which is a vector [x(1) x(2)] = [x y]
        BC.b = @(x,t)[];    

%% Initial Conditions
        BC.IC = @(t) 50*ones(Mesh.nn,1);
        
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
        Control.stress_calc = 'L2projection';

        % penalty parameter for solution of static problem with 
        % LinearSolver3
        Control.beta = 10^10;
        
        % parallel inversion
        % Use parallel processing to invert the matrix.
        % Usually more efficient at 2e5 dofs
        Control.parallel = 1;

        % method used for solving linear problem:
        % 'LinearSolver1': Partitioning
        % 'LinearSolver2': Zeroing DOFs in stiffness matrix 
        %                   corresponding to essential boundaries
        % 'LinearSolver3': Penalty method
        Control.LinearSolver = 'LinearSolver1';
 
        % time controls
        Control.StartTime = 0;
        Control.EndTime   = 10;
        NumberOfSteps     = 100;
        Control.TimeStep  = (Control.EndTime - Control.StartTime)/(NumberOfSteps);
        % save displacements and stresses at each timestep in matlab 
        % debugging and testing purposes only, vtk files are otherwise
        % recommended
        Control.dSave     = 0; 
        % Plot load vs displacement curve (if not included, default is set
        % to 0)
        Control.plotLoadDispl = 0;
        % DOF to plot (only necessary if Control.plotLoadDispl = 1)
        Control.plotAt = 0; % [add DOF number]
        
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