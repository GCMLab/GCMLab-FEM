function [Mesh, Material, BC, Control] = PatchTestC_Q8(config_dir, progress_on)
    global meshfilename quadorder
%PATCHTESTC_Q8 Mesh, material parameters, boundary conditions, 
%and control parameters
%   Mesh = PATCHTESTC_Q8() is a structure array with the
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
%   Mesh = PATCHTESTC_Q8(config_dir) defines the mesh using GMSH file 
%   import located in the directory config_dir
%
%   [Mesh, Material] = PATCHTESTC_Q8() also returns a
%   structure array with the following fields: 
%       .nmp:           number of material properties
%       .Prop:          Material properties
%       .Prop.E0:       Modulus of elasticity
%       .Prop.nu:       Poisson's ratio
%       .Prop.Dtype:    2D approximation ('PlaneStrain' or 'PlainStress')
%       .Prop.t:        Material thickness
% 
%   [Mesh, Material, BC] = PATCHTESTC_Q8() also returns a structure
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
%   [Mesh, Material, BC, Control] = PATCHTESTC_Q8() also returns a 
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
        % NLTH1 - Nonlinear thermal transient
    Material.Model = 'LE1';

    % number of material properties
    Material.nmp = 1;
    
    % Properties material 1
    Material.Prop(1).E0 = 2540; % Young's modulus [Pa]
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
            % Optional 5th input in case Q8 with reduced integration is desired
            Q8_reduced = 'Q8'; %Do not consider this input if a case different than Q8 with reduced integration is desired
            
            Mesh = BuildMesh_imported(meshFileName, nsd, config_dir, progress_on, Q8_reduced, Material.ProblemType);
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

    % Dirichlet boundary conditions (essential) according to exact solution
    % ux = (1-nu)*t/E*x
    % uy = (1-nu)*t/E*y
    % -----------------------------------------------------------------
        
        BC.traction = 3.495; % applied traction (both directions)
    
        BC.UU = @(x) (1-Material.Prop(1).nu)*BC.traction/Material.Prop(1).E0*x(:,1);
        BC.VV = @(x) (1-Material.Prop(1).nu)*BC.traction/Material.Prop(1).E0*x(:,2);
        
        % column vector of prescribed displacement dof  
        topleftnode = Mesh.left_nodes(find(Mesh.x(Mesh.left_nodes,2) == max(Mesh.x(:,2))));
        botleftnode = Mesh.left_nodes(find(Mesh.x(Mesh.left_nodes,2) == min(Mesh.x(:,2))));
        
        BC.fix_disp_dof1 = [Mesh.left_nodes*2-1];
        BC.fix_disp_dof2 = Mesh.bottom_nodes*2;
        
        BC.fix_disp_dof = [BC.fix_disp_dof1;BC.fix_disp_dof2];

        % prescribed displacement for each dof [u1; u2; ...] [m]
        BC.fix_disp_value = zeros(length(BC.fix_disp_dof),1);
        BC.fix_disp_value1 = BC.UU([Mesh.x(Mesh.left_nodes,1),Mesh.x(Mesh.left_nodes,2)]);
        BC.fix_disp_value2 = BC.VV([Mesh.x(Mesh.bottom_nodes,1),Mesh.x(Mesh.bottom_nodes,2)]);
        BC.fix_disp_value = @(t)  [BC.fix_disp_value1;BC.fix_disp_value2];  

    %% Neumann BC
    % -----------------------------------------------------------------
        % column vector of prescribed traction dofs
        BC.traction_force_dof = [];

        % magnitude of prescribed tractions [N]
        BC.traction_force_dof_value = [];

        % NOTE: this is slower than prescribing tractions at dofs
        % column vector of prescribed traction nodes 
        toprightnode = Mesh.right_nodes(Mesh.x(Mesh.right_nodes,2) == max(Mesh.x(:,2)));
        index_right = Mesh.right_nodes ~= toprightnode;
        index_top   = Mesh.top_nodes   ~= toprightnode;
        
        BC.traction_force_node = [Mesh.right_nodes(index_right);  Mesh.top_nodes(index_top); toprightnode];
        
        % prescribed traction [t1x t1y;t2x t2y;...] [N]
        %t = 4;
        if strcmp(Mesh.type, 'Q4') || strcmp(Mesh.type, 'T3')
            Fright = BC.traction*max(Mesh.x(:,2))/(length(Mesh.right_nodes) - 1); % traction * 1 element length (assumed evenly distributed)
            Ftop   = BC.traction*max(Mesh.x(:,1))/(length(Mesh.top_nodes)   - 1); % traction * 1 element length (assumed evenly distributed)
            BC.traction_force_value =       [   Fright*ones(size(Mesh.right_nodes(index_right))),     zeros(size(Mesh.right_nodes(index_right)));      % right side nodes
                                                zeros(size(Mesh.top_nodes(index_top))),               Ftop*ones(size(Mesh.top_nodes(index_top)));           % top side nodes
                                                Fright*1/2,                                           Ftop*1/2                                           ]; % top right node

            % find the nodes in the top left and bottom right corners
            botrightnode = find(Mesh.x(BC.traction_force_node,2) == min(Mesh.x(:,2)));
            topleftnode  = find(Mesh.x(BC.traction_force_node,1) == min(Mesh.x(:,1)));

            BC.traction_force_value(botrightnode,1) = BC.traction_force_value(botrightnode,1)/2;
            BC.traction_force_value(topleftnode,2) = BC.traction_force_value(topleftnode,2)/2;
        elseif strcmp(Mesh.type, 'Q9') || strcmp(Mesh.type, 'T6') || strcmp(Mesh.type, 'Q8')
            Fright = BC.traction*max(Mesh.x(:,2))/((length(Mesh.right_nodes) - 1)/2);
            Ftop   = BC.traction*max(Mesh.x(:,1))/((length(Mesh.top_nodes)   - 1)/2);
            BC.traction_force_value = zeros(length(BC.traction_force_node),2);
            for n = 1:length(BC.traction_force_node)
                if n <= length(Mesh.right_nodes(index_right)) % then node is on the right edge
                    if any( BC.traction_force_node(n) == Mesh.conn(:,1:4),'all') % then node is a corner node
                        BC.traction_force_value(n,:) = [Fright/3, 0];
                    else % then node is a midside node
                        BC.traction_force_value(n,:) = [Fright*2/3,0];
                    end
                elseif n == length(BC.traction_force_node) % then node is the top right node
                        BC.traction_force_value(n,:) = [Fright/6, Ftop/6];
                else % then node is on the top edge
                    if any( BC.traction_force_node(n) == Mesh.conn(:,1:4),'all') % then node is a corner node
                        BC.traction_force_value(n,:) = [0, Ftop/3];
                    else % then node is a midside node
                        BC.traction_force_value(n,:) = [0, Ftop*2/3];
                    end
                end
            end
            
            % find the nodes in the top left and bottom right corners
            botrightnode = find(Mesh.x(BC.traction_force_node,2) == min(Mesh.x(:,2)));
            topleftnode  = find(Mesh.x(BC.traction_force_node,1) == min(Mesh.x(:,1)));

            BC.traction_force_value(botrightnode,1) = BC.traction_force_value(botrightnode,1)/2;
            BC.traction_force_value(topleftnode,2) = BC.traction_force_value(topleftnode,2)/2;
            
        else
            error('Unsupported element type')
        end
  
        
    
        % NOTE: point loads at any of the element nodes can also be 
        % added as a traction.

        % magnitude of distributed body force [N/m] [bx;by]
            % 1D: [N/m], 2D: [N/m2]
        	% NOTE: if no body force, use '@(x)[]'
         	% NOTE: anonymous functions is defined with respect to the 
            %      variable x,  which is a vector [x(1) x(2)] = [x y]
        BC.b = @(x,t)[];    

%% Computation controls

        % quadrature order
        Control.qo = quadorder;

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
