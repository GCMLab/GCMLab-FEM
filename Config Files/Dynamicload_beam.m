function [Mesh, Material, BC, Control] = Dynamicload_beam(config_dir, progress_on)

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
        % ST1 - Stiffening model with 1st invariant of strain
        % LED1 - Linear elasticity dynamic
    Material.Model = 'LED1';
        
    % number of material properties
    Material.nmp = 1;

    % Properties material 1
    Material.Prop(1).E0 = 2e11; % Young's modulus [Pa]
    Material.Prop(1).nu = 0.3; % Poisson's ratio
    Material.Prop(1).C = 0; % Damping Coefficient
    Material.Prop(1).rho = 2400; % Poisson's ratio
    

    % Constitutive law: 'PlaneStrain' or 'PlaneStress' 
    Material.Dtype = 'PlaneStress'; 
    

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
    % 	'GMSH'  - Import .msh file from GMSH, structured or unstructured
    %   'EXCEL' - Import .xlsx file, structured or unstructured
    MeshType = 'MANUAL';        
    
    switch MeshType
        case 'MANUAL'
            % location of initial node [m] [x0;y0;z0] 
            x1 = [0;0;0];
            % number of space dimensions 
            nsd = 2;
            % size of domain [m] [Lx;Ly;Lz] 
            L = [5;1];
            % number of elements in each direction [nex; ney; nez] 
            nex = [5;1]*10;
            % element type ('Q4')
            type = 'Q4';
            
            Mesh = BuildMesh_structured(nsd, x1, L, nex, type, progress_on);
        case 'GMSH'
            % Allows input of files from GMSH
            % Note: the only currently supported .msh file formatting is
            % Version 2 ASCII
            % Ctrl + e to export the mesh, specify extension .msh, specify
            % format Version 2 ASCII
            meshFileName = '2DBarMesh.msh';
            % number of space dimensions 
            nsd = 2;
            
            Mesh = BuildMesh_GMSH(meshFileName, nsd, config_dir, progress_on); 
        case 'EXCEL'
            meshFileName = 'CricularInclusion.xlsx';
            % number of space dimensions
            nsd = 2;
            
            Mesh = BuildMesh_EXCEL(meshFileName, nsd, config_dir, progress_on);
    end    


    
%% Assign materials to mesh
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
        
    % Manufactured solution
    % u := (x1, x2, t) -> -1/1000*sin(1/2*pi*x1)*sin(1/2*pi*x2)*sin(2*pi*t)
    % v := (x1, x2, t) -> 1/1000*cos(1/2*pi*x1)*cos(1/2*pi*x2)*cos(2*pi*t)

    % Dirichlet boundary conditions (essential)
    % -----------------------------------------------------------------
        % prescribed displacement for each dof [u1; u2; ...] [m]
        % u := (0, 1, t) -> sin(pi*t/2);
        % v := (0, 1, t) -> 0

        BC.fix_disp_dof = [Mesh.left_dofy(end); Mesh.right_dofx;  Mesh.right_dofy];
        
        % Prescribed dispalcement
        
%         BC.fix_disp_value = @(t) [1*sin(pi*t/2)*(t < 8); zeros(length(BC.fix_disp_dof(2:end)),1)];
        BC.fix_disp_value = @(t) [t*(t < 1); zeros(length(BC.fix_disp_dof(2:end)),1)];

    %% Neumann BC
    % -----------------------------------------------------------------

        % Magnitude of amplitude of Fext applied to nodes at free-end of
        % beam
        BC.Fn = [];

        % column vector of prescribed traction dofs
        BC.traction_force_dof = [];

        % magnitude of prescribed tractions [N]
        BC.traction_force_dof_value = [];

        % NOTE: this is slower than prescribing tractions at dofs
        % column vector of prescribed traction nodes 
        BC.traction_force_node = [];  

        % prescribed traction [t1x t1y;t2x t2y;...] [N]
        BC.traction_force_value = @(t) [];
    
        % NOTE: point loads at any of the element nodes can also be 
        % added as a traction.

        % magnitude of distributed body force [N/m] [bx;by]
            % 1D: [N/m], 2D: [N/m2]
        	% NOTE: if no body force, use '@(x)[]'
         	% NOTE: anonymous functions is defined with respect to the 
            %      variable x,  which is a vector [x(1) x(2)] = [x y]
        BC.b = @(x,t) [ ];

%% Initial Conditions
        BC.IC_temp = zeros(Mesh.nDOF,1);
        
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
%         Control.EndTime   = 32; 
        Control.EndTime   = 4; 
        NumberOfSteps     = 100;
        Control.TimeStep  = (Control.EndTime - Control.StartTime)/(NumberOfSteps);
        Control.dSave     = 1;
        
        % transient controls
        Control.TimeCase = 'dynamic';    
                        % Static → Control.TimeCase = 'static;
                        % Transient → Control.TimeCase = 'transient';
                        % Dynamic (HHT method)→ Control.TimeCase = 'dynamic';
        Control.alpha = 0; %
        %   If Control.Timecase = 'transient'
        %           α = 1 Backward Euler, α = 1/2 Crank-Nicolson
        %   If Control.Timecase = 'dynamic'
        %           α [-1/3, 0]

        % Newton Raphson controls
        Control.r_tol = 1e-7; % Tolerance on residual forces
        Control.iter_max = 50; % Maximum number of iteration in Newton Raphson algorithm
        
        
end
