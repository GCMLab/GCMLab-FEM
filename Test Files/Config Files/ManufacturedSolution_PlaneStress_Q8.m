function [Mesh, Material, BC, Control] = ManufacturedSolution_PlaneStress_Q8(config_dir, progress_on)
global meshfilename quadorder E nu

%% Mesh Properties
    if progress_on
        disp([num2str(toc),': Building Mesh...']);
    end
    
    % Mesh formats: 
    %   'MANUAL'    - In-house structured meshing
    % 	'GMSH'      - Import .msh file from GMSH, structured or unstructured
    MeshType = 'GMSH';        

    
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
            
            Mesh = BuildMesh_structured(nsd, x1, L, nex, type, progress_on);
        case 'GMSH'
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
            
            Mesh = BuildMesh_GMSH(meshFileName, nsd, config_dir, progress_on, Q8_reduced);              
    end    
    

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

    % Specify stiffness matrix and stress/strain calculation files
    Material.ConstitutiveLawFile = 'getD';
    Material.StiffnessMatrixFile = 'getK_elastic';
    Material.StressStrainFile = 'getStrain';

    % number of material properties
    Material.nmp = 1;
    
    % Properties material 1
    Material.Prop(1).E = E; % Young's modulus [Pa]
    Material.Prop(1).nu = nu; % Poisson's ratio

    % type of material per element
    Mesh.MatList = zeros(Mesh.ne, 1, 'int8');
    
    % assign material type to elements
    Mesh.MatList(:) = 1;

    % Constitutive law: 'PlaneStrain' or 'PlaneStress' 
    Material.Dtype = 'PlaneStress'; 

    % Thickness (set as default to 1)
    Material.t = @(x) 1;

    % Alternatively, import a material file
    % Material = Material_shale();

%% Boundary Conditions
    % {TIPS}------------------------------------------------------------
        % TIP selecting edges:
        % bottom_nodes = find(Mesh.x(:,2)==0); 
        % top_nodes = find(Mesh.x(:,2)==2);
        % left_nodes = find(Mesh.x(:,1)==0);
        % right_nodes = find(Mesh.x(:,1)==4);
        % bottom_dof = [bottom_nodes*2 - 1; bottom_nodes*2];
        % top_dof = [top_nodes*2 - 1;top_nodes*2];

    % Dirichlet boundary conditions (essential) according to manufactured
    % solution
    % ux = x^5 + x*y^3 - y^6
    % uy = x^5 + x*y^3 - y^6
    % -----------------------------------------------------------------
        BC.UU = @(x) x(:,1).^5 + x(:,1).*x(:,2).^3 - x(:,2).^6;
        BC.VV = @(x) x(:,1).^5 + x(:,1).*x(:,2).^3 - x(:,2).^6;

        % column vector of prescribed displacement dof  
        BC.fix_disp_dof1 = Mesh.left_dof;
        BC.fix_disp_dof2 = Mesh.right_dof;
        BC.fix_disp_dof3 = Mesh.bottom_dof;
        BC.fix_disp_dof4 = Mesh.top_dof;
        BC.fix_disp_dof = unique([BC.fix_disp_dof1;BC.fix_disp_dof2;BC.fix_disp_dof3;BC.fix_disp_dof4]);

        % prescribed displacement for each dof [u1; u2; ...] [m]
        BC.fix_disp_value = zeros(length(BC.fix_disp_dof),1);
        BC.fix_disp_value(1:2:end) = BC.UU([Mesh.x(BC.fix_disp_dof(2:2:end)/2,1),Mesh.x(BC.fix_disp_dof(2:2:end)/2,2)]);
        BC.fix_disp_value(2:2:end) = BC.VV([Mesh.x(BC.fix_disp_dof(2:2:end)/2,1),Mesh.x(BC.fix_disp_dof(2:2:end)/2,2)]);
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
    
        % NOTE: point loads at any of the element nodes can also be 
        % added as a traction.

        % magnitude of distributed body force [N/m] [bx;by] according to
        % the manufactured solution:
        % bx = -E / (1-v^2) * ( 20x^3 + 3vy^2          + (1-v)/2*[ 6xy - 30y^4 + 3y^2 ])
        % by = -E/  (1-v^2) * ( (1-v)/2*[3y^2 + 20x^3] +           3vy^2 + 6xy - 30y^4 )
            % 1D: [N/m], 2D: [N/m2]
        	% NOTE: if no body force, use '@(x)[]'
         	% NOTE: anonymous functions is defined with respect to the 
            %      variable x,  which is a vector [x(1) x(2)] = [x y]
        BC.b = @(x,t)[-E / (1-nu^2)  * ( 20*x(1).^3 + 3*nu*x(2).^2              + (1-nu)/2*( 6*x(1).*x(2) - 30*x(2).^4 + 3*x(2).^2));
                    -E / (1-nu^2)  * ( (1-nu)/2*( 3*x(2).^2  + 20*x(1).^3)    + 3*nu*x(2).^2 + 6*x(1).*x(2) - 30*x(2).^4 )];

%% Computation controls

        % quadrature order
        Control.qo = quadorder;

        % Nodal averaging for discontinuous variables (stress/strain)
        % 'none', 'nodal'
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
 
        % Newton Raphson controls
        Control.r_tol = 1e-5; % Tolerance on residual forces
        Control.iter_max = 50; % Maximum number of iteration in Newton Raphson algorithm
        
end
