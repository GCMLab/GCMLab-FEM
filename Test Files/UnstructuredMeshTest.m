function [Mesh, Material, BC, Control] = UnstructuredMeshTest(config_dir, progress_on)

%% Mesh Properties
    if progress_on
        disp([num2str(toc),': Building Mesh...']);
    end
    
    % Mesh format - 'STRUCTURED', 'UNSTRUCTURED'
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
            
            Mesh = BuildMesh_structured(nsd, x1, L, nex, type);
        case 'UNSTRUCTURED'
        case 'GMSH'
            % Allows input of files from GMSH
            % Note: the only currently supported .msh file formatting is
            % Version 2 ASCII
            % Ctrl + e to export the mesh, specify extension .msh, specify
            % format Version 2 ASCII
            meshFileName = 'Unstructured_sample.msh';
            % number of space dimensions 
            nsd = 2;
          
            Mesh = BuildMesh_GMSH(meshFileName, nsd, config_dir, progress_on);            
    end    
    

%% Material Properties (Solid)

    % NOTES-------------------------------------------------------------
                                
        % NOTE: anonymous functions are defined with respect to the variable x,
        % which is a vector [x(1) x(2) x(3)] = [x y z]

        % NOTE: Material properties must be continuous along an element, 
        % otherwise, quadrature order must be increased significantly

    % Young's modulus [Pa]
    Material.E = @(x) 4;  

    % Constitutive law: 'PlaneStrain' or 'PlaneStress' 
    Material.Dtype = 'PlaneStrain'; 

    % Thickness (set as default to 1)
    Material.t = @(x) 1;

    % Poisson's ratio (set as default to 0.3)
    Material.nu = @(x) 0;

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

    % Dirichlet boundary conditions (essential)
    % -----------------------------------------------------------------
        % column vector of prescribed displacement dof  
        BC.fix_disp_dof = [Mesh.left_dof];

        % prescribed displacement for each dof [u1; u2; ...] [m]
        BC.fix_disp_value = zeros(length(BC.fix_disp_dof),1);  

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
        BC.traction_force_value = Fnode*[ones(size(BC.traction_force_node)), zeros(size(BC.traction_force_node))];
        
        % find the nodes in the top right and bottom right corners
        toprightnode = find(Mesh.x(BC.traction_force_node,2) == max(Mesh.x(:,2)));
        botrightnode = find(Mesh.x(BC.traction_force_node,2) == min(Mesh.x(:,2)));
        
        BC.traction_force_value(toprightnode,1) = BC.traction_force_value(toprightnode,1)/2;
        BC.traction_force_value(botrightnode,1) = BC.traction_force_value(botrightnode,1)/2;
    
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
        Control.qo = 2;

        % displacement magnification coefficient (for visualization)
        Control.MagCoef = 1;

        % Nodal averaging for discontinuous variables (stress/strain)
        % 'none', 'nodal', 'center'
        Control.contour = 'center';

        % penalty parameter for solution of static problem with 
        % LinearSolver3
        Control.beta = 10^10;

        % method used for solving linear problem:
            % 'LinearSolver1'
            % 'LinearSolver2'
            % 'LinearSolver3'
        Control.LinearSolver = 'LinearSolver1';
 
end