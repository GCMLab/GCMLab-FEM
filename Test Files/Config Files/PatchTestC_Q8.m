function [Mesh, Material, BC, Control] = PatchTestC_Q8(config_dir, progress_on)
    global E nu traction quadorder meshfilename
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

    % Dirichlet boundary conditions (essential) according to exact solution
    % ux = (1-nu)*t/E*x
    % uy = (1-nu)*t/E*y
    % -----------------------------------------------------------------
        
        BC.UU = @(x) (1-nu)*traction/E*x(:,1);
        BC.VV = @(x) (1-nu)*traction/E*x(:,2);
        
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
            Fright = traction*max(Mesh.x(:,2))/(length(Mesh.right_nodes) - 1); % traction * 1 element length (assumed evenly distributed)
            Ftop   = traction*max(Mesh.x(:,1))/(length(Mesh.top_nodes)   - 1); % traction * 1 element length (assumed evenly distributed)
            BC.traction_force_value =       [   Fright*ones(size(Mesh.right_nodes(index_right))),     zeros(size(Mesh.right_nodes(index_right)));      % right side nodes
                                                zeros(size(Mesh.top_nodes(index_top))),               Ftop*ones(size(Mesh.top_nodes(index_top)));           % top side nodes
                                                Fright*1/2,                                           Ftop*1/2                                           ]; % top right node

            % find the nodes in the top left and bottom right corners
            botrightnode = find(Mesh.x(BC.traction_force_node,2) == min(Mesh.x(:,2)));
            topleftnode  = find(Mesh.x(BC.traction_force_node,1) == min(Mesh.x(:,1)));

            BC.traction_force_value(botrightnode,1) = BC.traction_force_value(botrightnode,1)/2;
            BC.traction_force_value(topleftnode,2) = BC.traction_force_value(topleftnode,2)/2;
        elseif strcmp(Mesh.type, 'Q9') || strcmp(Mesh.type, 'T6') || strcmp(Mesh.type, 'Q8')
            Fright = traction*max(Mesh.x(:,2))/((length(Mesh.right_nodes) - 1)/2);
            Ftop   = traction*max(Mesh.x(:,1))/((length(Mesh.top_nodes)   - 1)/2);
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

        % Nodal averaging for discontinuous variables (stress/strain)
        % 'none', 'nodal', 'center'
        Control.stress_calc = 'nodal';

        % penalty parameter for solution of static problem with 
        % LinearSolver3
        Control.beta = 10^10;

        % method used for solving linear problem:
            % 'LinearSolver1'
            % 'LinearSolver2'
            % 'LinearSolver3'
        Control.LinearSolver = 'LinearSolver1';
 
end