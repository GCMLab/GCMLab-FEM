function [Mesh, Material, BC, Control] = ThermoElastic_Dirichlet(config_dir, progress_on)

global nex quadorder 

% nex:              number of elements in each direction [nex; ney; nez] 
% quadorder:        order of quadrature used in numerical integration

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
        % THLE1 - Thermoelasticity (Steady-State)
    Material.Model = 'THLE1';
    
    % number of material properties
    Material.nmp = 1;

    % Properties material 1
    Material.Prop(1).E0 = 2e11; % Young's modulus [Pa]
    Material.Prop(1).nu = 0.3; % Poisson's ratio
    Material.Prop(1).k1 = 50; % Conductivity in the x-direction [W/mK]
    Material.Prop(1).k2 = 25; % Conductivity in the y-direction [W/mK]
    Material.Prop(1).alpha = 140e-6; % Thermal expansion [1/K]
    Material.Prop(1).C = 5000;   % Heat Capacity (specific heat * density) = [J/K kg] * [kg/m^3] = [J/K m^3]

    % Constitutive law mechanic problem: 'PlaneStrain' or 'PlaneStress' 
    Material.Dtype_mech = 'PlaneStrain'; 

    % Constitutive law thermal problem: 'ISO' or 'ORTHO'
    Material.Dtype_therm = 'ISO'; 

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
    
    % location of initial node [m] [x0;y0;z0] 
    x1 = [0;0;0];
    % number of space dimensions 
    nsd = 2;
    % size of domain [m] [Lx;Ly;Lz] 
    L = [1;1];
    % element type ('Q4')
    type = 'Q4';

    Mesh = BuildMesh_structured(nsd, x1, L, nex, type, progress_on, Material.ProblemType);

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

    % Dirichlet boundary conditions (essential) - mechanical
    % -----------------------------------------------------------------
    
    % Manufacture solution
    % u:= (x1,x2) → -sin(pi*x1/2)*sin(pi*x2/2)
    % v:= (x1,x2) → cos(pi*x1/2)*cos(pi*x2/2)
        
        BC.UUx = @(x) - sin(pi.*x(:,1)./2).*sin(pi.*x(:,2)./2)./1000;
        BC.VVx = @(x)   cos(pi.*x(:,1)./2).*cos(pi.*x(:,2)./2)./1000;
    
        % Column vector of dofs where displacements are applied
        BC.fix_disp_dof1 = Mesh.left_dof_u;
        BC.fix_disp_dof2 = Mesh.right_dof_u;
        BC.fix_disp_dof3 = Mesh.bottom_dof_u;
        BC.fix_disp_dof4 = Mesh.top_dof_u;
        BC.fix_disp_dof = unique([BC.fix_disp_dof1;BC.fix_disp_dof2;BC.fix_disp_dof3;BC.fix_disp_dof4]);
        
        % Prescribed dispalcement
        BC.fix_disp_value = zeros(length(BC.fix_disp_dof),1);
        BC.fix_disp_value(1:2:end) = BC.UUx([Mesh.x(BC.fix_disp_dof(2:2:end)/2,1),Mesh.x(BC.fix_disp_dof(2:2:end)/2,2)]);
        BC.fix_disp_value(2:2:end) = BC.VVx([Mesh.x(BC.fix_disp_dof(2:2:end)/2,1),Mesh.x(BC.fix_disp_dof(2:2:end)/2,2)]);
        
        BC.UU_temp = zeros(length(BC.fix_disp_dof),1); 
        BC.UU_temp(1:2:end) = ones(length( BC.UU_temp(1:2:end)),1);
        BC.VV_temp = zeros(length(BC.fix_disp_dof),1); 
        BC.VV_temp(2:2:end) = ones(length( BC.VV_temp(1:2:end)),1);
        
        BC.fix_disp_value = @(t) BC.fix_disp_value.*BC.UU_temp + BC.fix_disp_value.*BC.VV_temp ;

    % Dirichlet boundary conditions (essential) - thermal
    % -----------------------------------------------------------------
    
    % Manufacture solution
    % T:= (x1,x2) → sin(pi*x1)*sin(pi*x2)
        
        BC.Tx = @(x) sin(pi.*x(:,1)).*sin(pi.*x(:,2))./1000;
        
        % Column vector of prescribed dofs where temperature is applied
        BC.fix_temp_dof1 = Mesh.left_dof_t;
        BC.fix_temp_dof2 = Mesh.right_dof_t;
        BC.fix_temp_dof3 = Mesh.bottom_dof_t;
        BC.fix_temp_dof4 = Mesh.top_dof_t;
        BC.fix_temp_dof = unique([BC.fix_temp_dof1;BC.fix_temp_dof2;BC.fix_temp_dof3;BC.fix_temp_dof4]);
        
        % Prescribed temperature
        BC.fix_temp_value = BC.Tx([Mesh.x(BC.fix_temp_dof,1),Mesh.x(BC.fix_temp_dof,2)]);

        % prescribed temperature for each dof [t1; t2; ...] 
        BC.fix_temp_value = @(t) BC.fix_temp_value;  
        
    %% Neumann BC
    % Neumann boundary conditions (natural) - mechanical
    % -----------------------------------------------------------------
        % column vector of prescribed traction dofs
        BC.traction_force_dof = [];

        % magnitude of prescribed tractions [N]
        BC.traction_force_dof_value = [];

        % NOTE: this is slower than prescribing tractions at dofs
        % column vector of prescribed traction nodes 
        BC.traction_force_node = Mesh.right_nodes;  

        % prescribed traction [t1x t1y;t2x t2y;...] [N]
        Fnode = 1e7/(length(BC.traction_force_node) - 1);
        BC.traction_force_value = Fnode*[ones(size(BC.traction_force_node)), zeros(size(BC.traction_force_node))];
        
        % find the nodes in the top right and bottom right corners
        toprightnode = find(Mesh.x(BC.traction_force_node,2) == max(Mesh.x(:,2)));
        botrightnode = find(Mesh.x(BC.traction_force_node,2) == min(Mesh.x(:,2)));
        
        BC.traction_force_value(toprightnode,1) = BC.traction_force_value(toprightnode,1)/2;
        BC.traction_force_value(botrightnode,1) = BC.traction_force_value(botrightnode,1)/2;
        
        % Make the vector into an anonymous function in time
        BC.traction_force_value = @(t) BC.traction_force_value*sin(t); 

    % Neumann boundary conditions (natural) - temperature
    % -----------------------------------------------------------------
        % column vector of prescribed traction dofs
        BC.flux_force_dof = [];

        % magnitude of prescribed tractions [N]
        BC.flux_force_dof_value = [];
        
        % NOTE: point loads at any of the element nodes can also be 
        % added as a traction.x(1
    
    % Distributed sources
    % -----------------------------------------------------------------
        E0 = Material.prop.E0;
        nu = Material.prop.nu;
        k1 = Material.Prop(1).k1;
        k2 = Material.Prop(1).k2;
        alpha = Material.Prop(1).alpha;
        C = Material.Prop(1).C;
        lambda = E0*nu/((1+nu)*(1-2*nu));
        mu = E0/(2*(1+nu));
        beta = -(3*lambda +2*mu)*alpha;
        
        % magnitude of distributed body force [N/m] [bx;by]
        BC.b = @(x,t)[-2*(E0*pi*sin(pi*x(2)/2)*(nu - 3/2)*sin(pi*x(1)/2) + ...
                       8*(-1/2 + nu)*(1 + nu)*sin(pi*x(2))*cos(pi*x(1))*beta)*pi/...
                       (16*nu^2 + 8*nu - 8) ; 
                       %
                       2*(E0*pi*cos(pi*x(2)/2)*(nu - 3/2)*cos(pi*x(1)/2) - ...
                       8*(1 + nu)*sin(pi*x(1))*cos(pi*x(2))*(-1/2 + nu)*beta)*pi/...
                       (16*nu^2 + 8*nu - 8) ]./1000;  
        
        % magnitude of distributed flux source 
        BC.s = @(x,t) sin(pi*x(1))*sin(pi*x(2))*pi^2*(k1 + k2);  

%% Initial Conditions
        BC.IC = @(t) zeros(Mesh.nsd*Mesh.nn,1);
        
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
        Control.stress_calc = 'L2projection';

        % penalty parameter for solution of static problem with 
        % LinearSolver3
        Control.beta = 10^10;
        
        % parallel inversion
        % Use parallel processing to invert the matrix. Specify number of
        % cores to use. Usually more efficient at >2e5 dofs
        Control.parallel = 1;

        % method used for solving linear problem:
        % 'LinearSolver1': Partitioning
        % 'LinearSolver2': Zeroing DOFs in stiffness matrix 
        %                   corresponding to essential boundaries
        % 'LinearSolver3': Penalty method
        Control.LinearSolver = 'LinearSolver1';
 
        % time controls
        Control.StartTime = 0;
        Control.EndTime   = 2*pi;
        NumberOfSteps     = 50;
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