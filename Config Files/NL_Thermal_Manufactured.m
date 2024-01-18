function [Mesh, Material, BC, Control] = NL_Thermal_Manufactured(config_dir, progress_on)

%% Material Properties (Solid)
        % LE1 - Linear elasticity
        % LET1 - Linear elastic with mass based damping
        % LED1 - Dynamic linear elasticity
        % ST1 - Stiffening model with 1st invariant of strain
        % ST2 - Softening model with 1st invariant of strain
        % TR2 - Stiffening model with mass based damping with 1st invariant of strain
        % VE1 - Viscoelaticity with stiffness based damping
        % TH1 - Thermal Conduction (Steady-State)
        % TH2 - Thermal Conduction (Transient)
        % THLE1 - Thermoelasticity (Steady-State)
        % THLE2 - Thermoelasticity (Transient)
        % NLTH1 - Nonlinear thermal transient
    Material.Model = 'NLTH1';
    
    % number of material properties
    Material.nmp = 1;

    % Diffusivty Properties
    % Properties Material 1
    kx0 = 10; % Diffusivity at reference temperature in the x-direction [W/m-K]
    ky0 = 1; % Diffusivity at reference temperature in the y-direction [W/m-K]
    ax = 10; % Diffusivity at slope in x-direction [W/m-K]
    ay = 1; % Diffusivity at slope in y-direction [W/m-K]
    n = 2; % Power Diffusivity Constant 

    % Nonlinear Transient Diffusion Inputs
    % 1D Heat Transfer- Define: k1, a1, n
    % 2D Heat Transfer Isotropic- Define: k1, a1, n
    % 2D Heat Transfer Orthotropic- Define: k1, k2 a1, a2 n

    Material.Prop(1).k1 = kx0; % Diffusivity at reference temperature in the x-direction [W/m-K]
    Material.Prop(1).k2 = ky0; % Diffusivity at reference temperature in the y-direction [W/m-K]
    Material.Prop(1).C = 1;   % Heat Capacity (specific heat * density) = [J/K kg] * [kg/m^3] = [J/K m^3]
    Material.Prop(1).a1 = ax; % Diffusivity at slope in x-direction [W/m-K]
    Material.Prop(1).a2 = ay; % Diffusivity at slope in y-direction [W/m-K]
    Material.Prop(1).n = n; % Power Diffusivity Constant 

    % Constitutive law: 'ISO' or 'ORTHO'
    Material.Dtype = 'ORTHO'; 
    
    % Thickness (set as default to 1)
    % 1D: [m2], 2D: [m]
    Material.t = @(x) 1;
    
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
            nex = [1;1]*10;
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
            meshFileName = 'Mesh Files\NLTH1_ABAQUS_Mesh.xlsx';
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
        BC.fix_disp_dof = [Mesh.bottom_dof; Mesh.top_dof; Mesh.left_dof; Mesh.right_dof];

        % prescribed displacement for each dof [u1; u2; ...] [m]
        BC.fix_disp_value = @(t) sqrt(Mesh.x(BC.fix_disp_dof,1).^2+Mesh.x(BC.fix_disp_dof,2).^2).*exp(-t);
        
        %% Neumann BC
    % -----------------------------------------------------------------
        % column vector of prescribed traction dofs
        BC.traction_force_dof = [];

        % magnitude of prescribed tractions [N]
        BC.traction_force_dof_value = [];

        % NOTE: this is slower than prescribing tractions at dofs
        % column vector of prescribed traction nodes 
        BC.traction_force_node = []; 

        % prescribed traction [t1x t1y;t2x t2y;...] [N]
        Fnode = 0;
        BC.traction_force_value = Fnode*[zeros(size(BC.traction_force_node)), zeros(size(BC.traction_force_node))];
       
        % Make the vector into an anonymous function in time
        BC.traction_force_value = @(t) BC.traction_force_value; 
    
        % NOTE: point loads at any of the element nodes can also be 
        % added as a traction.

        % magnitude of distributed body force [N/m] [bx;by]
            % 1D: [N/m], 2D: [N/m2]
        	% NOTE: if no body force, use '@(x)[]'
         	% NOTE: anonymous functions is defined with respect to the 

            %      variable x,  which is a vector [x(1) x(2)] = [x y]
        BC.b = @(x,t) -(exp(-t)*(ky0*x(1)^2 - x(2)^4 - x(1)^4 + kx0*x(2)^2 - 2*x(1)^2*x(2)^2 + ay*x(1)^2*(exp(-t)*(x(1)^2 + x(2)^2)^(1/2))^n + ax*x(2)^2*(exp(-t)*(x(1)^2 + x(2)^2)^(1/2))^n + ax*n*x(1)^2*(exp(-t)*(x(1)^2 + x(2)^2)^(1/2))^n + ay*n*x(2)^2*(exp(-t)*(x(1)^2 + x(2)^2)^(1/2))^n))/(x(1)^2 + x(2)^2)^(3/2);
 
%% Initial Conditions
        BC.IC = @(t) sqrt(Mesh.x(:,1).^2+Mesh.x(:,2).^2).*exp(-t);


%% Computation controls

        % quadrature order
        Control.qo = 10;

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
        Control.EndTime   = 1;
        NumberOfSteps     = 1;
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