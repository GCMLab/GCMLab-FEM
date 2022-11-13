function [Mesh, Material, BC, Control] = setDefaults(Mesh, Material, BC, Control)
%SETDEFAULTS sets the default value for variables not defined in the
%configuration file
%   [Mesh, Material, BC, Control] = SETDEFAULTS(Mesh, Material, BC, Control) 
%   returns the completed structure arrays BC and Control with the 
%   following fields,
%   

%% Mesh

    if ~isfield(Mesh, 'type')
        error('Element type is not defined - Define Mesh.type');
    end

    if ~isfield(Mesh, 'x')
        error('Nodal spatial locations are not defined - Define Mesh.x');
    end

    if ~isfield(Mesh, 'conn')
        error('Element connectivity is not defined - Define Mesh.conn');
    end

    if ~isfield(Mesh, 'DOF')
        error('DOF indices not defined - Define Mesh.DOF');
    end

    if ~isfield(Mesh, 'nodeconn')
        error('Nodal connectivity is not defined - Define Mesh.nodeconn');
    end

    if ~isfield(Mesh, 'eneighbours')
        error('Element neighbours not defined - Define Mesh.eneighbours');
    end

    if ~isfield(Mesh, 'nsd')
        Mesh.nsd = size(Mesh.x, 2);
    end

    if ~isfield(Mesh, 'ne')
        Mesh.ne = size(Mesh.conn, 1);
    end

    if ~isfield(Mesh, 'nne')
        Mesh.nne = size(Mesh.conn, 2);
    end

    if ~isfield(Mesh, 'nn')
        Mesh.nn = size(Mesh.x, 1);
    end

    if ~isfield(Mesh, 'nDOFe')
        Mesh.nDOFe = Mesh.nne*Mesh.nsd;
    end

    if ~isfield(Mesh, 'nDOF')
        Mesh.nDOF = Mesh.nn*Mesh.nsd;
    end

    if ~isfield(Mesh, 'left_nodes') ...
        || ~isfield(Mesh, 'left_dof') ...
        || ~isfield(Mesh, 'right_nodes') ...
        || ~isfield(Mesh, 'right_dof') ...
        || ~isfield(Mesh, 'xdofs') ...
        || ~isfield(Mesh, 'ydofs') ...
        || ~isfield(Mesh, 'zdofs')
        Mesh = NodeSets(Mesh);
    end

%% Material
 
    if ~isfield(Material, 'E')
        error('Elastic modulus is not defined - Define Material.E');
    end

    if ~isfield(Material, 'Dtype')
        error(['Two-dimensional approximation is not defined -', ...
                ' Define Material.Dtype']);
    end

    if Mesh.nsd==1 &&  ~isfield(Material, 'A')
        Material.A = @(x) 1;
    end

    if Mesh.nsd==2 &&  ~isfield(Material, 't')
        Material.t = @(x) 1;
    end

    if ~isfield(Material, 'nu')
        error('Poisson''s ratio is not defined - Define Material.nu');
    end
   
%% Boundary conditions

    if ~isfield(BC, 'fix_disp_dof')
	    BC.fix_disp_dof = [];
        BC.fix_disp_value = [];  
    end

    if ~isfield(BC, 'traction_force_dof')
        BC.traction_force_dof = [];
        BC.traction_force_dof_value = [];
    end

    if ~isfield(BC, 'traction_force_node')
        BC.traction_force_node = [];  
        BC.traction_force_value = [];
    end

    if ~isfield(BC, 'b')        
        BC.b = @(x)[];  
    end  

%% Control parameters

    if ~isfield(Control, 'qo')
        Control.qo = 2;
    end

    if ~isfield(Control, 'stress_calc')
        Control.stress_calc = 'center';
    end

    if ~isfield(Control, 'LinearSolver')
        Control.LinearSolver = 'LinearSolver1';
    end

    if strcmp(Control.LinearSolver, 'LinearSolver3') ...
    	&& ~isfield(Control, 'beta')
        Control.beta = 10^10;
    end
    
end