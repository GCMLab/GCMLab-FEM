function [Mesh, Material, BC, Control] = setDefaults(Mesh, Material, BC, Control)
%SETDEFAULTS sets the default value for variables not defined in the
%configuration file
%   [Mesh, Material, BC, Control] = SETDEFAULTS(Mesh, Material, BC, Control) 
%   returns the completed structure arrays BC and Control with the 
%   following fields,
%   
%% Output information
    err_message = sprintf('-------------------------------------------------\n');
    err_message = sprintf('%s\tfrom\t: setDefaults.m \n', err_message);
    war_message = sprintf('-------------------------------------------------\n');
    war_message = sprintf('%s\tfrom\t: setDefaults.m \n', war_message);
    err_count = 0; %Counts number of errors → counts missing data
    war_count = 0; %Counts number of warnings → shows when a default parameter has been set
%% Mesh
    err_mesh = sprintf('\t\tMesh \n');
    war_mesh = sprintf('\t\tMesh \n');
    
    if ~isfield(Mesh, 'left_nodes') ...
        || ~isfield(Mesh, 'left_dof') ...
        || ~isfield(Mesh, 'right_nodes') ...
        || ~isfield(Mesh, 'right_dof') ...
        || ~isfield(Mesh, 'xdofs') ...
        || ~isfield(Mesh, 'ydofs') ...
        || ~isfield(Mesh, 'zdofs')
        war_count = war_count+1;
        war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh not defined - set as Mesh = MeshNodeSets(Mesh)\n',war_mesh,war_count);
        Mesh = NodeSets(Mesh);
    end
    
    if ~isfield(Mesh, 'type')
        err_count = err_count+1;
        err_mesh = sprintf('%s\t\t\tError #%d\t:\t Element type is not defined - Define Mesh.type\n',err_mesh,err_count);
    end

    if ~isfield(Mesh, 'x')
        err_count = err_count+1;
        err_mesh = sprintf('%s\t\t\tError #%d\t:\t Nodal spatial locations are not defined - Define Mesh.x\n',err_mesh,err_count);
    end

    if ~isfield(Mesh, 'conn')
        err_count = err_count+1;
        err_mesh = sprintf('%s\t\t\tError #%d\t:\t Element connectivity is not defined - Define Mesh.conn\n',err_mesh,err_count);
    end

    if ~isfield(Mesh, 'DOF')
        err_count = err_count+1;
        err_mesh = sprintf('%s\t\t\tError #%d\t:\t DOF indices not defined - Define Mesh.DOF\n',err_mesh,err_count);
    end

    if ~isfield(Mesh, 'nodeconn')
        err_count = err_count+1;
        err_mesh = sprintf('%s\t\t\tError #%d\t:\t Nodal connectivity is not defined - Define Mesh.nodeconn\n',err_mesh,err_count);
    end

    if ~isfield(Mesh, 'eneighbours')
        err_count = err_count+1;
        err_mesh = sprintf('%s\t\t\tError #%d\t:\t Element neighbours not defined - Define Mesh.eneighbours\n',err_mesh,err_count);
    end

    if ~isfield(Mesh, 'MatList')
        err_count = err_count+1;
        err_mesh = sprintf('%s\t\t\tError #%d\t:\t Material assigned to elements not defined - Define Mesh.MatList\n',err_mesh,err_count);
    end
    
    if ~isfield(Mesh, 'nsd')
        war_count = war_count+1;
        war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh.nsd not defined - set as Mesh.nsd = size(Mesh.x, 2)\n',war_mesh,war_count);
        Mesh.nsd = size(Mesh.x, 2);
    end

    if ~isfield(Mesh, 'ne')
        war_count = war_count+1;
        war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh.ne not defined - set as Mesh.ne = size(Mesh.conn, 1)\n',war_mesh,war_count);
        Mesh.ne = size(Mesh.conn, 1);
    end

    if ~isfield(Mesh, 'nne')
        war_count = war_count+1;
        war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh.nne not defined - set as Mesh.nne = size(Mesh.conn, 2)\n',war_mesh,war_count);
        Mesh.nne = size(Mesh.conn, 2);
    end

    if ~isfield(Mesh, 'nn')
        war_count = war_count+1;
        war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh.nn not defined - set as Mesh.nn = size(Mesh.x, 1)\n',war_mesh,war_count);
        Mesh.nn = size(Mesh.x, 1);
    end

    if ~isfield(Mesh, 'nDOFe')
        war_count = war_count+1;
        war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh.nDOFe not defined - set as Mesh.nDOFe = Mesh.nne*Mesh.nsd\n',war_mesh,war_count);
        Mesh.nDOFe = Mesh.nne*Mesh.nsd;
    end

    if ~isfield(Mesh, 'nDOF')
        war_count = war_count+1;
        war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh.nDOF not defined - set as Mesh.nDOF = Mesh.nn*Mesh.nsd\n',war_mesh,war_count);
        Mesh.nDOF = Mesh.nn*Mesh.nsd;
    end

%% Material
    err_mat = sprintf('\t\tMat \n');
    war_mat = sprintf('\t\tMat \n');

    if ~isfield(Material, 'Prop')
        err_count = err_count+1;
        err_mat = sprintf('%s\t\t\tError #%d\t:\t Material properties are not defined - Define Material.Prop\n',err_mat,err_count);
    end

    if ~isfield(Material, 'Dtype')
        err_count = err_count+1;
        err_mat = sprintf('%s\t\t\tError #%d\t:\t Two-dimensional approximation is not defined - Define Material.Dtype\n',err_mat,err_count);
    end

    if ~isfield(Material, 't')
        err_count = err_count+1;
        err_mat = sprintf('%s\t\t\tError #%d\t:\t Model thickness is not defined - Define Material.t\n',err_mat,err_count);
    end
   
%% Boundary conditions
    err_BC = sprintf('\t\tBoundary conditions \n');
    war_BC = sprintf('\t\tBoundary conditions \n');

    if ~isfield(BC, 'fix_disp_dof')
        war_count = war_count+1;
        war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.fix_disp_dof and BC.fix_disp_value not defined - set as []\n',war_BC,war_count);
	    BC.fix_disp_dof = [];
        BC.fix_disp_value = [];  
    end

    if ~isfield(BC, 'traction_force_dof')
        war_count = war_count+1;
        war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.traction_force_dof and BC.traction_force_dof_value not defined - set as []\n',war_BC,war_count);
        BC.traction_force_dof = [];
        BC.traction_force_dof_value = [];
    end

    if ~isfield(BC, 'traction_force_node')
        war_count = war_count+1;
        war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.traction_force_node and BC.traction_force_node_value not defined - set as []\n',war_BC,war_count);
        BC.traction_force_node = [];  
        BC.traction_force_value = [];
    end

    if ~isfield(BC, 'b')    
        war_count = war_count+1;
        war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.b not defined - set as @(x)[]\n',war_BC,war_count);
        BC.b = @(x)[];  
    end  

%% Control parameters
    err_control = sprintf('\t\tControl Parameters \n');
    war_control = sprintf('\t\tControl Parameters \n');
    
    if ~isfield(Control, 'qo')
        war_count = war_count+1;
        war_control = sprintf('%s\t\t\tWarning #%d\t:\t Control.qo not defined - set as 2\n',war_control,war_count);
        Control.qo = 2;
    end

    if ~isfield(Control, 'stress_calc')
        war_count = war_count+1;
        war_control = sprintf('%s\t\t\tWarning #%d\t:\t Control.stress_calc not defined - set as ''center''\n',war_control,war_count);
        Control.stress_calc = 'center';
    end

    if ~isfield(Control, 'LinearSolver')
        war_count = war_count+1;
        war_control = sprintf('%s\t\t\tWarning #%d\t:\t Control.LinearSolver  not defined - set as ''LinearSolver1''\n',war_control,war_count);
        Control.LinearSolver = 'LinearSolver1';
    end

    if strcmp(Control.LinearSolver, 'LinearSolver3') ...
    	&& ~isfield(Control, 'beta')
        war_count = war_count+1;
        war_control = sprintf('%s\t\t\tWarning #%d\t:\t Control.beta  not defined - set as 10^10\n',war_control,war_count);
        Control.beta = 10^10;
    end
    
%% Output

war_message = sprintf('%s %s %s %s %s', war_message, war_mesh, war_mat, war_BC, war_control);

err_message = sprintf('%s %s %s %s %s', err_message, err_mesh, err_mat, err_BC, err_control);

if war_count > 0
    warning('\n%s', war_message)
end

if err_count > 0 
    error(err_message)
end

end