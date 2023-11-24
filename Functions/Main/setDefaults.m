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
    
    if isfield(Mesh,'ext') 
        switch Mesh.ext
            case '.fem'
                % The following variables are not prescribed when the data
                % is obtained from a .fem file
            case '.msh'
                if ~isfield(Mesh, 'left_nodes') ...
                    || ~isfield(Mesh, 'left_dof') ...
                    || ~isfield(Mesh, 'right_nodes') ...
                    || ~isfield(Mesh, 'right_dof') ...
                    || ~isfield(Mesh, 'xdofs') ...
                    || ~isfield(Mesh, 'ydofs') ...
                    || ~isfield(Mesh, 'zdofs')
                    war_count = war_count+1;
                    war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh not defined - has been set as Mesh = MeshNodeSets(Mesh)\n',war_mesh,war_count);
                    Mesh = NodeSets(Mesh);
                end
        end
    else
        % For structure mesh
        if ~isfield(Mesh, 'left_nodes') ...
                    || ~isfield(Mesh, 'left_dof') ...
                    || ~isfield(Mesh, 'right_nodes') ...
                    || ~isfield(Mesh, 'right_dof') ...
                    || ~isfield(Mesh, 'xdofs') ...
                    || ~isfield(Mesh, 'ydofs') ...
                    || ~isfield(Mesh, 'zdofs')
                    war_count = war_count+1;
                    war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh not defined - has been set as Mesh = MeshNodeSets(Mesh)\n',war_mesh,war_count);
                    Mesh = NodeSets(Mesh);
        end
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
        war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh.nsd not defined - has been set as Mesh.nsd = size(Mesh.x, 2)\n',war_mesh,war_count);
        Mesh.nsd = size(Mesh.x, 2);
    end

    if ~isfield(Mesh, 'ne')
        war_count = war_count+1;
        war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh.ne not defined - has been set as Mesh.ne = size(Mesh.conn, 1)\n',war_mesh,war_count);
        Mesh.ne = size(Mesh.conn, 1);
    end

    if ~isfield(Mesh, 'nne')
        war_count = war_count+1;
        war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh.nne not defined - has been set as Mesh.nne = size(Mesh.conn, 2)\n',war_mesh,war_count);
        Mesh.nne = size(Mesh.conn, 2);
    end

    if ~isfield(Mesh, 'nn')
        war_count = war_count+1;
        war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh.nn not defined - has been set as Mesh.nn = size(Mesh.x, 1)\n',war_mesh,war_count);
        Mesh.nn = size(Mesh.x, 1);
    end

    if ~isfield(Mesh, 'nDOFe')
        war_count = war_count+1;
        war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh.nDOFe not defined - has been set as Mesh.nDOFe = Mesh.nne*Mesh.nsd\n',war_mesh,war_count);
        Mesh.nDOFe = Mesh.nne*Mesh.nsd;
    end

    if ~isfield(Mesh, 'nDOF')
        war_count = war_count+1;
        war_mesh = sprintf('%s\t\t\tWarning #%d\t:\t Mesh.nDOF not defined - has been set as Mesh.nDOF = Mesh.nn*Mesh.nsd\n',war_mesh,war_count);
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
        % exception: coupled problem where Dtype_mech and Dtype_therm need
        % to be prescribed
        if ~isfield(Material, 'Dtype_mech') && ~isfield(Material, 'Dtype_therm')
            err_count = err_count+1;
            err_mat = sprintf('%s\t\t\tError #%d\t:\t Two-dimensional approximation is not defined - Define Material.Dtype\n',err_mat,err_count);
        end
    end

    if ~isfield(Material, 't')
        err_count = err_count+1;
        err_mat = sprintf('%s\t\t\tError #%d\t:\t Model thickness is not defined - Define Material.t\n',err_mat,err_count);
    end
    
    if ~isfield(Material, 'Model')
        Material.Model = 'LE1';
        war_count = war_count+1;
        war_mat = sprintf('%s\t\t\tError #%d\t:\t Material model not defined, set to linear elastic\n',err_mat,err_count);
    end
    
    if ~isfield(Material, 'ExternalForceFile')
        % Assume uncoupled problme unless otherwise specified.
        Material.ExternalForceFile = 'getFext';
    end    
   
%% Boundary conditions
    err_BC = sprintf('\t\tBoundary conditions \n');
    war_BC = sprintf('\t\tBoundary conditions \n');

    if ~isfield(BC, 'fix_disp_dof')
        war_count = war_count+1;
        war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.fix_disp_dof and BC.fix_disp_value not defined - has been set as []\n',war_BC,war_count);
	    BC.fix_disp_dof = [];
        BC.fix_disp_value = [];  
    end

    if ~isfield(BC, 'traction_force_dof')
        war_count = war_count+1;
        war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.traction_force_dof and BC.traction_force_dof_value not defined - has been set as []\n',war_BC,war_count);
        BC.traction_force_dof = [];
        BC.traction_force_dof_value = [];
    end

    if ~isfield(BC, 'traction_force_node')
        war_count = war_count+1;
        war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.traction_force_node and BC.traction_force_node_value not defined - has been set as []\n',war_BC,war_count);
        BC.traction_force_node = [];  
        BC.traction_force_value = [];
    end
    
    if ~isfield(BC, 'c_N_t_f')
        war_count = war_count+1;
        war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.c_N_t_f  not defined - has been set as @(x,t)[]\n',war_BC,war_count);
        BC.c_N_t_f = @(x,t)[];  
    end

    if ~isfield(BC, 'c_N_t_flag')
        war_count = war_count+1;
        war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.c_N_t_flag not defined - has been set as []\n',war_BC,war_count);
        BC.c_N_t_flag = [];  
    end

    if ~isfield(BC, 'b')    
        war_count = war_count+1;
        war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.b not defined - has been set as @(x,t)[]\n',war_BC,war_count);
        BC.b = @(x,t)[];  
    end  
    
    if ~isfield(BC, 'IC')
        % Assume start from rest unless otherwise specified. 
        %war_count = war_count+1;
        %war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.IC not defined - set as Mesh.nsd*Mesh.nn\n',war_BC,war_count);
        BC.IC = @(t) zeros(Mesh.nDOFn*Mesh.nn,1);
    end

    % set boundary conditions defaults for coupled thermoelasticity problem
    if Material.ProblemType == 3
        if ~isfield(BC, 'fix_temp_dof')
            war_count = war_count+1;
            war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.fix_temp_dof and BC.fix_temp_value not defined - has been set as []\n',war_BC,war_count);
    	    BC.fix_temp_dof = [];
            BC.fix_temp_value = [];
        end

        if ~isfield(BC, 'flux_dof')
            war_count = war_count+1;
            war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.flux_dof and BC.flux_dof_value not defined - has been set as []\n',war_BC,war_count);
            BC.flux_dof = [];
            BC.flux_dof_value = [];
        end

        if ~isfield(BC, 'flux_node')
            war_count = war_count+1;
            war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.flux_node and BC.flux_node_value not defined - has been set as []\n',war_BC,war_count);
            BC.flux_node = [];
            BC.flux_value = [];
        end

        if ~isfield(BC, 's')
            war_count = war_count+1;
            war_BC = sprintf('%s\t\t\tWarning #%d\t:\t BC.s not defined - has been set as @(x)[]\n',war_BC,war_count);
            BC.s = @(x)[];
        end
        
        % define unique set of fixed DOFs
        BC.fix_dof = [BC.fix_disp_dof; BC.fix_temp_dof];
        BC.fix_value = @(t) [BC.fix_disp_value(t); BC.fix_temp_value(t)];

    else
        BC.fix_dof = BC.fix_disp_dof;
        BC.fix_value = BC.fix_disp_value;
    end

%% Control parameters
    err_control = sprintf('\t\tControl Parameters \n');
    war_control = sprintf('\t\tControl Parameters \n');
    
    if ~isfield(Control, 'qo')
        war_count = war_count+1;
        war_control = sprintf('%s\t\t\tWarning #%d\t:\t Control.qo not defined - has been set as 2\n',war_control,war_count);
        Control.qo = 2;
    end

    if ~isfield(Control, 'stress_calc')
        war_count = war_count+1;
        war_control = sprintf('%s\t\t\tWarning #%d\t:\t Control.stress_calc not defined - has been set as ''center''\n',war_control,war_count);
        Control.stress_calc = 'center';
    end

    if ~isfield(Control, 'LinearSolver')
        war_count = war_count+1;
        war_control = sprintf('%s\t\t\tWarning #%d\t:\t Control.LinearSolver  not defined - has been set as ''LinearSolver1''\n',war_control,war_count);
        Control.LinearSolver = 'LinearSolver1';
    end
    
    if ~isfield(Control, 'parallel')
        Control.parallel = 1;
    end
    

    if strcmp(Control.LinearSolver, 'LinearSolver3') ...
    	&& ~isfield(Control, 'beta')
        war_count = war_count+1;
        war_control = sprintf('%s\t\t\tWarning #%d\t:\t Control.beta  not defined - has been set as 10^10\n',war_control,war_count);
        Control.beta = 10^10;
    end
    
    % Aitken Relaxation Parameters (optional)
    if ~isfield(Control, 'aitkenON')
        Control.aitkenON = 0; % Off by default
    end
    
    if ~isfield(Control, 'aitkenRange')
        Control.aitkenRange = [0.01,2];
    end
    
    if ~isfield(Control, 'aitkenNeg')
        Control.aitkenNeg = 0; % Off by default
    end
    
    if ~isfield(Control, 'relaxDOFs')
        Control.relaxDOFs   = setdiff(1:Mesh.nDOF,BC.fix_disp_dof);
    else
        % Ensure that the relaxDOFs do not contain any fixed BCs, as this
        % interferes with the calculation of the relaxation parameter
        Control.relaxDOFs   = setdiff(Control.relaxDOFs,BC.fix_disp_dof);
    end
    
    % If these variables are not defined, assume a static analysis
    if ~isfield(Control, 'StartTime')
        Control.StartTime = 0;
    end
    
    if ~isfield(Control, 'EndTime')
        Control.EndTime = 1;
    end
    
    if ~isfield(Control, 'TimeStep')
        Control.TimeStep = 1;
    end
    
    if ~isfield(Control, 'dSave')
        Control.dSave = 0;
    end
    
    if ~isfield(Control, 'plotLoadDispl')
        Control.plotLoadDispl = 0;
    end
    
    if ~isfield(Control, 'TimeCase')
        Control.TimeCase = 'static';
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