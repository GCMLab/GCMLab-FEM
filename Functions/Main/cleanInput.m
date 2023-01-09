function [Mesh, Material, BC, Control] = cleanInput(Mesh, Material, BC, Control)
%cleanInput checks the user inputs and returns errors for invalid inputs
%   [Mesh, Material, BC, Control] = CLEANINPUT(Mesh, Material, BC, Control)
%   returns the clean input structure arrays BC 
%   
%% Output information
    err_message = sprintf('-------------------------------------------------\n');
    err_message = sprintf('%s\tfrom\t: cleanInput.m \n', err_message);
    war_message = sprintf('-------------------------------------------------\n');
    war_message = sprintf('%s\tfrom\t: cleanInput.m \n', war_message);
    err_count = 0; %Counts number of errors → counts missing data
    war_count = 0; %Counts number of warnings → shows when a default parameter has been set
%% Mesh
    err_mesh = sprintf('\t\tMesh \n');
    war_mesh = sprintf('\t\tMesh \n');

%% Material
    err_mat = sprintf('\t\tMat \n');
    war_mat = sprintf('\t\tMat \n');
    
    if isfield(Material.Prop, 'E')
        if any(Material.Prop(:).E < 0, 'all')
            err_count = err_count+1;
            err_mat = sprintf('%s\t\t\tError #%d\t:\t The elastic modulus is non-positive inside the domain\n',err_mat,err_count);
        end
    end

    if isfield(Material.Prop, 'nu')
        if any(Material.Prop(:).nu >= 0.5, 'all')
            err_count = err_count+1;
            err_mat = sprintf('%s\t\t\tError #%d\t:\t Poisson''s ratio must fall within the range -1 < nu < 0.5 within the domain - nu > 0.5 was found\n',err_mat,err_count);
        end
        if any(Material.Prop(:).nu <= -1, 'all')
            err_count = err_count+1;
            err_mat = sprintf('%s\t\t\tError #%d\t:\t Poisson''s ratio must fall within the range -1 < nu < 0.5 within the domain - nu < -1 was found\n',err_mat,err_count);
        end
    end
        
    if isfield(Material, 't')
        t_nodes = Material.t(Mesh.x);
        if any(t_nodes < 0, 'all')
            err_count = err_count+1;
            err_mat = sprintf('%s\t\t\tError #%d\t:\t The thickness is non-positive inside the domain\n',err_mat,err_count);
        end
    end


%% Boundary conditions
    err_BC = sprintf('\t\tBoundary conditions \n');
    war_BC = sprintf('\t\tBoundary conditions \n');
    
    % Prescribed dofs
    if ~any(BC.fix_disp_dof(1:Mesh.nsd:end) == Mesh.xdofs, 'all')
        war_count = war_count+1;
        war_BC = sprintf('%s\t\t\tWarning #%d\t:\t There may be insufficient prescribed degrees of freedom to prevent rigid body motion in the x-direction\n',war_BC,war_count);
    end

    if Mesh.nsd >= 2 && ~any(BC.fix_disp_dof(2:Mesh.nsd:end) == Mesh.ydofs, 'all') 
        war_count = war_count+1;
        war_BC = sprintf('%s\t\t\tWarning #%d\t:\t There may be insufficient prescribed degrees of freedom to prevent rigid body motion in the y-direction\n',war_BC,war_count);
    end

    if Mesh.nsd >= 3 && ~any(BC.fix_disp_dof(3:Mesh.nsd:end) == Mesh.zdofs, 'all') 
        war_count = war_count+1;
        war_BC = sprintf('%s\t\t\tWarning #%d\t:\t There may be insufficient prescribed degrees of freedom to prevent rigid body motion in the z-direction\n',war_BC,war_count);
    end

    if any(BC.fix_disp_dof < 0)
        err_count = err_count+1;
        err_BC = sprintf('%s\t\t\tError #%d\t:\t Prescribed degrees of freedom are negative\n',err_BC,err_count);
    end

    if any(BC.fix_disp_dof > Mesh.nDOF)
        err_count = err_count+1;
        err_BC = sprintf('%s\t\t\tError #%d\t:\t Prescribed degrees of freedom exceed the total number of degrees of freedom\n',err_BC,err_count);
    end

    % Tractions
    if any(BC.traction_force_dof < 0)
        err_count = err_count+1;
        err_BC = sprintf('%s\t\t\tError #%d\t:\t Tractions are prescribed at negative degrees of freedom\n',err_BC,err_count);
    end

    if any(BC.traction_force_node < 0)
        err_count = err_count+1;
        err_BC = sprintf('%s\t\t\tError #%d\t:\t Tractions are prescribed at negative nodes\n',err_BC,err_count);
    end

    if any(BC.traction_force_dof > Mesh.nDOF)
        err_count = err_count+1;
        err_BC = sprintf('%s\t\t\tError #%d\t:\t Tractions are prescribed at degrees of freedom which exceed the total number of degrees of freedom\n',err_BC,err_count);
    end

    if any(BC.traction_force_node > Mesh.nn)
        err_count = err_count+1;
        err_BC = sprintf('%s\t\t\tError #%d\t:\t Tractions are prescribed at nodes which exceed the total number of nodes\n',err_BC,err_count);
    end
%% Control parameters
    err_control = sprintf('\t\tControl Parameters \n');
    war_control = sprintf('\t\tControl Parameters \n');

    
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