function [Mesh, Material, BC, Control] = cleanInput(Mesh, Material, BC, Control)
%CLEANINPUT checks the user inputs and returns errors for invalid inputs
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
    
    if isfield(Mesh, 'c_BC_N_t') && ~isfield(Mesh,'c_BC_N_t_n_m')
        err_count = err_count+1;
        err_mesh = sprintf('%s\t\t\tError #%d\t:\t Edge elements have been defined but the normals have not been computed as input. Includ Mesh.c_BC_N_t_n_m in input.\n',err_mesh,err_count);
    end
    
    if isfield(Mesh, 'c_BC_N_t') && ~isfield(Mesh,'c_BC_N_t_t_m')
        err_count = err_count+1;
        err_mesh = sprintf('%s\t\t\tError #%d\t:\t Edge elements have been defined but the tangents have not been computed as input. Includ Mesh.c_BC_N_t_t_m in input.\n',err_mesh,err_count);
    end
    
    if ~isfield(Mesh, 'c_BC_N_t') && isfield(Mesh,'c_BC_N_t_n_m')
        err_count = err_count+1;
        err_mesh = sprintf('%s\t\t\tError #%d\t:\t Normals for traction computation have been defined, but edge elements are not included in input. Include Mesh.c_BC_N_t in input.\n',err_mesh,err_count);
    end
    
    if isfield(Mesh, 'c_BC_N_t') && isfield(Mesh,'c_BC_N_t_n_m')
        if length(Mesh.c_BC_N_t) ~= length(Mesh.c_BC_N_t_n_m)
            err_count = err_count+1;
            err_mesh = sprintf('%s\t\t\tError #%d\t:\t The number of sets of edge elements is not compatible with the number of sets given for the normal vectors. Review Mesh.BC_N_t and Mesh.c_BC_N_t_n_m.\n',err_mesh,err_count);
        end
    end
    
    if isfield(Mesh, 'c_BC_N_t') && isfield(Mesh,'c_BC_N_t_t_m')
        if length(Mesh.c_BC_N_t) ~= length(Mesh.c_BC_N_t_t_m)
            err_count = err_count+1;
            err_mesh = sprintf('%s\t\t\tError #%d\t:\t The number of sets of edge elements is not compatible with the number of sets given for the tangent vectors. Review Mesh.BC_N_t and Mesh.c_BC_N_t_t_m.\n',err_mesh,err_count);
        end
    end
    
    if isfield(Mesh, 'c_BC_N_t') && isfield(Mesh,'c_BC_N_t_n_m')
        if length(Mesh.c_BC_N_t) == length(Mesh.c_BC_N_t_n_m)
            for i = 1:length(Mesh.c_BC_N_t_n_m)
                BC_N_t_temp = Mesh.c_BC_N_t{i};
                BC_N_t_n_m_temp = Mesh.c_BC_N_t_n_m{i};
                if length(BC_N_t_temp) ~= length(BC_N_t_n_m_temp)
                    err_count = err_count+1;
                    err_mesh = sprintf('%s\t\t\tError #%d\t:\t The number of edge elements and normals in set %d is not compatible. Check size of Mesh.c_BC_N_t{%d} and Mesh.c_BC_N_t_n_m{%d}.\n',err_mesh,err_count, i, i , i);    
                end
            end
        end
        
        if length(Mesh.c_BC_N_t) == length(Mesh.c_BC_N_t_t_m)
            for i = 1:length(Mesh.c_BC_N_t_t_m)
                BC_N_t_temp = Mesh.c_BC_N_t{i};
                BC_N_t_t_m_temp = Mesh.c_BC_N_t_n_m{i};
                if length(BC_N_t_temp) ~= length(BC_N_t_t_m_temp)
                    err_count = err_count+1;
                    err_mesh = sprintf('%s\t\t\tError #%d\t:\t The number of edge elements and tangents in set %d is not compatible. Check size of Mesh.c_BC_N_t{%d} and Mesh.c_BC_N_t_t_m{%d}.\n',err_mesh,err_count, i, i , i);    
                end
            end
        end
    end
    
    if ~strcmp(Mesh.type, 'Q4')
        if ~strcmp(Mesh.type, 'T3')
            if isfield(Mesh, 'c_BC_N_t')
                err_count = err_count+1;
                err_mesh = sprintf('%s\t\t\tError #%d\t:\tThe computation of tractions using edge elements only works for Q4 and T3 elements.\n',err_mesh,err_count);
            end
        end
    end
    
%% Material
    err_mat = sprintf('\t\tMat \n');
    war_mat = sprintf('\t\tMat \n');
    
    % Elasticity problems
    
    if isfield(Material.Prop, 'E')
        for ii = 1 : Material.nmp
            if any(Material.Prop(ii).E < 0, 'all')
                err_count = err_count+1;
                err_mat = sprintf('%s\t\t\tError #%d\t:\t The elastic modulus is non-positive inside the domain\n',err_mat,err_count);
            end
        end
    end

    if isfield(Material.Prop, 'nu')
        for ii = 1 : Material.nmp
            if any(Material.Prop(ii).nu >= 0.5, 'all')
                err_count = err_count+1;
                err_mat = sprintf('%s\t\t\tError #%d\t:\t Poisson''s ratio must fall within the range -1 < nu < 0.5 within the domain - nu > 0.5 was found\n',err_mat,err_count);
            end
            if any(Material.Prop(ii).nu <= -1, 'all')
                err_count = err_count+1;
                err_mat = sprintf('%s\t\t\tError #%d\t:\t Poisson''s ratio must fall within the range -1 < nu < 0.5 within the domain - nu < -1 was found\n',err_mat,err_count);
            end
        end
    end
        
    if isfield(Material, 't')
        t_nodes = Material.t(Mesh.x);
        if any(t_nodes < 0, 'all')
            err_count = err_count+1;
            err_mat = sprintf('%s\t\t\tError #%d\t:\t The thickness is non-positive inside the domain\n',err_mat,err_count);
        end
    end
    
    if ~isfield(Material, 'ProblemType')
        Material.ProblemType = 1;
        war_count = war_count+1;
        war_mat = sprintf('%s\t\t\tWarning #%d\t:\t Outdated config file format, problem type not defined. Deformation problem is assumed. \n',war_mat,war_count);
    end
    
    % Diffusion problems
    if isfield(Material.Prop, 'k1')
        for ii = 1 : Material.nmp
            if any(Material.Prop(ii).k1 < 0, 'all')
                err_count = err_count+1;
                err_mat = sprintf('%s\t\t\tError #%d\t:\t The conductivity is non-positive inside the domain\n',err_mat,err_count);
            end
        end
    end

    if isfield(Material.Prop, 'C')
        for ii = 1 : Material.nmp
            if any(Material.Prop(ii).C < 0, 'all')
                err_count = err_count+1;
                err_mat = sprintf('%s\t\t\tError #%d\t:\t The heat capacity is non-positive inside the domain\n',err_mat,err_count);
            end
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

    if Material.ProblemType ~= 2
        if Mesh.nsd >= 2 && ~any(BC.fix_disp_dof(2:Mesh.nsd:end) == Mesh.ydofs, 'all') 
            war_count = war_count+1;
            war_BC = sprintf('%s\t\t\tWarning #%d\t:\t There may be insufficient prescribed degrees of freedom to prevent rigid body motion in the y-direction\n',war_BC,war_count);
        end

        if Mesh.nsd >= 3 && ~any(BC.fix_disp_dof(3:Mesh.nsd:end) == Mesh.zdofs, 'all') 
            war_count = war_count+1;
            war_BC = sprintf('%s\t\t\tWarning #%d\t:\t There may be insufficient prescribed degrees of freedom to prevent rigid body motion in the z-direction\n',war_BC,war_count);
        end
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
    
    if isfield(Mesh, 'c_BC_N_t') && isfield(Mesh,'c_BC_N_t_n_m')
        if ~isfield(BC, 'c_N_t_f') 
            err_count = err_count+1;
            err_BC = sprintf('%s\t\t\tError #%d\t:\t Tractions functions compatible with the use of edge elements for integration have not been applied. Include BC.c_N_t_f in input.\n',err_BC,err_count);
        end
        
        if ~isfield(BC, 'c_N_t_flag') 
            err_count = err_count+1;
            err_BC = sprintf('%s\t\t\tError #%d\t:\t Vector related to orientation of traction application has not been defined. Include BC.c_N_t_flag in input.\n',err_BC,err_count);
        end
        
        if length(BC.c_N_t_f) ~= length(Mesh.c_BC_N_t)
            err_count = err_count+1;
            err_BC = sprintf('%s\t\t\tError #%d\t:\t Number of functions BC.c_N_t_f is not compatible with the number of sets of edge elements.\n',err_BC,err_count);
        end
        
        if length(BC.c_N_t_flag) ~= length(Mesh.c_BC_N_t)
            err_count = err_count+1;
            err_BC = sprintf('%s\t\t\tError #%d\t:\t Number of elements in BC.c_N_t_flag is not compatible with the number of functions defined in Mesh.c_BC_N_t.\n',err_BC,err_count);
        end
        
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