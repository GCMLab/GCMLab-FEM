function [Mesh, Material, BC, Control] = cleanInput(Mesh, Material, BC, Control)
%cleanInput checks the user inputs and returns errors for invalid inputs
%   [Mesh, Material, BC, Control] = CLEANINPUT(Mesh, Material, BC, Control)
%   returns the clean input structure arrays BC 
%   

%% Mesh


%% Material

if isfield(Material, 'E')
    E_nodes = Material.E(Mesh.x);
    if any(E_nodes < 0, 'all')
        error('The elastic modulus is non-positive inside the domain.')
    end
end

if isfield(Material, 'nu')
    nu_nodes = Material.nu(Mesh.x);
    if any(nu_nodes >= 0.5, 'all')
        error('Poisson''s ratio must fall within the range -1 < nu < 0.5 within the domain.')
    end
    if any(nu_nodes <= -1, 'all')
        error('Poisson''s ratio must fall within the range -1 < nu < 0.5 within the domain.')
    end
end

if isfield(Material, 't')
    t_nodes = Material.t(Mesh.x);
    if any(t_nodes < 0, 'all')
        error('The thickness is non-positive inside the domain.')
    end
end


%% Boundary conditions
% Prescribed dofs
if ~any(BC.fix_disp_dof == Mesh.xdofs, 'all')
    warning('There may be insufficient prescribed degrees of freedom to prevent rigid body motion in the x-direction.')
end

if Mesh.nsd >= 2 && ~any(BC.fix_disp_dof == Mesh.ydofs, 'all') 
    warning('There may be insufficient prescribed degrees of freedom to prevent rigid body motion in the y-direction.')
end

if Mesh.nsd >= 3 && ~any(BC.fix_disp_dof == Mesh.zdofs, 'all') 
    warning('There may be insufficient prescribed degrees of freedom to prevent rigid body motion in the z-direction.')
end

if any(BC.fix_disp_dof < 0)
    error('Prescribed degrees of freedom are negative.')
end

if any(BC.fix_disp_dof > Mesh.nDOF)
   error('Prescribed degrees of freedom exceed the total number of degrees of freedom.') 
end

% Tractions
if any(BC.traction_force_dof < 0)
    error('Tractions are prescribed at negative degrees of freedom.')
end

if any(BC.traction_force_node < 0)
    error('Tractions are prescribed at negative nodes.')
end

if any(BC.traction_force_dof > Mesh.nDOF)
    error('Tractions are prescribed at degrees of freedom which exceed the total number of degrees of freedom.')
end

if any(BC.traction_force_node > Mesh.nn)
    error('Tractions are prescribed at nodes which exceed the total number of nodes.')
end
%% Control parameters

    
end