function [Mesh, Material, BC, Control] = cleanInput(Mesh, Material, BC, Control)
%cleanInput checks the user inputs and returns errors for invalid inputs
%   [Mesh, Material, BC, Control] = CLEANINPUT(Mesh, Material, BC, Control)
%   returns the clean input structure arrays BC 
%   
adsf= 0;
%% Mesh


%% Material

if isfield(Material.E)
    E_nodes = Material.E(Mesh.x);
    if any(E_nodes < 0, 'all')
        error('The elastic modulus is non-positive inside the domain.')
    end
end

if isfield(Material.nu)
    nu_nodes = Material.nu(Mesh.x);
    if any(nu_nodes >= 0.5, 'all')
        error('Poisson''s ratio must fall within the range -1 < nu < 0.5 within the domain.')
    end
    if any(nu_nodes <= -1, 'all')
        error('Poisson''s ratio must fall within the range -1 < nu < 0.5 within the domain.')
    end
end

if isfield(Material.t)
    t_nodes = Material.t(Mesh.x);
    if any(t_nodes < 0, 'all')
        error('The thickness is non-positive inside the domain.')
    end
end


%% Boundary conditions


%% Control parameters

    
end