function Mesh = Normal_traction(Mesh)
%Normal_traction Get the normals of the boundary elements to apply
%tractions and rearranges the elements of BC_N_e_t so that they match the
%location of BC_N_t
%
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   Mesh                
%       x          Spatial coordinates
%       conn       Mesh conectivity
%       BC_N_t     Nodes of neumann boundary conditions (natural) - tractions
%       BC_N_e_t   List of elements where tractions are applied    
%
%   --------------------------------------------------------------------
%   Output
%   --------------------------------------------------------------------
%
%   Mesh
%       n_m         list of normals to edge elements
%       BC_N_e_t_arranged       List of elements where tractions are applied    
%                               arranged in the same order as BC_N_t

% Acknowledgements: Nils Betancourt

n_m = zeros(length(Mesh.BC_N_t),2);
BC_N_e_t_arranged = zeros(length(Mesh.BC_N_t),1);

conn_t = Mesh.conn(Mesh.BC_N_e_t,:);

% Loop over all edge elements
for i = 1:length(Mesh.BC_N_t)
    
    % Get corresponding element
        [h1,~] = find(conn_t == Mesh.BC_N_t(i,1));
        [h2,~] = find(conn_t == Mesh.BC_N_t(i,2));
        temp = intersect(h1,h2);
        BC_N_e_t_arranged(i) = Mesh.BC_N_e_t(temp(1)); %element
        
    % Get normal
        element_i = BC_N_e_t_arranged(i);
        nodes_edge = Mesh.BC_N_t(i,:);
        nodes_interior = setdiff(Mesh.conn(element_i,:), nodes_edge);
        x_edge = Mesh.x(nodes_edge,:); 
        x_interior = Mesh.x(nodes_interior,:);
                
        % Equation of line in edge
        x0 = x_edge(1,1);
        y0 = x_edge(1,2);
        xe = x_edge(2,1);
        ye = x_edge(2,2);
        m = (ye-y0)/(xe-x0);
        
        n_i = abs([m, -1])./(m^2 + 1^2)^0.5;
        
        if isnan(n_i(1))
            n_i(1) = 1;
        end
        
        if (x_interior(:,2) <= m*(x_interior(:,1) - x0) + y0) %&& ...
%                 (x_interior(2,2) <= m*(x_interior(2,1) - x0) + y0)
            if max(x_interior(:,1)) < max(x_edge(:,1))
                n_m(i,:) = n_i;
            else
                n_m(i,:) = [-n_i(1), n_i(2)];
            end
        else
            if max(x_interior(:,1)) < max(x_edge(:,1))
                n_m(i,:) = [n_i(1), -n_i(2)];
            else
                n_m(i,:) = -n_i;
            end
        end   
end
    
Mesh.n_m = n_m;
Mesh.BC_N_e_t_arranged = BC_N_e_t_arranged;
    