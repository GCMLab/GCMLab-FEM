function nodeconn = NodalConn(Mesh)
%NODALCONN defines nodal connectivity in a mesh
%   nodeconn = NODALCONN(Mesh) is the nodal connectivity matrix
%   (size Mesh.nn x 8), providing the indices of elements connected to 
%   each node
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   Mesh: structure array with the fields,
%       .ne:     number of elements
%       .nn:     number of nodes
%       .conn:   element connectivity matrix (size ne x nne in which nne
%                is the number of nodes per element)

% Acknowledgements: Chris Ladubec, Matin Parchei-Esfahani

    % initialize nodal connectivity (list of elements connected to each node)
    nodeconn = zeros(Mesh.nn, 8); 

    % temp nodal counter
    temp = zeros(Mesh.nn, 1); 

    % loop through all elements
    for e = 1:Mesh.ne
        % element nodes
        enodes = Mesh.conn(e,:);

        % loop through all nodes in the element
        for n = 1:Mesh.nne
            % local node number
            nID = enodes(n);
            % increment the temp counter for the node 
            temp(nID) = temp(nID) + 1;
            % the element is added to the node connectivity
            nodeconn(nID, temp(nID)) = e;
        end
    end
end