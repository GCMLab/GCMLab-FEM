function nodeconn = NodalConn(Mesh)
    % initialize nodal connectivity (list of elements connected to each node)
    nodeconn = zeros(Mesh.nn,8); 

    % temp nodal counter
    temp = zeros(Mesh.nn,1); 

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
            nodeconn(nID,temp(nID)) = e;
        end
    end
end