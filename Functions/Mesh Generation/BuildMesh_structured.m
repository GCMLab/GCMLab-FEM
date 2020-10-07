function Mesh = BuildMesh_structured(Mesh)
%BUILDMESH Build a mesh
%   BuildMesh()
%
%   ----------------------------------------------------------------------
%   Created by Endrina Rivas
%       endrina.rivas@uwaterloo.ca
%       Department of Civil Engineering
%       University of Waterloo
%   Last Updated: January 2016
%   Acknowledgments: Chris Ladubec, Matin Parchei Esfahani
%   ----------------------------------------------------------------------

%% Number of nodes per element
    Mesh.nnex = NodesPerElement(Mesh.type); 
            % [nnex;nney] nnex: number of nodes per element in the x-direction
            % nney: number of nodes per element in the y-direction
    Mesh.nne = prod(Mesh.nnex);

%% Nodal spatial locations

    if ~isfield(Mesh,'x')
        % compute nodal locations
        disp([num2str(toc),': Computing nodal locations...']);
        dL = Mesh.L./(Mesh.nex.*(Mesh.nnex-1));  % node spacing
        x_temp = cell(1,Mesh.nsd);
        for i = 1:Mesh.nsd
            % nodal spatial locations
            x_temp{i} = Mesh.x1(i):dL(i):(Mesh.x1(i)+Mesh.L(i)); 
            % total number of nodes per direction
            Mesh.nnx(i) = length(x_temp{i});   
        end

        Mesh.nn = prod(Mesh.nnx);               % total number of nodes

        Mesh.x = zeros(Mesh.nn,Mesh.nsd);
        switch Mesh.nsd
            case 1
                Mesh.x = x_temp{1}';   % nodal spatial locations
            case 2
                node_count = 1;
                for nx = 1:Mesh.nnx(1)
                    for ny = 1:Mesh.nnx(2)
                        Mesh.x(node_count,1) = x_temp{1}(nx);
                        Mesh.x(node_count,2) = x_temp{2}(ny);
                        node_count = node_count + 1;
                    end
                end
            case 3
                node_count = 1;
                for nx = 1:Mesh.nnx(1)
                    for ny = 1:Mesh.nnx(2)
                        for nz = 1:Mesh.nnx(3)
                            Mesh.x(node_count,1) = x_temp{1}(nx);
                            Mesh.x(node_count,2) = x_temp{2}(ny);
                            Mesh.x(node_count,3) = x_temp{3}(nz);
                            node_count = node_count + 1;
                        end
                    end
                end
        end
    end

    % total number of nodes
    Mesh.nn = size(Mesh.x,1);

%% Nodal matrix

    switch Mesh.nsd
        case 1
            % [nodex1 nodex2 ... nodexnn]
            Mesh.nodes = 1:Mesh.nnx;
        case 2
            % Matrix of nodes in the x- and y-direction
            % [nodex1y1     nodex2y1    ... nodex_nnxy1;
            % [nodex1y2     nodex2y2    ... nodex_nnxy2;
            % [...;
            % [nodex1y_nny  nodex2y_nny ... nodex_nnx y_nny];
            Mesh.nodes = zeros(Mesh.nnx(2),Mesh.nnx(1));
            % for nx = 1:Mesh.nnx(1)
            %     Mesh.nodes(:,nx) = (Mesh.nnx(2)*(nx-1)+1):(Mesh.nnx(2)*nx);
            % end
            
            node_count = 1;
            for nx = 1:Mesh.nnx(1)
                for ny = 1:Mesh.nnx(2)
                    Mesh.nodes(ny,nx) = node_count;
                    node_count = node_count + 1;
                end
            end
            
        case 3
            % 3D Array of nodes in the x-, y- and z-directions
            % [nodex1y1z1     nodex2y1z1    ... nodex_nnxy1z1;
            % [nodex1y2z1     nodex2y2z1    ... nodex_nnxy2z1;
            % [...;
            % [nodex1y_nnyz1  nodex2y_nnyz1 ... nodex_nnx y_nnyz1];
            % ... extended in the third dimension by z2, z3...
            Mesh.nodes = zeros(Mesh.nnx(2),Mesh.nnx(1),Mesh.nnx(3));
            node_count = 1;
            for nx = 1:Mesh.nnx(1)
                for ny = 1:Mesh.nnx(2)
                    for nz = 1:Mesh.nnx(3)
                        Mesh.nodes(ny,nx,nz) = node_count;
                        node_count = node_count + 1;
                    end
                end
            end
    end   

%% Nodal DOFs
    Mesh.nDOFe = Mesh.nne*Mesh.nsd;         % number of DOF per element
    Mesh.nDOF = Mesh.nn*Mesh.nsd;           % total number of DOF
    Mesh.DOF = zeros(Mesh.nn,Mesh.nsd); 

    for sd = 1:Mesh.nsd
       Mesh.DOF(:,sd) = (sd : Mesh.nsd : (Mesh.nDOF-(Mesh.nsd-sd)))';
    end

    Mesh.ne = prod(Mesh.nex);               % total number of elements

%% Element connectivity and neighbours

    disp([num2str(toc),': Defining element connectivity...']);

    % initialize member connectivity
    Mesh.conn = zeros(Mesh.ne,Mesh.nne);

    % number of neighbors
    if Mesh.nsd == 1
        nneighbors = 2;
    else
        nneighbors = 4;
    end

    % initialize matrix of neighboring elements
    % [e1 neighbor 1, e1 neighbor 2, e1 neighbor 3, e1 neighbor 4 ...;]
    % [e2 neighbor 1, e2 neighbor 2, e2 neighbor 3, e2 neighbor 4 ...;]
    % [...]
    % [e_ne neighbor 1, e_ne neighbor 2, e_ne neighbor 3, e_ne neighbor 4 ...;]

    Mesh.eneighbours = zeros(Mesh.ne,nneighbors);

    if Mesh.nsd == 1
        for e = 1:Mesh.ne   % loop through all elements
            % if it is the first element
            if e == 1
                Mesh.conn(1,:) = 1:Mesh.nne;
            else
                Mesh.conn(e,:) = Mesh.conn(e-1,end): ...
                                (Mesh.conn(e-1,end) + Mesh.nnex - 1); 
            end
            
            Mesh.eneighbours(e,1) = e - 1;
            Mesh.eneighbours(e,2) = e + 1;
            
        end
    elseif strcmp(Mesh.type,'Q4')
        e = 1; % element counter
        % run through the elements in the x direction
        for nex = 1:Mesh.nex(1)
            % run through the elements in the y direction
            for ney = 1:Mesh.nex(2)
                % starting node number for the element
                L1 = [ney,nex];            
                L2 = L1 + [0,1];
                L3 = L1 + [1,1];
                L4 = L1 + [1,0];

                Mesh.conn(e,1) = Mesh.nodes(L1(1),L1(2));
                Mesh.conn(e,2) = Mesh.nodes(L2(1),L2(2));
                Mesh.conn(e,3) = Mesh.nodes(L3(1),L3(2));
                Mesh.conn(e,4) = Mesh.nodes(L4(1),L4(2));
                
                Mesh.eneighbours(e,1) = e - Mesh.nex(2);
                Mesh.eneighbours(e,2) = e + 1;
                Mesh.eneighbours(e,3) = e + Mesh.nex(2);
                Mesh.eneighbours(e,4) = e - 1;
            
                badleftneighbours = Mesh.eneighbours(e,1) < ney;
                badrightneighbours = Mesh.eneighbours(e,3)> ...
                                            Mesh.ne - (Mesh.nex(2) - ney);
                Mesh.eneighbours(e,1) = Mesh.eneighbours(e,1)*(1-badleftneighbours);
                Mesh.eneighbours(e,3) = Mesh.eneighbours(e,3)*(1-badrightneighbours);
                e = e + 1;

            end
            % Set non-existent neighbours to zero (for edge elements)
            % elements in this column
            colel = (nex-1)*Mesh.nex(2)+1:nex*Mesh.nex(2);
            badbotneighbours = find(Mesh.eneighbours(colel,4)<(nex-1)*Mesh.nex(2)+1);
            badtopneighbours = find(Mesh.eneighbours(colel,2)>nex*Mesh.nex(2));
            Mesh.eneighbours(colel(badbotneighbours),4) = zeros(size(badbotneighbours));
            Mesh.eneighbours(colel(badtopneighbours),2) = zeros(size(badtopneighbours));
        end
    elseif strcmp(Mesh.type,'Q9')
        e = 1; % element counter
        % run through the elements in the y-direction
        for ney = 1:Mesh.nex(2)
            % run through the elements in the x-direction
            for nex = 1:Mesh.nex(1)
                % starting node
                L1 = 2*(ney-1)*Mesh.nnx(1) + 2*(nex-1) + 1;
                
                Mesh.conn(e,1) = L1;
                Mesh.conn(e,2) = L1 + 2;
                Mesh.conn(e,3) = L1 + 2 + 2*Mesh.nnx(1);
                Mesh.conn(e,4) = L1 + 2*Mesh.nnx(1);
                Mesh.conn(e,5) = L1 + 1;
                Mesh.conn(e,6) = L1 + 2 + Mesh.nnx(1);
                Mesh.conn(e,7) = L1 + 1 + 2*Mesh.nnx(1);
                Mesh.conn(e,8) = L1 + Mesh.nnx(1);
                Mesh.conn(e,9) = L1 + 1 + Mesh.nnx(1);
                
                Mesh.eneighbours(e,1) = e - Mesh.nex(1);
                Mesh.eneighbours(e,2) = e + 1;
                Mesh.eneighbours(e,3) = e + Mesh.nex(1);
                Mesh.eneighbours(e,4) = e - 1;
                
                e = e + 1;
                
                % Set non-existent neighbours to zero (for edge elements)
                Mesh.eneighbours(Mesh.eneighbours <= 0) = ...
                zeros(size(find(Mesh.eneighbours <= 0)));

                Mesh.eneighbours(Mesh.eneighbours > Mesh.ne) = ...
                zeros(size(find(Mesh.eneighbours > Mesh.ne)));

            end
            % Set non-existent neighbours to zero (for edge elements)
            Mesh.eneighbours(Mesh.eneighbours <= 0) = ...
            zeros(size(find(Mesh.eneighbours <= 0)));

            Mesh.eneighbours(Mesh.eneighbours > Mesh.ne) = ...
            zeros(size(find(Mesh.eneighbours > Mesh.ne)));

        end             
    elseif strcmp(Mesh.type,'B8')
        e = 1;
        for nex = 1:Mesh.nex(1)
            for ney = 1:Mesh.nex(2)
                for nez = 1:Mesh.nex(3)
                    L1 = [ney,nex,nez];
                    L2 = L1 + [0,1,0];
                    L3 = L1 + [1,1,0];
                    L4 = L1 + [1,0,0];
                    L5 = L1 + [0,0,1];
                    L6 = L1 + [0,1,1];
                    L7 = L1 + [1,1,1];
                    L8 = L1 + [1,0,1];
                    
                    Mesh.conn(e,1) = Mesh.nodes(L1(1),L1(2),L1(3));
                    Mesh.conn(e,2) = Mesh.nodes(L2(1),L2(2),L2(3));
                    Mesh.conn(e,3) = Mesh.nodes(L3(1),L3(2),L3(3));
                    Mesh.conn(e,4) = Mesh.nodes(L4(1),L4(2),L4(3));
                    Mesh.conn(e,5) = Mesh.nodes(L5(1),L5(2),L5(3));
                    Mesh.conn(e,6) = Mesh.nodes(L6(1),L6(2),L6(3));
                    Mesh.conn(e,7) = Mesh.nodes(L7(1),L7(2),L7(3));
                    Mesh.conn(e,8) = Mesh.nodes(L8(1),L8(2),L8(3));
                    
                    e = e + 1;
                end
            end
        end     
    end

%% Nodal Connectivity

    % initialize nodal connectivity (list of elements connected to each node)
    Mesh.nodeconn = zeros(Mesh.nn,8); 

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
            Mesh.nodeconn(nID,temp(nID)) = e;
        end
    end

%% Node sets
       
    switch Mesh.nsd
        case 1
            Mesh.left_nodes = find(Mesh.x(:,1)==Mesh.x1(1));
            Mesh.right_nodes = find(Mesh.x(:,1)==Mesh.L(1));

            Mesh.left_dof = Mesh.DOF(Mesh.left_nodes);
            Mesh.right_dof = Mesh.DOF(Mesh.right_nodes);

            Mesh.xdofs = 1:Mesh.nDOF;
            Mesh.ydofs = [];
            Mesh.zdofs = [];

        case 2
            % Nodes and dofs on the left face of the domain
            % (left and right are defined along the x-direction)
            Mesh.left_nodes = find(Mesh.x(:,1)==Mesh.x1(1));
            Mesh.left_dof = Mesh.DOF(Mesh.left_nodes,:);        
            Mesh.left_dofx = Mesh.left_dof(:,1);
            Mesh.left_dofy = Mesh.left_dof(:,2);
            Mesh.left_dof = reshape(Mesh.left_dof',numel(Mesh.left_dof),1);

            % Nodes and dofs on the right face of the domain
            Mesh.right_nodes = find(Mesh.x(:,1)==Mesh.L(1));
            Mesh.right_dof = Mesh.DOF(Mesh.right_nodes,:);
            Mesh.right_dofx = Mesh.right_dof(:,1);
            Mesh.right_dofy = Mesh.right_dof(:,2);
            Mesh.right_dof = reshape(Mesh.right_dof',numel(Mesh.right_dof),1);
            
            % Nodes and dofs on the top face of the domain
            % (top and bottom are defined along the y-direction)
            Mesh.top_nodes = find(Mesh.x(:,2)==Mesh.L(2));
            Mesh.top_dof = Mesh.DOF(Mesh.top_nodes,:);
            Mesh.top_dofx = Mesh.top_dof(:,1);
            Mesh.top_dofy = Mesh.top_dof(:,2);
            Mesh.top_dof = reshape(Mesh.top_dof',numel(Mesh.top_dof),1);

            % Nodes and dofs on the bottom face of the domain
            Mesh.bottom_nodes = find(Mesh.x(:,2)==Mesh.x1(2));
            Mesh.bottom_dof =  Mesh.DOF(Mesh.bottom_nodes,:);
            Mesh.bottom_dofx = Mesh.bottom_dof(:,1);
            Mesh.bottom_dofy = Mesh.bottom_dof(:,2);
            Mesh.bottom_dof = reshape(Mesh.bottom_dof',numel(Mesh.bottom_dof),1);
            
            % DOFs in the x- and y- directions
            Mesh.xdofs = 1:2:Mesh.nDOF;
            Mesh.ydofs = 2:2:Mesh.nDOF;
            Mesh.zdofs = [];
            
        case 3
            
            % Nodes and dofs on the left face of the domain
            % (left and right are defined along the x-direction)
            Mesh.left_nodes = find(Mesh.x(:,1)==Mesh.x1(1));
            Mesh.left_dof = Mesh.DOF(Mesh.left_nodes,:,:);        
            Mesh.left_dofx = Mesh.left_dof(:,1);
            Mesh.left_dofy = Mesh.left_dof(:,2);
            Mesh.left_dofz = Mesh.left_dof(:,3);
            Mesh.left_dof = reshape(Mesh.left_dof',numel(Mesh.left_dof),1);

            % Nodes and dofs on the right face of the domain
            Mesh.right_nodes = find(Mesh.x(:,1)==Mesh.L(1));
            Mesh.right_dof = Mesh.DOF(Mesh.right_nodes,:);
            Mesh.right_dofx = Mesh.right_dof(:,1);
            Mesh.right_dofy = Mesh.right_dof(:,2);
            Mesh.right_dofz = Mesh.right_dof(:,3);
            Mesh.right_dof = reshape(Mesh.right_dof',numel(Mesh.right_dof),1);
            
            % Nodes and dofs on the top face of the domain
            % (top and bottom are defined along the z-direction)
            Mesh.top_nodes = find(Mesh.x(:,3)==Mesh.L(3));
            Mesh.top_dof = Mesh.DOF(Mesh.top_nodes,:);
            Mesh.top_dofx = Mesh.top_dof(:,1);
            Mesh.top_dofy = Mesh.top_dof(:,2);
            Mesh.top_dofz = Mesh.top_dof(:,3);
            Mesh.top_dof = reshape(Mesh.top_dof',numel(Mesh.top_dof),1);
            
            % Nodes and dofs on the bottom face of the domain
            Mesh.bottom_nodes = find(Mesh.x(:,3)==Mesh.x1(3));
            Mesh.bottom_dof =  Mesh.DOF(Mesh.bottom_nodes,:);
            Mesh.bottom_dofx = Mesh.bottom_dof(:,1);
            Mesh.bottom_dofy = Mesh.bottom_dof(:,2);
            Mesh.bottom_dofz = Mesh.bottom_dof(:,3);
            Mesh.bottom_dof = reshape(Mesh.bottom_dof',numel(Mesh.bottom_dof),1);
            
            % Nodes and dofs on the near face of the domain
            % (near and far are defined along the y-direction)
            Mesh.near_nodes = find(Mesh.x(:,2)==Mesh.x1(2));
            Mesh.near_dof =  Mesh.DOF(Mesh.near_nodes,:);
            Mesh.near_dofx = Mesh.near_dof(:,1);
            Mesh.near_dofy = Mesh.near_dof(:,2);
            Mesh.near_dofz = Mesh.near_dof(:,3);
            Mesh.near_dof = reshape(Mesh.near_dof',numel(Mesh.near_dof),1);
            
            % Nodes and dofs on the far face of the domain
            Mesh.far_nodes = find(Mesh.x(:,2)==Mesh.L(2));
            Mesh.far_dof =  Mesh.DOF(Mesh.far_nodes,:);
            Mesh.far_dofx = Mesh.far_dof(:,1);
            Mesh.far_dofy = Mesh.far_dof(:,2);
            Mesh.far_dofz = Mesh.far_dof(:,3);
            Mesh.far_dof = reshape(Mesh.far_dof',numel(Mesh.far_dof),1);

            Mesh.xdofs = 1:3:Mesh.nDOF;
            Mesh.ydofs = 2:3:Mesh.nDOF;
            Mesh.zdofs = 3:3:Mesh.nDOF;
    end

disp([num2str(toc),': Done generating mesh...']);
end