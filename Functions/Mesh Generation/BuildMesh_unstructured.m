function Mesh = BuildMesh_unstructured(Mesh, Control)
%BUILDMESH_unstructured Build an unstructured mesh from a GMSH file
%   BuildMesh()
%
%   ----------------------------------------------------------------------
%   Created by Bruce Gee
%       b3gee@uwaterloo.ca
%       Department of Civil Engineering
%       University of Waterloo
%   Last Updated: October 2020 by Bruce Gee
%   Acknowledgments: Chris Ladubec, Matin Parchei Esfahani, Endrina Rivas
%   ----------------------------------------------------------------------


% Load GMSH file
[nodes, conn] = LoadMesh(Mesh.MeshFileName, Mesh.nsd, Mesh.nsd, Control);



% total number of nodes
Mesh.ne = size(conn,1);
% number of nodes per element
Mesh.nne = size(conn,2);
% total number of nodes
Mesh.nn = size(nodes,1);

switch Mesh.nne
    case 3
        Mesh.type = 'T3';
    case 4
        Mesh.type = 'Q4';
    case 6
        Mesh.type = 'T6';
    case 9
        Mesh.type = 'Q9';
    case 8
        Mesh.type = 'B8';
end
        

%% Nodal matrix
    Mesh.x = nodes;
    
%% Nodal DOFs
    Mesh.nDOFe = Mesh.nne*Mesh.nsd;         % number of DOF per element
    Mesh.nDOF = Mesh.nn*Mesh.nsd;           % total number of DOF
    Mesh.DOF = zeros(Mesh.nn,Mesh.nsd); 

    for sd = 1:Mesh.nsd
       Mesh.DOF(:,sd) = (sd : Mesh.nsd : (Mesh.nDOF-(Mesh.nsd-sd)))';
    end


%% Element connectivity and neighbours

    disp([num2str(toc),': Defining element connectivity...']);
    
    % initialize member connectivity
    Mesh.conn = conn;



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
            Mesh.left_nodes = find(Mesh.x(:,1)==min(Mesh.x(:,1)));
            Mesh.right_nodes = find(Mesh.x(:,1)==max(Mesh.x(:,1)));

            Mesh.left_dof = Mesh.DOF(Mesh.left_nodes);
            Mesh.right_dof = Mesh.DOF(Mesh.right_nodes);

            Mesh.xdofs = 1:Mesh.nDOF;
            Mesh.ydofs = [];
            Mesh.zdofs = [];

        case 2
            % Nodes and dofs on the left face of the domain
            % (left and right are defined along the x-direction)
            Mesh.left_nodes = find(Mesh.x(:,1)==min(Mesh.x(:,1)));
            Mesh.left_dof = Mesh.DOF(Mesh.left_nodes,:);        
            Mesh.left_dofx = Mesh.left_dof(:,1);
            Mesh.left_dofy = Mesh.left_dof(:,2);
            Mesh.left_dof = reshape(Mesh.left_dof',numel(Mesh.left_dof),1);

            % Nodes and dofs on the right face of the domain
            Mesh.right_nodes = find(Mesh.x(:,1)==max(Mesh.x(:,1)));
            Mesh.right_dof = Mesh.DOF(Mesh.right_nodes,:);
            Mesh.right_dofx = Mesh.right_dof(:,1);
            Mesh.right_dofy = Mesh.right_dof(:,2);
            Mesh.right_dof = reshape(Mesh.right_dof',numel(Mesh.right_dof),1);
            
            % Nodes and dofs on the top face of the domain
            % (top and bottom are defined along the y-direction)
            Mesh.top_nodes = find(Mesh.x(:,2)==max(Mesh.x(:,2)));
            Mesh.top_dof = Mesh.DOF(Mesh.top_nodes,:);
            Mesh.top_dofx = Mesh.top_dof(:,1);
            Mesh.top_dofy = Mesh.top_dof(:,2);
            Mesh.top_dof = reshape(Mesh.top_dof',numel(Mesh.top_dof),1);

            % Nodes and dofs on the bottom face of the domain
            Mesh.bottom_nodes = find(Mesh.x(:,2)==min(Mesh.x(:,2)));
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
            Mesh.left_nodes = find(Mesh.x(:,1)==min(Mesh.x(:,1)));
            Mesh.left_dof = Mesh.DOF(Mesh.left_nodes,:,:);        
            Mesh.left_dofx = Mesh.left_dof(:,1);
            Mesh.left_dofy = Mesh.left_dof(:,2);
            Mesh.left_dofz = Mesh.left_dof(:,3);
            Mesh.left_dof = reshape(Mesh.left_dof',numel(Mesh.left_dof),1);

            % Nodes and dofs on the right face of the domain
            Mesh.right_nodes = find(Mesh.x(:,1)==max(Mesh.x(:,1)));
            Mesh.right_dof = Mesh.DOF(Mesh.right_nodes,:);
            Mesh.right_dofx = Mesh.right_dof(:,1);
            Mesh.right_dofy = Mesh.right_dof(:,2);
            Mesh.right_dofz = Mesh.right_dof(:,3);
            Mesh.right_dof = reshape(Mesh.right_dof',numel(Mesh.right_dof),1);
            
            % Nodes and dofs on the top face of the domain
            % (top and bottom are defined along the z-direction)
            Mesh.top_nodes = find(Mesh.x(:,2)==max(Mesh.x(:,2)));
            Mesh.top_dof = Mesh.DOF(Mesh.top_nodes,:);
            Mesh.top_dofx = Mesh.top_dof(:,1);
            Mesh.top_dofy = Mesh.top_dof(:,2);
            Mesh.top_dofz = Mesh.top_dof(:,3);
            Mesh.top_dof = reshape(Mesh.top_dof',numel(Mesh.top_dof),1);
            
            % Nodes and dofs on the bottom face of the domain
            Mesh.bottom_nodes = find(Mesh.x(:,2)==min(Mesh.x(:,2)));
            Mesh.bottom_dof =  Mesh.DOF(Mesh.bottom_nodes,:);
            Mesh.bottom_dofx = Mesh.bottom_dof(:,1);
            Mesh.bottom_dofy = Mesh.bottom_dof(:,2);
            Mesh.bottom_dofz = Mesh.bottom_dof(:,3);
            Mesh.bottom_dof = reshape(Mesh.bottom_dof',numel(Mesh.bottom_dof),1);
            
            % Nodes and dofs on the near face of the domain
            % (near and far are defined along the y-direction)
            Mesh.near_nodes = find(Mesh.x(:,3)==min(Mesh.x(:,3)));
            Mesh.near_dof =  Mesh.DOF(Mesh.near_nodes,:);
            Mesh.near_dofx = Mesh.near_dof(:,1);
            Mesh.near_dofy = Mesh.near_dof(:,2);
            Mesh.near_dofz = Mesh.near_dof(:,3);
            Mesh.near_dof = reshape(Mesh.near_dof',numel(Mesh.near_dof),1);
            
            % Nodes and dofs on the far face of the domain
            Mesh.far_nodes = find(Mesh.x(:,3)==max(Mesh.x(:,3)));
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