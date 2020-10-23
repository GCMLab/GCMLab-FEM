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


%% Load GMSH file
    % nodal matrix and member connectivity
    [Mesh.x, Mesh.conn] = LoadMesh(meshFileName, nsd, nsd, config_dir);

%% Mesh properties
    % total number of elements
    Mesh.ne = size(Mesh.conn,1);
    % number of nodes per element
    Mesh.nne = size(Mesh.conn,2);
    % total number of nodes
    Mesh.nn = size(Mesh.x,1);

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
            
%% Nodal DOFs
    Mesh = NodeDOFs(Mesh);

%% Element connectivity and neighbours
    disp([num2str(toc),': Defining element connectivity...']);

%% Nodal Connectivity
    % list of elements connected to each node
    Mesh.nodeconn = NodalConn(Mesh);

%% Element neighbours 
    % NOTE: Only works for Q4 elements at the moment
    Mesh.eneighbours = zeros(Mesh.ne,4);  % element neighbours (share an edge)
    for e = 1:Mesh.ne
        % list of elements which share at least one node with element e       
        elist = Mesh.nodeconn(Mesh.conn(e,:),:);
        % reshape into a vector      
        elist = reshape(elist,[numel(elist),1]);
        % remove element e from the vector
        elist = setdiff(unique(elist),[e,0]);  

        switch Mesh.type
            case 'Q4'
                pattern = [1,2,3,4,1];
        end

        for i = 1:Mesh.nne
            [r1,~] = find(Mesh.conn(elist,:) == Mesh.conn(e, pattern(i)));
            [r2,~] = find(Mesh.conn(elist,:) == Mesh.conn(e, pattern(i+1)));

            if isempty(intersect(r1,r2))
                Mesh.eneighbours(e,i) = 0;
            else
                Mesh.eneighbours(e,i) = elist(intersect(r1,r2));
            end
        end
    end

%% Node sets
    Mesh = NodeSets(Mesh);

disp([num2str(toc),': Done generating mesh...']);
end