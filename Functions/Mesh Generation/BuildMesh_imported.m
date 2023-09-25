function Mesh = BuildMesh_imported(meshFileName, nsd, config_dir, progress_on, Q8_reduced, problemtype)
%BUILDMESH_IMPORTED Import GMSH file or HYPERMESH file
%    Mesh = BUILDMESH_IMPORTED(meshFileName, nsd, config_dir) is a structure 
%   array with the mesh description. The mesh is imported from the GMSH 
%   or HYPERMESH file meshFileName in the directory config_dir. The mesh has spatial 
%   dimension nsd. The import is capable of handling unstructured meshes. 
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   meshFileName:   Name of the mesh file exported from GMSH
%   nsd:            Number of spatial dimensions
%   config_dir:     Directory where mesh file is stored
%   Q8_reduced:     Optional input to be defined when Q8 with reduced
%                   integration is desired.This case is oibtained from que
%                   Q9 mesh
% 
%   --------------------------------------------------------------------
%   Output
%   --------------------------------------------------------------------
%   Mesh: structure array with the fields,
%       .type:          the topological class of finite element; it is in 
%                       the general form 'topology-#of nodes' ie a three 
%                       node triangle is T3 a four node quadralateral is 
%                       Q4 a 4 node tetrahedra is H4 a 27 node brick is 
%                       B27 etc. Presently defined are L2, Q4, and Q9. 
%       .nsd:           Number of spatial dimensions
%       .ne:            Total number of elements in the mesh
%       .nne:           Vector of number of nodes per element (size nsd x 1)
%       .nn:            Total number of nodes 
%       .nDOFe:         Number of DOFs per element
%       .nDOF:          Total number of DOFs
%       .x:             Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)
%       .conn:          Array of element connectivity (size ne x nne)
%       .eneighbours:   Array of neighbouring elements (size ne x nneighbours
%                       in which nneighbours is 1 for 1D elements and 4
%                       for 2D elements)
%       .DOF:           Array of DOF indices (size nn x nsd)
%       .nodeconn:      Array of nodal connectivity (size nn x 8)
%                       contains the indices of elements connected to 
%                       each node
%       .left_nodes     Nodes on the left edge of the domain
%       .left_dof       DOFs on the left edge of the domain
%       .right_nodes    Nodes on the right edge of the domain
%       .right_dof      DOFs on the right edge of the domain
%       .xdofs          DOFs in the x-direction
%       .ydofs          DOFs in the y-direction
%       .zdofs          DOFs in the z-direction
%   Two-dimensional meshes also contain the fields,
%       .top_nodes      Nodes on the top edge of the domain
%       .top_dof        DOFs on the top edge of the domain
%       .top_dofx       DOFs on the top boundary in the x-direction
%       .top_dofy       DOFs on the top boundary in the y-direction
%       .bottom_nodes   Nodes on the bottom edge of the domain
%       .bottom_dof     DOFs on the bottom edge of the domain
%       .bottom_dofx    DOFs on the bottom boundary in the x-direction
%       .bottom_dofy    DOFs on the bottom boundary in the y-direction
%       .left_dofx      DOFs on the left boundary in the x-direction
%       .left_dofy      DOFs on the left boundary in the y-direction
%       .right_dofx     DOFs on the right boundary in the x-direction
%       .right_dofy     DOFs on the right boundary in the y-direction
%   Three-dimensional meshes also contain the fields, 
%       .near_nodes     nodes on the nearest face of the domain
%       .near_dof       DOFs on the nearest face of the domain
%       .near_dofx      DOFs on the near face in the x-direction
%       .near_dofy      DOFs on the near face in the y-direction
%       .near_dofz      DOFs on the near face in the z-direction
%       .far_nodes      Nodes on the farthest face of the domain
%       .far_dof        DOFs on the farthest face of the domain
%       .far_dofx       DOFs on the far face in the x-direction
%       .far_dofy       DOFs on the far face in the y-direction
%       .far_dofz       DOFs on the far face in the z-direction
%       .left_dofz      DOFs on the left face in the z-direction
%       .right_dofz     DOFs on the right face in the z-direction
%       .top_dofz       DOFs on the top face in the z-direction
%       .bottom_dofz    DOFs on the bottom face in the z-direction

% Acknowledgments: Matin Parchei Esfahani

%% ProblemType
if nargin < 6 
   problemtype =  1; % Assume solid mechanics equilibrium problem by default
end

%% Load GMSH file
    % nodal matrix and member connectivity
    [Mesh] = LoadMesh(meshFileName, nsd, config_dir);

%% Mesh properties
    % total number of elements
    Mesh.ne = size(Mesh.conn,1);
    % number of nodes per element
    Mesh.nne = size(Mesh.conn,2);
    % total number of nodes
    Mesh.nn = size(Mesh.x,1);
    % number of spatial dimensions
    Mesh.nsd = nsd;

    switch Mesh.nne
        case 3
            Mesh.type = 'T3';
        case 4
            Mesh.type = 'Q4';
        case 6
            Mesh.type = 'T6';
        case 9
            switch nargin
                case 5
                    Mesh.type = Q8_reduced; %Q8 with reduced integration;
                otherwise
                    Mesh.type = 'Q9';
            end
        case 8
            Mesh.type = 'B8';
    end
            
%% Nodal DOFs
    Mesh = NodeDOFs(Mesh, problemtype);

%% Element connectivity and neighbours
    if progress_on
        disp([num2str(toc),': Defining element connectivity...']);
    end

%% Element neighbours 
    if strcmp(Mesh.type,'Q4')% NOTE: Only works for Q4 elements at the moment
        
        % list of elements connected to each node
        Mesh.nodeconn = NodalConn(Mesh);
        
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
    else
        Mesh.eneighbours = 'NA - currently only available for Q4 element';
    end
    
%% Transform Q9 → Q8
if strcmp(Mesh.ext, '.msh')
    switch Mesh.type 
        case 'Q8'
            % GMSH does not generate Q8 elements - remove column 9

            middlenodes = Mesh.conn(:,9);
            Mesh.conn(:,9) = [];

            % remove middle nodes from node list
            Mesh.nn = Mesh.nn - length(middlenodes);
            Mesh.nne = 8;
            Mesh.nDOFe = 16;
            Mesh.nDOF = Mesh.nDOF - length(middlenodes)* nsd;
            Mesh.x(middlenodes,:) = [];

            % Update Mesh.DOF
            Mesh.DOF = zeros(Mesh.nn, Mesh.nsd); 
            for sd = 1:Mesh.nsd
               Mesh.DOF(:,sd) = (sd : Mesh.nsd : (Mesh.nDOF-(Mesh.nsd-sd)))';
            end

            % remove middle nodes from connectivity
            middlenodes = sort(middlenodes,'ascend');
            conn = Mesh.conn; 
            for i = 1:length(middlenodes)
                conn(conn>middlenodes(i)) = conn(conn>middlenodes(i))-1;
                middlenodes(middlenodes>middlenodes(i)) = middlenodes(middlenodes>middlenodes(i))-1;
            end

            % Update Mesh.nodeconn
            Mesh.conn = conn;
        otherwise
            % No changes to Mesh. structure

    end
end
%% Nodal Connectivity
    % list of elements connected to each node
    switch Mesh.type
        case 'Q4'
            %This was already done in line 125
        case 'Q8'
            %This has not been implemented for Q8s with reduced integration
            Mesh.nodeconn = [];
        otherwise
            Mesh.nodeconn = NodalConn(Mesh);
    end

%% Node sets
    Mesh = NodeSets(Mesh, problemtype);
   
if progress_on
    disp([num2str(toc),': Done generating mesh...']);
end
end