function Mesh = LoadMesh(meshfile, nsd, config_dir)
%LOADMESH Load mesh data from a mesh file generated by GMSH
%   x = LOADMESH(meshfile, nsd, config_dir) is a matrix containing the 
%   spatial coordinates of the nodes (size nn x nsd in which nn is the 
%   number of nodes in the mesh and nsd is the number of spatial dimensions).
%   The mesh is defined in the file exported from GMSH (meshfile) or Hyperworks, 
%   located in the directory config_dir.  
%
%   Note: Hyperworks files only for 2D meshes
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   meshfile:       Name of the mesh file exported from GMSH (.msh) or
%                   Hypermesh (.fem)
%   nsd:            Number of spatial dimensions
%   config_dir:     Directory where mesh file is stored
% 
%   [x, conn] = LOADMESH(meshfile, nsd, config_dir) also 
%   returns the element connectivity matrix (size ne x nne in which  
%   ne is the number of elements in the mesh and nne is the number of 
%   nodes per element).
%
%   --------------------------------------------------------------------
%   Output
%   --------------------------------------------------------------------
%
%   Mesh.ext             Extension of the input file
%   
%   Outpuf for .msh files
%       Mesh.x          Spatial coordinates of mesh
%       Mesh.conn       Nodal conectivity of elements
%
%   Output for .fem files
%       Mesh.x          Spatial coordinates of mesh
%       Mesh.conn       Nodal conectivity of elements (ne x nne)
%       Mesh.BC_E       Direction of dirichlet boundary conditions (essential)
%       Mesh.BC_nE      Nodes of dirichlet boundary conditions (essential)
%       Mesh.BC_N_n     Direction of neumann boundary conditions (natural) - points
%       Mesh.BC_nN_n    Nodes of dirichlet boundary conditions (essential)
%       Mesh.BC_N_t     Nodes of neumann boundary conditions (natural) - tractions
%       Mesh.MatList    Material type of each element (1 x ne)

% Acknowledgements: Matin Parchei Esfahani

% Define if it is a .msh or .fem file
temp = strfind(meshfile,'.');
Mesh.ext = meshfile(temp:end); % Extension of the file
                        
filename = fullfile(config_dir, meshfile);
fileID = fopen(filename,'r');

s = textscan(fileID, '%s', 'delimiter', '\n');

switch Mesh.ext
    case '.msh' %GMSH file
        % start of nodes section
        n_str = find(strcmp(s{1}, '$Nodes'), 1, 'first');       
        % start of elements section
        e_str = find(strcmp(s{1}, '$Elements'), 1, 'first');    

        % number of nodes
        nnode = str2double( s{1}(n_str+1) );                   
        % number of elements
        nelem = str2double( s{1}(e_str+1) );                    

        x = zeros(nnode, nsd);

        for i = 1:nnode
            temp = s{1}(n_str+1+i);
            temp = sscanf(temp{1}(1,:), '%f');
            % nodal coordinates
            x(i,:) = temp(2:nsd+1)';                        
        end

        conn = [];

        for i = 1:nelem
            temp = s{1}(e_str+1+i);
            temp = sscanf(temp{1}(1,:), '%d');
            % equivalent type number in GMSH
            elmtyp = temp(2);                                  

            switch elmtyp
                case 5          % B8
                    nne = 8;
                    if nsd == 3
                        edg = 0;
                    else
                        edg = 1;
                    end
                case 1          % L2
                    nne = 2;    % number of nodes per element
                    edg = 1;
                case 8          % L3
                    nne = 3;    % number of nodes per element
                    edg = 1;
                case 3          % Q4
                    nne = 4;    % number of nodes per element
                    if nsd == 2
                        edg = 0;
                    else
                        edg = 1;
                    end
                case 10         % Q9
                    nne = 9;    % number of nodes per element
                    if nsd == 2
                        edg = 0;
                    else
                        edg = 1;
                    end
                case 15         % single node element
                    nne = 1;    % number of nodes per element
                    edg = 1;    
                case 12         % B27
                    nne = 27;   % number of nodes per element
                    if nsd == 3
                        edg = 0;
                    else
                        edg = 1;
                    end
                case 2          % T3 element
                    nne = 3;
                    if nsd == 2
                        edg = 0;
                    else
                        edg = 1;
                    end
                case 9          % T6 element
                    nne = 6;
                    if nsd == 2
                        edg = 0;
                    else
                        edg = 1;
                    end
            end

            if ~edg
                conn = [conn; temp(end-nne+1:end)'];
            end

        end

        clear s
        fclose(fileID);
    case '.fem' %Hypermesh file
        if nsd ~= 2
            error('Hypermesh file reader only supports 2D mesh')
        end
        
        % start of nodes section
        n_str = find(strcmp(s{1}, '$$  GRID Data'), 1, 'first') + 2;
        % end of node section
        temp_m = find(strcmp(s{1}, '$$')); % Positions where $$ is found
        n_str_e = temp_m(find(find(strcmp(s{1}, '$$'))>n_str,1, 'first')) - 1;
        % number of nodes 
        nnode = n_str_e - n_str + 1;
        
        % start of element section 
        %       Note: supports only one element type per mesh
        e_str = n_str_e + 5;
        % end of node section
        e_str_e = temp_m(find(find(strcmp(s{1}, '$$'))>e_str,1, 'first')) - 1;
        % number of elements
        nelem = e_str_e - e_str + 1;
        
        % start of CROD section for tractions Mesh.BC_N_t
        %       Note: supports only one element type per mesh
        t_str = e_str_e + 5;
        % end of section
        t_str_e = temp_m(find(find(strcmp(s{1}, '$$'))>t_str,1, 'first')) - 1;
        % number of surfaces
        nt = t_str_e - t_str + 1;
        
        % start of SPC (single - point constrains) for BC_E
        ebc_str = find(strcmp(s{1}, '$$  SPC Data'), 1, 'first') + 2;
        % end of section
        ebc_str_e = temp_m(find(find(strcmp(s{1}, '$$'))>ebc_str,1, 'first')) - 1;
        % number of supports
        nebc = ebc_str_e - ebc_str + 1;        
        
        % start of FORCE for Mesh.BC_N_n
        nbc_str = ebc_str_e + 4;
        % end of section
        nbc_str_e = temp_m(find(find(strcmp(s{1}, '$$'))>nbc_str,1, 'first')) - 1;
        % number of forces
        nnbc = nbc_str_e - nbc_str;      
        
        % Initialization of relevant matrices/vectors
        x = zeros(nnode, nsd);
        BC_N_t = zeros(nt,2);
        BC_E = zeros(nebc,2);
        BC_nE = zeros(nebc,1);
        BC_N_n = zeros(nnbc,2);
        BC_nN_n = zeros(nnbc,1);
        
        % Get grid coordinates
        for i = 1:nnode
            temp = s{1}(n_str+i-1);
            temp = char(split(temp,','));
            x(i,1) = sscanf(temp(end-3,:), '%f'); %coordinate X
            x(i,2) = sscanf(temp(end-2,:), '%f'); %coordinate Y                      
        end
        
        % Get nodal conectivity
        temp = char(split(s{1}(e_str),','));
        elmtyp =  temp(1,:);
        switch elmtyp
            case 'CQUAD4' % Q4
                nne = 4; % number of nodes per element
            case 'CTRIA3' % T3
                nne = 3; % number of nodes per element
            otherwise
                error('Mesh reader only supports Q4 and T3 for hypermesh')
        end
          
        % Initialization of conectivity matrix
        conn = zeros(nelem, nne);
        for i = 1:nelem
            temp = s{1}(e_str+i-1);
            temp = char(split(temp,','));
            for j = 1 : nne
               conn(i,j) =  sscanf(temp(end-nne-1+j,:), '%f');
            end
        end
        
        % Get natural boundary conditions for Mesh.BC_N_t
        for i = 1:nt
            temp = s{1}(t_str+i-1);
            temp = char(split(temp,','));
            BC_N_t(i,1) =  sscanf(temp(end-2,:), '%f');
            BC_N_t(i,2) =  sscanf(temp(end-1,:), '%f');
        end
        
        % Get essential boundary conditions for BC_E
        for i = 1:nebc
            temp = s{1}(ebc_str+i-1);
            temp = char(split(temp,','));
            supp = sscanf(temp(end-2,:), '%f');
            if supp == 1        % constrain on X
                BC_E(i,1) = 1;  
            elseif supp == 2    % constrain on Y
                BC_E(i,2) = 1;  
            elseif supp == 12   % constrain on X and Y
                BC_E(i,1) = 1;  
                BC_E(i,2) = 1;  
            else
                error('Essential boundary conditions should be applied only on X and Y directions in .fem file')
            end
            BC_nE(i,1) = sscanf(temp(end-3,:), '%f'); % Node label of constrain
        end       
        
        % Get natural boundary conditions for Mesh.BC_N_n
        for i = 1:nnbc
            temp = s{1}(nbc_str+i-1);
            temp = char(split(temp,','));
            BC_nN_n(i,:) =  sscanf(temp(3,:), '%f'); % Node label of force
            if abs(sscanf(temp(end-3,:), '%f')) && abs(sscanf(temp(end-1,:), '%f'))
                error('Forces are applied in a single node both in X and Y \n Define the forces separately')
            elseif abs(sscanf(temp(end-3,:), '%f')) % force on X
               BC_N_n(i,1) = 1; 
            elseif abs(sscanf(temp(end-2,:), '%f')) % force on Y
               BC_N_n(i,2) = 1;
            elseif abs(sscanf(temp(end-1,:), '%f'))
                error('Nodal forces should not be applied in the Z direction in .fem file')
            end
        end
end


% Store data in Mesh.structure
Mesh.x = x;
Mesh.conn = conn;

switch Mesh.ext
    case '.fem'
        Mesh.BC_N_t = BC_N_t;
        Mesh.BC_E = BC_E;
        Mesh.BC_nE = BC_nE;
        Mesh.BC_nN_n = BC_nN_n;
        Mesh.BC_N_n = BC_N_n;
end