function Mesh = LoadMesh(meshfile, nsd, config_dir)
%LOADMESH Load mesh data from a mesh file generated by GMSH or HYPERMESH
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
%
%       Mesh.BC_E       Direction of dirichlet boundary conditions (essential - displacements)
%       Mesh.BC_nE      Nodes of dirichlet boundary conditions (essential)
%                                   OR in case sets are used
%       Mesh.c_BC_E       Direction of dirichlet boundary conditions (essential - displacements)
%       Mesh.c_BC_nE      Nodes of dirichlet boundary conditions (essential)
%
%       Mesh.BC_N_n     Direction of neumann boundary conditions (natural) - points
%       Mesh.BC_nN_n    Nodes of dirichlet boundary conditions (essential)
%                                   OR in case sets are used
%       Mesh.c_BC_N_n     Direction of neumann boundary conditions (natural) - points
%       Mesh.c_BC_nN_n    Nodes of dirichlet boundary conditions (essential)
%
%       Mesh.c_BC_N_t     Nodes of neumann boundary conditions (natural) -
%                       tractions (related to edge elements)
%       Mesh.c_BC_N_e_t   List of elements where tractions are applied
%       Mesh.c_BC_N_t_n_m        List of normals to edge elements
%       Mesh.c_BC_N_t_t_m        List of tangents to edge elements
%
%       Mesh.MatList    Material type of each element (1 x ne)


% Acknowledgements: Matin Parchei Esfahani & Nils Betancourt

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
        
        end_data_location = find(strcmp(s{1}, 'ENDDATA'));
        
        % start of nodes section
        n_str = find(strcmp(s{1}, '$$  GRID Data'), 1, 'first') + 2;
        % end of node section
        temp_m = find(strcmp(s{1}, '$$')); % Positions where $$ is found
        temp_m0 = find(strcmp(s{1}, '$')); % Positions where $ is found
        n_str_e = temp_m(find(find(strcmp(s{1}, '$$'))>n_str,1, 'first')) - 1;
        % number of nodes 
        nnode = n_str_e - n_str + 1;
        
        % start of element section 
        %       Note: supports only one element type per mesh
        e_str = n_str_e + 5;
        % end of node section
        e_str_e = temp_m(find(find(strcmp(s{1}, '$$'))>e_str,1, 'first')) - 1;
        e_str_e_1 = temp_m0(find(find(strcmp(s{1}, '$'))>e_str,1, 'first')) - 1; %Addition due to changes in Hypermesh 2022.2
        if e_str_e > e_str_e_1
            e_str_e = e_str_e_1;
        end
        % number of elements
        nelem = e_str_e - e_str + 1;

        % start of CROD section for tractions Mesh.BC_N_t
        %       Note: supports only one element type per mesh
        if ~isempty(find(strcmp(s{1}, '$$  CROD Elements'), 1, 'first'))
            t_str = e_str_e + 6;
            % end of section
            t_str_e = temp_m(find(find(strcmp(s{1}, '$$'))>t_str,1, 'first')) - 1;
            t_str_e_1 = temp_m0(find(find(strcmp(s{1}, '$'))>t_str,1, 'first')) - 1; %Addition due to changes in Hypermesh 2022.2
            if t_str_e > t_str_e_1
                t_str_e = t_str_e_1;
            end
            % number of surfaces
            nt = t_str_e - t_str + 1;
        else
            nt = 0;
        end
        
        % start of SPC (single - point constrains) for BC_E 
        spc_text = find(strcmp(s{1}, '$$  SPC Data'), 1, 'first');
        % Obtained location for Essential boundary conditions
        if ~isempty(spc_text)
            ebc_str = find(strcmp(s{1}, '$$  SPC Data'), 1, 'first') + 2;
            % end of section
            ebc_str_e = temp_m(find(find(strcmp(s{1}, '$$'))>ebc_str,1, 'first')) - 1;
            % number of supports
            nebc = ebc_str_e - ebc_str + 1*(ebc_str_e ~=  end_data_location);        
        else
            nebc = 0;
        end
        
        % start of PLOAD4 for Mesh.BC_N_e_t
        %       Note: Gets elements where tractions are applied
        if ~isempty(find(strcmp(s{1}, '$$  PLOAD4 Data'), 1, 'first'))
            t_e_str = ebc_str_e + 4;
            % end of section
            t_e_str_e = temp_m(find(find(strcmp(s{1}, '$$'))>t_e_str,1, 'first')) - 1;
            % number of elements (must be the same as nt)
            nt_e = t_e_str_e - t_e_str + 1*(t_e_str_e ~=  end_data_location);
            if nt_e ~= nt
                error('The number of edge elements is not compatible with the number of elements where tractions are applied\nNumber of PLOAD4 and CROD are different')
            end
        end
        
        % start of FORCE for Mesh.BC_N_n
        if ~isempty(find(strcmp(s{1}, '$$  FORCE Data'), 1, 'first'))
            % CONTINUE FROM HERE - PROBLEM IF BC_N_e_t does not exists!!
            if exist('t_e_str_e', 'var')
                % Tractions where applied as distributed boundary conditions
                nbc_str = t_e_str_e + 4;
            else
                nbc_str = ebc_str_e + 4;
            end
            % end of section
            nbc_str_e = temp_m(find(find(strcmp(s{1}, '$$'))>nbc_str,1, 'first')) - 1;
            % number of forces
            nnbc = nbc_str_e - nbc_str;
        else
            nnbc = 0;
        end
        
        % Initialization of relevant matrices/vectors
        x = zeros(nnode, nsd);
        
        BC_N_t = zeros(nt,2);
        BC_N_e_t = zeros(nt,1);
        BC_N_t_set = zeros(nt,1);
        n_m = zeros(nt,1);
        t_m = zeros(nt,1);
        
%         BC_E = zeros(nebc,4); % For rotated DOF case
        BC_E = zeros(nebc,2);
        BC_nE = zeros(nebc,1);
        BC_E_set = zeros(nebc,1);
        
        BC_N_n = zeros(nnbc,2);
        BC_nN_n = zeros(nnbc,1);
        BC_N_n_set = zeros(nnbc,1);
        
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
        % Correction for node labeling;
        min_conn = min(min(conn));
        conn(:,:)  = conn(:,:) - min_conn + 1;

        % Get natural boundary conditions for Mesh.BC_N_t 
        % Matrix of nodes that define were tractions are applied
        % and matrix of elements where nodes are appleid
        for i = 1:nt
            temp = s{1}(t_str+i-1);
            temp = char(split(temp,','));
            BC_N_t(i,1) =  sscanf(temp(end-2,:), '%f');
            BC_N_t(i,2) =  sscanf(temp(end-1,:), '%f');
            
            temp = s{1}(t_e_str+i-1);
            temp = char(split(temp,','));
            BC_N_e_t(i) =  sscanf(temp(3,:), '%f');
            BC_N_t_set(i,1) =  sscanf(temp(4,:), '%f'); % Set label of distributed loads
        end
        BC_N_t(:,:) = BC_N_t(:,:) - min_conn + 1;            

        % Get essential boundary conditions for BC_E
        for i = 1:nebc
            temp = s{1}(ebc_str+i-1);
            temp = char(split(temp,','));
            supp = sscanf(temp(end-2,:), '%f');
            % Case for DOF aligned with cartesian system
            if supp == 1        % constrain on X
                BC_E(i,1) = 1;         
            elseif supp == 2    % constrain on Y
                BC_E(i,2) = 1;  
            elseif supp == 12   % constrain on X and Y
                BC_E(i,1) = 1;  
                BC_E(i,2) = 1;  
            % Case for rotated DOF
            elseif supp == 4        % constrain on X'
%                 BC_E(i,3) = 1;         
            elseif supp == 5    % constrain on Y'
%                 BC_E(i,4) = 1;  
            elseif supp == 45   % constrain on X' and Y'
%                 BC_E(i,3) = 1;  
%                 BC_E(i,4) = 1;    
            else
                error('Essential boundary conditions should be applied only on X and Y directions in .fem file')
            end
            BC_nE(i,1) = sscanf(temp(end-3,:), '%f'); % Node label of constrain
            BC_E_set(i,1) = sscanf(temp(end-1,:), '%f'); % Set label of constrain
        end   
        % Correction for node labeling;
        BC_nE(:,:) = BC_nE(:,:) - min_conn + 1;
        % Collect sets from data
        u_E = unique(BC_E_set);
        if length(u_E) ~= 1  && ~isempty(u_E) %Sets are used
            % Create cell collector
            c_BC_E = cell(1,length(u_E));
            c_BC_nE = cell(1,length(u_E));      
            for i = 1:length(u_E) % loop over all sets
                temp = find(BC_E_set == i);
                c_BC_E{i} = BC_E(temp,:);
                c_BC_nE{i} = BC_nE(temp,:);
            end
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
               BC_N_n_set(i,1) =  sscanf(temp(end-3,:), '%f'); % Set label of force
            elseif abs(sscanf(temp(end-2,:), '%f')) % force on Y
               BC_N_n(i,2) = 1;
               BC_N_n_set(i,1) =  sscanf(temp(end-2,:), '%f'); % Set label of force
            elseif abs(sscanf(temp(end-1,:), '%f'))
                error('Nodal forces should not be applied in the Z direction in .fem file')
            end
        end
        % Correction for node labeling;
        BC_nN_n(:,:) = BC_nN_n(:,:) - min_conn + 1;
        % Collect sets from data
        u_N_n = unique(BC_N_n_set);
        if length(u_N_n) ~= 1 && ~isempty(u_N_n)%Sets are used
            % Create cell collector
            c_BC_N_n = cell(1,length(u_N_n));
            c_BC_nN_n = cell(1,length(u_N_n));      
            for i = 1:length(u_N_n) % loop over all sets
                temp = find(BC_N_n_set == i);
                c_BC_N_n{i} = BC_N_n(temp,:);
                c_BC_nN_n{i} = BC_nN_n (temp,:);
            end
        end
end


% Store data in Mesh.structure
Mesh.x = x;
Mesh.conn = conn;

switch Mesh.ext
    case '.fem'
                        
        % Get normals to the tractions
        if ~isempty(BC_N_t)
             [n_m, t_m, BC_N_e_t,  BC_N_t_set] = Normal_traction(Mesh, BC_N_t, BC_N_e_t, BC_N_t_set);
        end
        
        u_N_t = unique(BC_N_t_set);
        if length(u_N_t) ~= 1 && ~isempty(u_N_t) %Sets are used on tractions applied as distributed forces
            % Collect sets from data
            if length(u_N_t) ~= 1 %Sets are used
                % Create cell collector
                c_BC_N_e_t = cell(1,length(u_N_t));   
                c_BC_n_m = cell(1,length(u_N_t));   
                c_BC_t_m = cell(1,length(u_N_t));
                c_BC_N_t = cell(1,length(u_N_t));
                for i = 1:length(u_N_t) % loop over all sets
                    temp = find(BC_N_t_set == i);
                    c_BC_N_e_t{i} = BC_N_e_t(temp,:);
                    c_BC_n_m{i} = n_m(temp,:);
                    c_BC_t_m{i} = t_m(temp,:);
                    c_BC_N_t{i} = BC_N_t(temp,:);
                end
            end
            Mesh.c_BC_N_t = c_BC_N_t;
            Mesh.c_BC_N_e_t = c_BC_N_e_t;
            Mesh.c_BC_N_t_n_m = c_BC_n_m;
            Mesh.c_BC_N_t_t_m = c_BC_t_m;
        else
            Mesh.c_BC_N_t = cell(1,1); Mesh.c_BC_N_t{1} = BC_N_t;
            Mesh.c_BC_N_e_t = cell(1,1); Mesh.c_BC_N_e_t{1} = BC_N_e_t;
            Mesh.c_BC_N_t_n_m = cell(1,1); Mesh.c_BC_N_t_n_m{1} = n_m;
            Mesh.c_BC_N_t_t_m = cell(1,1); Mesh.c_BC_N_t_t_m{1} = t_m;
        end
        
        if length(u_E) ~= 1 && ~isempty(u_E) %Sets are used on essential boundary conditions
            Mesh.c_BC_E = c_BC_E;
            Mesh.c_BC_nE = c_BC_nE;
        else 
            Mesh.BC_E = BC_E;
            Mesh.BC_nE = BC_nE;
        end
        
        if length(u_N_n) ~= 1 && ~isempty(u_N_n) %Sets are used on tractions applied at nodes directly 
            Mesh.c_BC_N_n = c_BC_N_n;
            Mesh.c_BC_nN_n = c_BC_nN_n;
        else
            Mesh.BC_nN_n = BC_nN_n;
            Mesh.BC_N_n = BC_N_n;
        end

end
