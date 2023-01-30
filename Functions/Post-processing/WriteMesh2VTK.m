function WriteMesh2VTK(filename, description, x, conn, ...
                            nodal_data, cell_data)
%WRITEMESH2VTK Exports mesh data to a file in VTK format
%   WRITEMESH2VTK(filename, description, x, conn, nodal_data, cell_data)
%   produces a vtk file (filename), with scalar nodal data (nodal_data), 
%   and scalar element data (cell_data). The nodal spatial locations are 
%   described in x, and the element connectivity in conn. 
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   filename:       Full path for the filename of the vtk file produced
%   description:    Text description of the problem 
%   x:              spatial locations of the nodes in the mesh, 
%                   specified as a matrix of size n x nsd in which n
%                   is the number of nodes and nsd is the number of 
%                   spatial dimensions.
%   conn:           Element connectivity of the mesh, specified as a 
%                   matrix of size ne x nne in which ne is the number 
%                   of elements in the mesh and nne is the number of 
%                   nodes per element. Each row of the connectivity 
%                   matrix corresponds to an element in the mesh, and
%                   the value in each column correspond to the global 
%                   node numbers associated with the element.
%   nodal_data:    Scalar nodal data, specified as a structure array 
%                   of size n x 1 in which n is the number of data arrays
%                   to output. The structure array has the following
%                   fields,
%                   .data:  column vector (size nn x 1 in which nn is 
%                           the number of nodes) containing the scalar 
%                           data 
%                   .name:  name of the data field
%                   .type:  numerical data type ('float' or 'int'); used 
%                           to determine number of significant digits to 
%                           export. Default is 'float'.
%   cell_data:      Scalar element data, specified as a structure array 
%                   of size n x 1 in which n is the number of data arrays
%                   to output. The structure array has the following
%                   fields,
%                   .data:  column vector (size ne x 1 in which ne is 
%                           the number of elements) containing the scalar 
%                           data 
%                   .name:  name of the data field
%                   .type:  numerical data type ('float' or 'int'); used 
%                           to determine number of significant digits to 
%                           export. Default is 'float'.

% Acknowledgements: Robert Gracie

min_sig_digs = 5;
max_sig_digs = 20;

if nargin < 6
    cell_data = [];
end
if nargin < 5
    nodal_data = [];
end

nn = size(x,1);             % number of nodes (points)
nsd = size(x,2);            % number of space dimensions
ne = size(conn,1);          % number of elements (cells)
nne = size(conn,2);         % number of nodes per element 
nd = length(nodal_data);   % number of scalar datum defined at each node (point);
nc = length(cell_data);     % number of scalar datum defined at each cell (element)

% fix length of nodal vector
if nsd == 1 && size(x,2)==1
    x = [x, zeros(nn,1), zeros(nn,1)];
elseif nsd == 2 && size(x,2)==2
    x = [x, zeros(nn,1)];
elseif size(x,2)<3
    error('WriteMesh2VTK: The size of the vector x is incorrect.')
end

% open the file in write mode
fid = fopen(filename, 'w');

%% header
    fprintf(fid,'%s\n','# vtk DataFile Version 2.0');
    fprintf(fid,'%s\n',['HFX mesh description: ',description]);
    fprintf(fid,'%s\n','ASCII');

%% Geometry/topology
if nne == 1
    fprintf(fid,'%s\n','DATASET POLYDATA');

    % print out the number of nodes 
    %(POINTS = nodes, float bc data is float type)
    text = ['POINTS ',num2str(nn),' float'];
    fprintf(fid,'%s\n',text);

    % print out the nodal locations
    fprintf(fid,'%f %f %f\n',x');

else
    if nsd==1 && nne == 2 % L2
        outputformat = '%d %d %d  \n';
        conn = [nne*ones(ne,1),conn-ones(ne,nne)];
        cell_type = 3;
    elseif nsd==1 && nne == 3 && nsd == 1 % L3
        outputformat = '%d %d %d %d \n';
        conn = [nne*ones(ne,1),conn-ones(ne,nne)];
        cell_type = 21;
    elseif nsd==1 && nne == 4 && nsd == 1 % L4
        outputformat = '%d %d %d %d %d \n';
        conn = [nne*ones(ne,1),conn-ones(ne,nne)];
        cell_type = 4;
    elseif nsd==2 && nne == 3 % T3
        outputformat = '%d %d %d %d \n';
        conn = [nne*ones(ne,1),conn-ones(ne,nne)];
        cell_type = 5;
    elseif nsd==2 && nne == 6 % T6
        outputformat = '%d %d %d %d %d %d %d \n';
        conn = [nne*ones(ne,1),conn-ones(ne,nne)];
        cell_type = 22;
    elseif nsd == 2 && nne == 4 % Q4 
        outputformat = '%d %d %d %d %d \n';
        conn = [4*ones(ne,1),conn(:,1:4)-ones(ne,4)];
        cell_type = 9;
    elseif nsd == 2 && nne == 9 % Q9
        % Paraview does not render Q9 elements, only Q8 - remove column 9,
        % element middle nodes
        outputformat = '%d %d %d %d %d %d %d %d %d \n';  
        middlenodes = conn(:,9);
        conn(:,9) = [];
        
        % remove middle nodes from node list
        nn = nn - length(middlenodes);
        nne = 8;
        x(middlenodes,:) = [];
        
        % remove middle nodes from nodal_data
        ndata = length(nodal_data);
        for field = 1:ndata
            fielddata = nodal_data(field).data;
            fielddata(middlenodes,:) = [];
            nodal_data(field).data = fielddata;
        end
        
        % remove middle nodes from connectivity
        middlenodes = sort(middlenodes,'ascend');
        for i = 1:length(middlenodes)
            conn(conn>middlenodes(i)) = conn(conn>middlenodes(i))-1;
            middlenodes(middlenodes>middlenodes(i)) = middlenodes(middlenodes>middlenodes(i))-1;
        end
     
        conn = [nne*ones(ne,1),conn-ones(ne,nne) ];
        cell_type = 23;   
    elseif nsd == 2 && nne == 8 %Q8
        conn = [nne*ones(ne,1),conn-ones(ne,nne)];
        outputformat = '%d %d %d %d %d %d %d %d %d \n';
        cell_type = 23;
    elseif nsd == 3 && nne == 8 %B8
        conn = [nne*ones(ne,1),conn-ones(ne,nne)];
        outputformat = '%d %d %d %d %d %d %d %d %d \n';
        cell_type = 12;
    else
        error('WriteMesh2VTK: Element not supported')
    end
    
    % Print out to VTK
     fprintf(fid,'%s\n','DATASET UNSTRUCTURED_GRID');   
    % print out the number of nodes
    %(POINTS = nodes, float bc data is float type)
    text = ['POINTS ',num2str(nn),' float'];
    fprintf(fid,'%s\n',text);

    % print out the nodal locations
    fprintf(fid, '%f %f %f\n', x');

    % print out 'CELLS', # of elements, # of nodes (+1 for each element)
    text = ['CELLS ',num2str(ne), ' ',num2str((nne+1)*ne)];
    fprintf(fid,'%s\n',text);
    
    % output the connectivity of the mesh
    % in vtk: cells are elements
    fprintf(fid,outputformat,conn');

    % print 'CELL_TYPES #of elements'
    text = ['CELL_TYPES ',num2str(ne)];
    fprintf(fid,'%s\n',text);

    % print the number of nodes per cell for each element
    fprintf(fid, '%d\n', cell_type*ones(1,ne));
end

%% Dataset Attributes

    %% Point data
    if nd > 0 % if point data is defined
        % print POINT_DATA and # of nodes for each data set
        text = ['POINT_DATA ',num2str(nn)];
        fprintf(fid,'%s\n',text);
    end

    for i = 1:nd % for each data set
        if isfield(nodal_data,'type')
            type = nodal_data(i).type;
        else
            type = 'float';
        end

        if strcmp(type,'float')
            sig_digs = real(floor(log10(nodal_data(i).data)));
            lowest_exp = min(sig_digs(~isinf(sig_digs)));
            numsigdigs = min(max_sig_digs,max([abs(lowest_exp),min_sig_digs]));
        elseif strcmp(type,'int')
            numsigdigs = '0';
        end
         
        % number of components
        numcomp = size(nodal_data(i).data,2);
        if numcomp == 3
            text = ['VECTORS ', nodal_data(i).name,' ' type ' '];
            fprintf(fid,'%s\n',text);
        else
            % print 'SCALARS data_name data_type'
            text = ['SCALARS ', nodal_data(i).name,' ' type ' ', num2str(numcomp)];
            fprintf(fid,'%s\n',text);
            text = 'LOOKUP_TABLE default';
            fprintf(fid,'%s\n',text);
        end

        % print data
        fprintf(fid,[ repmat(['%.' num2str(numsigdigs) 'f '],1,numcomp), '\n'], nodal_data(i).data');

    end

    %% Cell data
    if nc > 0
        text = ['CELL_DATA ',num2str(ne)];
        fprintf(fid,'%s\n',text);
    end

    for i = 1:nc
        if isfield(cell_data,'type')
            type = cell_data(i).type;
        else
            type = 'float';
        end

        if strcmp(type,'float')
            % find smallest significant digit
            sig_digs = real(floor(log10(cell_data(i).data)));
            lowest_exp = min(sig_digs(~isinf(sig_digs)));
            numsigdigs = min(max_sig_digs,max([abs(lowest_exp),min_sig_digs]));
        elseif strcmp(type,'int')
            numsigdigs = '0';
        end
        
        % number of components
        numcomp = size(cell_data(i).data,2);
        if numcomp == 3
            text = ['VECTORS ', cell_data(i).name,' ' type ' '];
            fprintf(fid,'%s\n',text);
        else
            % print 'SCALARS data_name data_type'
            text = ['SCALARS ', cell_data(i).name,' ' type ' ', num2str(numcomp)];
            fprintf(fid,'%s\n',text);
            text = 'LOOKUP_TABLE default';
            fprintf(fid,'%s\n',text);
        end

        % print data
        fprintf(fid,[ repmat(['%.' num2str(numsigdigs) 'f '],1,numcomp), '\n'], cell_data(i).data');
    end

fclose(fid);

end