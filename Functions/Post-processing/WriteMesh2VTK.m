function WriteMesh2VTK(filename,description,nodes,conn,...
                            scalar_data,cell_data)
% This function writes FEM mesh data to a file in VTK format
%
% INPUTS:
% nodes = array of nodal coordinates (n x nsd)
% conn = connectivity matrix (ne x nne)
%
% Acknowledgements: Robert Gracie, Endrina Rivas

min_sig_digs = 5;

if nargin < 6
    cell_data = [];
end
if nargin < 5
    scalar_data = [];
end

nn = size(nodes,1);         % number of nodes (points)
nsd = size(nodes,2);        % number of space dimensions
ne = size(conn,1);          % number of elements (cells)
nne = size(conn,2);         % number of nodes per element 
nd = length(scalar_data);   % number of scalar datum defined at each node (point);
nc = length(cell_data);     % number of scalar datum defined at each cell (element)

% open the file in write mode
fid = fopen(filename, 'w');

% header
fprintf(fid,'%s\n','# vtk DataFile Version 2.0');
fprintf(fid,'%s\n',['HFX mesh description: ',description]);
fprintf(fid,'%s\n','ASCII');
fprintf(fid,'%s\n','DATASET UNSTRUCTURED_GRID');

% print out the number of nodes 
%(POINTS = nodes, float bc data is float type)
text = ['POINTS ',num2str(nn),' float'];
fprintf(fid,'%s\n',text);

% if it's a 2D model, add a row of zeros to the end of nodes that is the
% same length as the number of nodes (to signify zero in the z-direction), 
% and if it's a 1D model, add two rows of zeros to end of nodes (for y- and
% z- directions)
if nsd == 1
    nodes = [nodes,zeros(nn,1),zeros(nn,1)];
elseif nsd == 2
    nodes = [nodes,zeros(nn,1)];
end

% print out the nodal locations
fprintf(fid,'%f %f %f\n',nodes');

% print out 'CELLS', #of elements, # of nodes (+1 for each element)
text = ['CELLS ',num2str(ne), ' ',num2str((nne+1)*ne)];
fprintf(fid,'%s\n',text);

% output the connectivity of the mesh
% in vtk: cells are elements
if nsd == 2 && (nne == 4 || nne ==9) % Q4 or Q9
      outputformat = '%d %d %d %d %d \n';
      conn = [4*ones(ne,1),conn(:,1:4)-ones(ne,4)];
      cell_type = 9;
elseif nne == 2 % L2
    outputformat = '%d %d %d  \n';
    conn = [nne*ones(ne,1),conn-ones(ne,nne)];
    cell_type = 3;
elseif nne == 3 && nsd == 1 % L3
    outputformat = '%d %d %d %d \n';
    conn = [nne*ones(ne,1),conn-ones(ne,nne)];
    cell_type = 21;
elseif nne == 4 && nsd == 1 % L4
    outputformat = '%d %d %d %d %d \n';
    conn = [nne*ones(ne,1),conn-ones(ne,nne)];
    cell_type = 4;
elseif nne == 3 % T3
    outputformat = '%d %d %d %d \n';
    conn = [nne*ones(ne,1),conn-ones(ne,nne)];
    cell_type = 5;
else %B8
      conn = [nne*ones(ne,1),conn-ones(ne,nne)];
      outputformat = '%d %d %d %d %d %d %d %d %d \n';
      cell_type = 12;
end
   
fprintf(fid,outputformat,conn');

% print 'CELL_TYPES #of elements'
text = ['CELL_TYPES ',num2str(ne)];
fprintf(fid,'%s\n',text);

% print the number of nodes per cell for each element
fprintf(fid,'%d\n',cell_type*ones(1,ne));

if nd > 0 % if point data is defined
    % print POINT_DATA and #of nodes for each data set
    text = ['POINT_DATA ',num2str(nn)];
    fprintf(fid,'%s\n',text);
end

for i = 1:nd % for each data set
    if isfield(scalar_data,'type')
        type = scalar_data(i).type;
    else
        type = 'float';
    end

    if strcmp(type,'float')
        % find smallest significant digit
        mindata = num2str(abs(scalar_data(i).data),'%e');
        test = min(str2num(mindata(:,end-2:end)));
        numsigdigs = num2str(max([abs(test),min_sig_digs]));
    elseif strcmp(type,'int')
        numsigdigs = '0';
    end
    
    if str2double(numsigdigs) > 30
%         disp(['NUM SIG DIGS = ' num2str(numsigdigs)])
        numsigdigs = '30';
    end
    
    % print 'SCALARS data_name data_type'
    text = ['SCALARS ',scalar_data(i).name,' ' type ' 1'];
    fprintf(fid,'%s\n',text);
    text = 'LOOKUP_TABLE default';
    fprintf(fid,'%s\n',text);
    % print data
    fprintf(fid,['%.' numsigdigs 'f\n'],scalar_data(i).data);
end

if nc > 0
text = ['CELL_DATA ',num2str(ne)];
fprintf(fid,'%s\n',text);
end

for i = 1:nc
    if isfield(scalar_data,'type')
        type = scalar_data(i).type;
    else
        type = 'float';
    end

    if strcmp(type,'float')
        % find smallest significant digit
        mindata = num2str(abs(scalar_data(i).data),'%e');
        test = min(str2num(mindata(:,end-2:end)));
        numsigdigs = num2str(max([abs(test),min_sig_digs]));
    elseif strcmp(type,'int')
        numsigdigs = '0';
    end
    
    text = ['SCALARS ',cell_data(i).name,' ' type ' 1'];
    fprintf(fid,'%s\n',text);
    text = ['LOOKUP_TABLE default'];
    fprintf(fid,'%s\n',text);
    fprintf(fid,['%.' numsigdigs 'f\n'],cell_data(i).data');
end

fclose(fid);