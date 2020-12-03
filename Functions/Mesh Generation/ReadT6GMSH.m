function [nodes,ncrds,elmConn,oneDelmConn]= ReadGMSHGeometryFile ()
%   ----------------------------------------------------------------------
%   Created by Ali Ghavidel
%       aghavidel@uwaterloo.ca
%       Department of Civil Engineering
%       University of Waterloo
%   Last Updated: Oct 2020
%% Inputs
% .msh format file from GMSH with T6 elements

%% Outputs
%nodes: a vector nn*3 (nodes number, x, y)
%ncrds: coordinates of nodes (x,y)
%elmConn: element connectivity (ne*6)
%oneDelmConn: 1D elements on the boundaries

%% Import data from text file.
% Script for importing data from the following text file:

%% Initialize variables.
filename = 'C:\A DRIVE D\UW\Git\GCMLab-FEM\Test Files\GeometryT6.msh';
delimiter = ' ';
startRow = 2;

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
THMMainGeometry = table;
THMMainGeometry.MeshFormat = cell2mat(raw(:, 1));
THMMainGeometry.VarName2 = cell2mat(raw(:, 2));
THMMainGeometry.VarName3 = cell2mat(raw(:, 3));
THMMainGeometry.VarName4 = cell2mat(raw(:, 4));
THMMainGeometry.VarName5 = cell2mat(raw(:, 5));
THMMainGeometry.VarName6 = cell2mat(raw(:, 6));
THMMainGeometry.VarName7 = cell2mat(raw(:, 7));
THMMainGeometry.VarName8 = cell2mat(raw(:, 8));
THMMainGeometry.VarName9 = cell2mat(raw(:, 9));
THMMainGeometry.VarName10 = cell2mat(raw(:, 10));
THMMainGeometry.VarName11 = cell2mat(raw(:, 11));
THMMainGeometry.VarName12 = cell2mat(raw(:, 12));
THMMainGeometry.VarName13 = cell2mat(raw(:, 13));
THMMainGeometry.VarName14 = cell2mat(raw(:, 14));
THMMainGeometry.VarName15 = cell2mat(raw(:, 15));

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;



%%

a = table2array(THMMainGeometry);

for i=1:(size(a,1))
    if ~isnan(a(i,:))
        RowAllNumber = i;
    end
end

a(1:RowAllNumber,:)=[];

nodescell = find (isnan(a(:,5))& a(:,4)==0);
nodes = sortrows(a(nodescell,1:3));
ncrds = nodes(:,2:3);

a(nodescell,:)=[];

elmconncell = find (~isnan(a(:,6))& ~isnan(a(:,7)));
elmConn = a(elmconncell,2:7);

a(elmconncell,:)=[];


oneDelm = find (~isnan(a(:,1)) & ~isnan(a(:,4)));
b = a(oneDelm,1:4);
kk = 0;
for k=1:size(b,1)
    if b(k,1)>9
        kk = kk+1;
        oneDelmConn(kk,:) = b (k,2:4);
    end
end



end % function