function nne = NodesPerElement( etype )
%NODESPERELEMENT 
%   
%   ----------------------------------------------------------------------
%   Created by Endrina Rivas
%       endrina.rivas@uwaterloo.ca
%       Department of Civil Engineering
%       University of Waterloo
%       July 2015
%   ----------------------------------------------------------------------

switch etype
    case 'pt'
        nne = 1;
    case 'L2'
        nne = 2;
    case 'L3'
        nne = 3;
    case 'L4'
        nne = 4;
    case 'Q4'
        nne = [2;2];
    case 'Q9'
        nne = [3;3];
    case 'B8'
        nne = [2;2;2];
end

end

