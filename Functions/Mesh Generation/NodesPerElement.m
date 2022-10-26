function nne = NodesPerElement( etype )
%NODESPERELEMENT Number of nodes per element type
%   nne = NODESPERELEMENT(etype) is the number of nodes for the specified 
%   element type (etype)
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   etype:      the toplogical class of finite element; it is in the 
%               general form 'topology-#of nodes' ie a three node 
%               triangle is T3 a four node quadralateral is Q4 a 4 node 
%               tetrahedra is H4 a 27 node brick is B27 etc. Presently 
%               defined are L2, L3, L4, T3, T6, Q4, Q9, and B8.  

switch etype
    case 'pt'
        nne = 1;
    case 'L2'
        nne = 2;
    case 'L3'
        nne = 3;
    case 'L4'
        nne = 4;
    case 'T3'
        nne = 3;
    case 'T6'
        nne = 6;
    case 'Q4'
        nne = [2;2];
    case 'T6'
        nne = 6;
    case 'Q9'
        nne = [3;3];
    case 'B8'
        nne = [2;2;2];
end

end

