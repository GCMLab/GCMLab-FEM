function xi_node = getXI(node_num, type)
%GETXI Parent coordinates of the nodal position specified.
%   xi_node = GETXI(node_num, type) returns the parent coordinates at the  
%   node position (node_num) for the specified element type (type). 
%   For multi-dimensional elements, the coordinates are given as a 
%   vector of size 1 x nsd where nsd is the number of spatial dimensions.
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   node_num:   node number, based on standard element numbering
%   type:       the topological class of finite element; it is in the 
%               general form 'topology-#of nodes' ie a three node 
%               triangle is T3 a four node quadralateral is Q4 a 4 node 
%               tetrahedra is H4 a 27 node brick is B27 etc. Presently 
%               defined are L2, L3, L4, T3, Q4, Q9, and B8 B27.  

switch type
    case 'L2'
        if node_num == 1
            xi_node = -1;
        elseif node_num == 2
            xi_node = 1; 
        end
    case 'L3'
        if node_num == 1
            xi_node = -1;
        elseif node_num == 2
            xi_node = 0;
        elseif node_num == 3
            xi_node = 1;
        end
    case 'L4'
        if node_num == 1
            xi_node = -1;
        elseif node_num == 2
            xi_node = -1/3;
        elseif node_num == 3
            xi_node = 1/3;
        elseif node_num == 4
            xi_node = 1;
        end
    case 'T3'
        if node_num == 1
        elseif node_num == 2
        elseif node_num == 3
        end
    case 'Q4'
        if node_num == 1
            xi_node = [-1 -1];
        elseif node_num == 2
            xi_node = [1 -1];
        elseif node_num == 3
            xi_node = [1 1];
        elseif node_num == 4
            xi_node = [-1 1];
        end
    case 'Q9'
        if node_num == 1
            xi_node = [-1 -1];
        elseif node_num == 2
            xi_node = [1 -1];
        elseif node_num == 3
            xi_node = [1 1];
        elseif node_num == 4
            xi_node = [-1 1];
        elseif node_num == 5
            xi_node = [0 -1];
        elseif node_num == 6
            xi_node = [1 0];
        elseif node_num == 7
            xi_node = [0 1];
        elseif node_num == 8
            xi_node = [-1 0];
        elseif node_num == 9
            xi_node = [0 0];
        end
    case 'B8'
        if node_num == 1
            xi_node = [-1 -1 -1];
        elseif node_num == 2
            xi_node = [1 -1 -1];
        elseif node_num == 3
            xi_node = [1 1 -1];
        elseif node_num == 4
            xi_node = [-1 1 -1];
        elseif node_num == 5
            xi_node = [-1 -1 1];
        elseif node_num == 6
            xi_node = [1 -1 1];
        elseif node_num == 7
            xi_node = [1 1 1];
        elseif node_num == 8
            xi_node = [-1 1 1];
        end
end

end

