function xi_node = getXI( x_node,type )
%GETXI returns the parent coordinates of the nodal position specified.
%   ----------------------------------------------------------------------
%   Created by Endrina Rivas
%       endrina.rivas@uwaterloo.ca
%       Department of Civil Engineering
%       University of Waterloo
%       November 2015
%   ----------------------------------------------------------------------

switch type
    case 'L2'
        if x_node == 1
            xi_node = -1;
        elseif x_node == 2
            xi_node = 1; 
        end
    case 'L3'
        if x_node == 1
            xi_node = -1;
        elseif x_node == 2
            xi_node = 0;
        elseif x_node == 3
            xi_node = 1;
        end
    case 'L4'
        if x_node == 1
            xi_node = -1;
        elseif x_node == 2
            xi_node = -1/3;
        elseif x_node == 3
            xi_node = 1/3;
        elseif x_node == 4
            xi_node = 1;
        end
    case 'T3'
        if x_node == 1
        elseif x_node == 2
        elseif x_node == 3
        end
    case 'Q4'
        if x_node == 1
            xi_node = [-1 -1];
        elseif x_node == 2
            xi_node = [1 -1];
        elseif x_node == 3
            xi_node = [1 1];
        elseif x_node == 4
            xi_node = [-1 1];
        end
        
    case 'Q9'
        if x_node == 1
            xi_node = [-1 -1];
        elseif x_node == 2
            xi_node = [1 -1];
        elseif x_node == 3
            xi_node = [1 1];
        elseif x_node == 4
            xi_node = [-1 1];
        elseif x_node == 5
            xi_node = [0 -1];
        elseif x_node == 6
            xi_node = [1 0];
        elseif x_node == 7
            xi_node = [0 1];
        elseif x_node == 8
            xi_node = [-1 0];
        elseif x_node == 9
            xi_node = [0 0];
        end
    case 'B8'
        if x_node == 1
            xi_node = [-1 -1 -1];
        elseif x_node == 2
            xi_node = [1 -1 -1];
        elseif x_node == 3
            xi_node = [1 1 -1];
        elseif x_node == 4
            xi_node = [-1 1 -1];
        elseif x_node == 5
            xi_node = [-1 -1 1];
        elseif x_node == 6
            xi_node = [1 -1 1];
        elseif x_node == 7
            xi_node = [1 1 1];
        elseif x_node == 8
            xi_node = [-1 1 1];
        end
end

end

