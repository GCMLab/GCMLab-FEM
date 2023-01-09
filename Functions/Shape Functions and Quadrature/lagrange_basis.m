function [N, dNdxi, Nv] = lagrange_basis(type, coord, dim)
%LAGRANGE_BASIS Returns the lagrange interpolant basis and its gradients 
%with respect to the parent coordinate system
%   N = LAGRANGE_BASIS(type, coord) is a vector of shape functions
%   (size n x 1 where n is the number of nodes) evaluated at the 
%   specified parent coordinates (coord), for the specified element 
%   type (type).
%
%   [N, dNdxi] = LAGRANGE_BASIS(type, coord) also returns a vector of 
%   shape function derivatives (size n x dim where n is the number of 
%   nodes and dim is the number of model dimensions) evaluated at the 
%   specified parent coordinates (coord), for specified element type (type). 
% 
%   [N, dNdxi, Nv] = LAGRANGE_BASIS(type, coord) also returns a vector 
%   of shape functions in Voigt notation (size n*dim x dim where n is 
%   the number of nodes and dim is the number of spatial dimensions) 
%   evaluated at the specified parent coordinates (coord).
% 
%   --------------------------------------------------------------------
% 	Input
%   --------------------------------------------------------------------
%   type:   the topological class of finite element; it is in the general
%           form 'topology-#of nodes' ie a three node triangle is T3 a 
%           four node quadralateral is Q4 a 4 node tetrahedra is H4 a 27 
%           node brick is B27 etc. Presently defined are L2, L3, L4, T3, 
%           T4(cubic bubble), T6, Q4, Q9, H4, H10, B8 and B27.  
%   coord:  the parent coordinates at which the basis and its gradients are 
%           to be evaluated at. If enough coordinates are not specified, 
%           code will return error. If too many coordinates are specified, 
%           code will ignore the extra.
%   dim:    Spatial dimensions. The default is unity. 
%
%   --------------------------------------------------------------------
% 	Output
%   --------------------------------------------------------------------
%   N:      shape functions evaluated at the specified parent coordinates 
%           (coord). These are given as nx1 where n is the number of nodes
%   dNdxi:  derivatives of the shape functions evaluated at the specified 
%           parent coordinates (coord). These are given as n x dim, where 
%           n is the number of nodes and dim is the number of model 
%           dimensions.
%   Nv:     shape functions evaluated at the specified parent coordinates 
%           (coord). There are given in Voigt notation.
%
%   --------------------------------------------------------------------
%   References
%   --------------------------------------------------------------------
%   [1] Fish, J., & Belytschko, T. (2007). A first course in finite 
%       elements. Chichester; Hoboken, NJ: John Wiley & Sons.
%   [2] Liu, G. R., & Quek, S. S. (Eds.). (2003). Finite Element Method. 
%       Oxford: Butterworth-Heinemann. 
%   [3] Zienkiewicz, O. C., & Taylor, R. L. (2005). The Finite Element 
%       Method: Volume 1 (6th ed., Vol. 1). Butterworth-Heinemann. 
%   [4] Belytschko, T., Liu, W. K., Moran, B., & Elkhodary, K. (2014). 
%       Nonlinear Finite Elements for Continua and Structures. John Wiley &
%       Sons.
%   [5] J. N. Reddy  (2006). AN Introduction to the Finite Element Meethod.

% Acknowledgements: Jack Chessa

% -------------------------------------------------------------------
%% Set to one-dimension if dimension has not been specified
% -------------------------------------------------------------------
if nargin == 2
    dim = 1;
end

% -------------------------------------------------------------------
%% Define the element shape functions in 1D parent coordinates for
% selected element type
% -------------------------------------------------------------------    
switch type
    case 'L2'  
    %%%%%%%%%%%%%%%%%%% L2 TWO NODE LINE ELEMENT %%%%%%%%%%%%%%%%%%% 
    %  
    %    1---------2
    %
        if size(coord,2) < 1
            error('Error: coordinate needed for the L2 element')
        else
            xi = coord(1);
            N =([1-xi, 1+xi]/2)';   %Ref. [1] pg. 87
            dNdxi = [-1; 1]/2;
        end

    case 'L3' 
    %%%%%%%%%%%%%%%% L3 THREE NODE LINE ELEMENT %%%%%%%%%%%%%%%%%% 
    %  
    %    1---------2----------3
    %
        if size(coord,2) < 1
            error('Error: coordinate needed for the L3 element')
        else
            xi = coord(1);
            %N = [(1-xi)*xi/(-2); (1+xi)*xi/2; 1-xi^2]; 
            % Seems to be in the wrong order, but Ref. [2] pg. 88 also
            % has it this way.
            % Compare with Ref. [1] pg. 82
            N = [   xi*(xi-1)/2; 
                    1-xi^2; 
                    xi*(1+xi)/2];
            dNdxi = [xi-0.5; -2*xi; xi+0.5];
        end

    case 'L4'
    %%%%%%%%%%%%%%%% L4 FOUR NODE LINE ELEMENT %%%%%%%%%%%%%%%%%% 
    %  
    %    1---------2----------3----------4
    %
        if size(coord,2) < 1
            error('Error: coordinate needed for the L4 element')
        else
            xi = coord(1);
            N = [ -1/16*(1-xi)*(1-9*xi^2);  % Modified from Ref. [2] pg. 88 
                  9/16*(1-3*xi)*(1-xi^2);
                  9/16*(1+3*xi)*(1-xi^2);
                  -1/16*(1+xi)*(1-9*xi^2)];
            dNdxi = [-1/16*(27*xi^2-18*xi-1);
                     9/16*(9*xi^2-2*xi-3);
                     9/16*(-9*xi^2-2*xi+3);
                     -1/16*(-27*xi^2-18*xi+1)];
        end

    case 'T3'
    %%%%%%%%%%%%% T3 THREE NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%% 
    %   
    %               3
    %             /  \
    %            /    \
    %           /      \
    %          /        \
    %         /          \
    %        /            \
    %       /              \
    %      /                \
    %     /                  \
    %    1--------------------2
    %
        if size(coord,2) < 2
            error('Error: two coordinates needed for the T3 element')
        else
            xi = coord(1); 
            eta = coord(2);
            N = [1-xi-eta; xi; eta];
            dNdxi = [-1,-1; 1,0; 0,1];
        end

    case 'T4'
    %%%%%%% T4 FOUR NODE TRIANGULAR CUBIC BUBBLE ELEMENT %%%%%%%%% 
    %   
    %               3
    %             /  \
    %            /    \
    %           /      \
    %          /        \
    %         /          \
    %        /      4     \
    %       /              \
    %      /                \
    %     /                  \
    %    1--------------------2
    %
        if size(coord,2) < 2
            error('Error: two coordinates needed for the T4 element')
        else
            xi = coord(1); 
            eta = coord(2);
            N = [1-xi-eta-3*xi*eta; 
                xi*(1-3*eta); 
                eta*(1-3*xi); 
                9*xi*eta];
            dNdxi = [-1-3*eta, -1-3*xi;
                     1-3*eta, -3*xi;
                     -3*eta,   1-3*xi;
                     9*eta,   9*xi ];
        end

    case 'T6'
    %%%%%%%%%%%%%%% T6 SIX NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%
    %   
    %               3
    %             /  \
    %            /    \
    %           /      \
    %          /        \
    %         6          5
    %        /            \
    %       /              \
    %      /                \
    %     /                  \
    %    1---------4----------2
    %

        if size(coord,2) < 2
            error('Error: two coordinates needed for the T6 element')
        else
            xi = coord(1); 
            eta = coord(2);
            N = [1-3*(xi+eta)+4*xi*eta+2*(xi^2+eta^2);
                                          xi*(2*xi-1);
                                        eta*(2*eta-1);
                                      4*xi*(1-xi-eta);
                                             4*xi*eta;
                                    4*eta*(1-xi-eta)];

            dNdxi = [4*(xi+eta)-3   4*(xi+eta)-3;
                       4*xi-1              0; 
                            0        4*eta-1;
                    4*(1-eta-2*xi)       -4*xi;
                        4*eta           4*xi;
                       -4*eta  4*(1-xi-2*eta)];
        end

    case 'Q4'
    %%%%%%%%%%%% Q4 FOUR NODE QUADRILATERAL ELEMENT %%%%%%%%%%%%%
    %
    %    4--------------------3
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    1--------------------2
    %
        if size(coord,2) < 2
            error('Error: two coordinates needed for the Q4 element')
        else
            xi = coord(1); 
            eta = coord(2);
            N = 1/4*[(1-xi)*(1-eta);        % Ref. [2] pg. 143
                    (1+xi)*(1-eta);         % Ref. [1] pg. 165
                    (1+xi)*(1+eta);
                    (1-xi)*(1+eta)];
            dNdxi = 1/4*[-(1-eta), -(1-xi);
                        1-eta,    -(1+xi);
                        1+eta,      1+xi;
                        -(1+eta),   1-xi];
        end
        
    case 'Q8'
    %%%%%%%%%%%% Q8 NINE NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%
    %
    %    4---------7----------3
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    8                    6
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    1----------5---------2   
        if size(coord,2) < 2
            error('Error: two coordinates needed for the Q8 element')
        else
            xi = coord(1); 
            eta = coord(2);
            % In this case, nodes are numbered as in % Ref. [5] pg. 537 
            % Ref [5] → sketch
            % 1 → 1
            % 2 → 5
            % 3 → 2
            % 4 → 8
            % 5 → 6
            % 6 → 4
            % 7 → 7
            % 8 → 3
%             N = [-1/4*(1 - xi)*(1 - eta)*(1 + xi - eta);        
%                     1/2*(-xi^2 + 1)*(1 - eta);         
%                     1/4*(1 + xi)*(1 - eta)*(-1 + xi - eta);
%                     1/2*(1 - xi)*(-eta^2 + 1);
%                     1/2*(1 + xi)*(-eta^2 + 1);
%                     1/4*(1 - xi)*(1 + eta)*(-1 - xi + eta);
%                     1/2*(-xi^2 + 1)*(1 + eta);
%                     1/4*(1 + xi)*(1 + eta)*(-1 + xi + eta)];
%             dNdxi = [((1 - eta)*(1 + xi - eta))/4 - (1/4 - xi/4)*(1 - eta), (1/4 - xi/4)*(1 + xi - eta) + (1/4 - xi/4)*(1 - eta);        
%                     -xi*(1 - eta), xi^2/2 - 1/2;  
%                     ((1 - eta)*(-1 + xi - eta))/4 + ((1 + xi)*(1 - eta))/4, -((1 + xi)*(-1 + xi - eta))/4 - ((1 + xi)*(1 - eta))/4 ;  
%                     eta^2/2 - 1/2, -(1 - xi)*eta;  
%                     -eta^2/2 + 1/2, -(1 + xi)*eta ;  
%                     -((1 + eta)*(-1 - xi + eta))/4 - ((1 - xi)*(1 + eta))/4, ((1 - xi)*(-1 - xi + eta))/4 + ((1 - xi)*(1 + eta))/4;  
%                     -xi*(1 + eta), -xi^2/2 + 1/2;  
%                     ((1 + eta)*(-1 + xi + eta))/4 + ((1 + xi)*(1 + eta))/4,  ((1 + xi)*(-1 + xi + eta))/4 + ((1 + xi)*(1 + eta))/4 ];
            % In this case, nodes are numbered as in the sketch 
            
            N = [-1/4*(1 - xi)*(1 - eta)*(1 + xi - eta); 
                    1/4*(1 + xi)*(1 - eta)*(-1 + xi - eta);
                    1/4*(1 + xi)*(1 + eta)*(-1 + xi + eta);
                    1/4*(1 - xi)*(1 + eta)*(-1 - xi + eta);
                    1/2*(-xi^2 + 1)*(1 - eta);         
                    1/2*(1 + xi)*(-eta^2 + 1);
                    1/2*(-xi^2 + 1)*(1 + eta);
                    1/2*(1 - xi)*(-eta^2 + 1)];
            dNdxi = [((1 - eta)*(1 + xi - eta))/4 - (1/4 - xi/4)*(1 - eta), (1/4 - xi/4)*(1 + xi - eta) + (1/4 - xi/4)*(1 - eta);        
                    ((1 - eta)*(-1 + xi - eta))/4 + ((1 + xi)*(1 - eta))/4, -((1 + xi)*(-1 + xi - eta))/4 - ((1 + xi)*(1 - eta))/4 ;  
                    ((1 + eta)*(-1 + xi + eta))/4 + ((1 + xi)*(1 + eta))/4,  ((1 + xi)*(-1 + xi + eta))/4 + ((1 + xi)*(1 + eta))/4;
                    -((1 + eta)*(-1 - xi + eta))/4 - ((1 - xi)*(1 + eta))/4, ((1 - xi)*(-1 - xi + eta))/4 + ((1 - xi)*(1 + eta))/4;
                    -xi*(1 - eta), xi^2/2 - 1/2;  
                    -eta^2/2 + 1/2, -(1 + xi)*eta ;
                    -xi*(1 + eta), -xi^2/2 + 1/2;
                    eta^2/2 - 1/2, -(1 - xi)*eta];      
        end
    

    case 'Q9'
    %%%%%%%%%%%% Q9 NINE NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%
    %
    %    4---------7----------3
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    8          9         6
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    1----------5---------2
    %
        if size(coord,2) < 2
            error('Error: two coordinates needed for the Q9 element')
        else
            xi = coord(1); 
            eta = coord(2);
            N = 1/4*[xi*eta*(xi-1)*(eta-1); % Ref. [2] pg. 157
                     xi*eta*(xi+1)*(eta-1); % Ref. [1] pg. 168
                     xi*eta*(xi+1)*(eta+1);
                     xi*eta*(xi-1)*(eta+1);
                     -2*eta*(xi+1)*(xi-1)*(eta-1);
                     -2*xi*(xi+1)*(eta+1)*(eta-1);
                     -2*eta*(xi+1)*(xi-1)*(eta+1);
                     -2*xi*(xi-1)*(eta+1)*(eta-1);
                     4*(xi+1)*(xi-1)*(eta+1)*(eta-1)];
          dNdxi = 1/4*[eta*(2*xi-1)*(eta-1),    xi*(xi-1)*(2*eta-1);
                        eta*(2*xi+1)*(eta-1),   xi*(xi+1)*(2*eta-1);
                        eta*(2*xi+1)*(eta+1),   xi*(xi+1)*(2*eta+1);
                        eta*(2*xi-1)*(eta+1),   xi*(xi-1)*(2*eta+1);
                        -4*xi*eta*(eta-1),  -2*(xi+1)*(xi-1)*(2*eta-1);
                        -2*(2*xi+1)*(eta+1)*(eta-1),  -4*xi*eta*(xi+1);
                        -4*xi*eta*(eta+1),  -2*(xi+1)*(xi-1)*(2*eta+1);
                        -2*(2*xi-1)*(eta+1)*(eta-1),  -4*xi*eta*(xi-1);
                        8*xi*(eta^2-1), 8*eta*(xi^2-1)];
        end

    case 'H4'
    %%%%%%%%%%%% H4 FOUR NODE TETRAHEDRAL ELEMENT %%%%%%%%%%%%%%
    %
    %             4
    %           / | \
    %          /  |  \
    %         /   |   \ 
    %        /    |    \ 
    %       /     |     \
    %      1 -----|------3
    %         -   2  -
        if size(coord,2) < 3
            error('Error: three coordinates needed for the H4 element')
        else % Ref.(3) pg. 165
            xi = coord(1); % L1
            eta = coord(2); % L2
            zeta = coord(3); % L3
            N = [1-xi-eta-zeta; % L4
                            xi;
                           eta;
                          zeta];
            dNdxi = [-1  -1  -1;
                      1   0   0;
                      0   1   0;
                      0   0   1];
        end

    case 'H10'
    %%%%%%%%%%%%% H10 TEN NODE TETRAHEDRAL ELEMENT %%%%%%%%%%%%%%%
    %
    %             4
    %           / | \
    %          /  |  \
    %         /   |   \ 
    %        7    |    9 
    %       /     10     \
    %      /      |      \
    %     /       |       \
    %    1- - - -6|- - - --3
    %      -      |      -
    %        5    |    8
    %          -  |  -
    %             2
        if size(coord,2) < 3
            error('Error: three coordinates needed for the H10 element')
        else
            xi = coord(1); 
            eta = coord(2); 
            zeta = coord(3);
            phi = [1-xi-eta-zeta;  xi;  eta;  zeta];
            N = [phi(1)*(2*phi(1)-1);
                phi(2)*(2*phi(2)-1);
                phi(3)*(2*phi(3)-1);
                phi(4)*(2*phi(4)-1);
                4*phi(1)*phi(2);
                4*phi(1)*phi(3);
                4*phi(1)*phi(4);
                4*phi(2)*phi(3);
                4*phi(3)*phi(4);
                4*phi(2)*phi(4)];
            dNdxi = 4*[-phi(1)+.25,   -phi(1)+.25,   -phi(1)+.25;
                        phi(2)-.25,    0,             0;
                        0,             phi(3)-.25,    0;
                        0,             0,             phi(4)-.25;
                        phi(1)-phi(2), -phi(2),       -phi(2);
                        -phi(3),       phi(1)-phi(3), -phi(3);
                        -phi(4),       -phi(4),       phi(1)-phi(4);
                        phi(3),        phi(2),        0;
                        0,             phi(4),        phi(3);
                        phi(4),        0,             phi(2) ];
        end

    case 'B8'
    %%%%%%%%%%%%%%% B8 EIGHT NODE BRICK ELEMENT %%%%%%%%%%%%%%%%
    % 
    %                  8 
    %               /    \    
    %            /          \
    %         /                \
    %      5                     \
    %      |\                     7
    %      |   \                / |
    %      |     \     4    /     |
    %      |        \    /        |
    %      |           6          |
    %      1           |          |
    %       \          |          3
    %          \       |        /
    %            \     |     /
    %               \  |  /
    %                  2
    %                
        if size(coord,2) < 3
            error('Error: three coordinates needed for the B8 element')
        else
            I1 = 1/2-coord/2;
            I2 = 1/2+coord/2;
            N = [I1(1)*I1(2)*I1(3);
                 I2(1)*I1(2)*I1(3);
                 I2(1)*I2(2)*I1(3);
                 I1(1)*I2(2)*I1(3);
                 I1(1)*I1(2)*I2(3);
                 I2(1)*I1(2)*I2(3);
                 I2(1)*I2(2)*I2(3);
                 I1(1)*I2(2)*I2(3)]; % Ref.[4] pg. 626
            dNdxi = [-I1(2)*I1(3)   -I1(1)*I1(3)    -I1(1)*I1(2);
                     I1(2)*I1(3)    -I2(1)*I1(3)    -I2(1)*I1(2);
                     I2(2)*I1(3)    I2(1)*I1(3)     -I2(1)*I2(2);
                     -I2(2)*I1(3)   I1(1)*I1(3)     -I1(1)*I2(2)
                     -I1(2)*I2(3)   -I1(1)*I2(3)    I1(1)*I1(2)
                     I1(2)*I2(3)    -I2(1)*I2(3)    I2(1)*I1(2)
                     I2(2)*I2(3)    I2(1)*I2(3)     I2(1)*I2(2)
                     -I2(2)*I2(3)   I1(1)*I2(3)     I1(1)*I2(2)]*0.5;
        end

    case 'B27'
    %%%%%%%%%%% B27 TWENTY SEVEN NODE BRICK ELEMENT %%%%%%%%%%%%%%
    %
    %                  19 
    %               /    \    
    %            20         \
    %         /                22
    %      21                    \
    %      |\         23          25
    %      |   \                / |
    %     12     24         26    |
    %      |        \    /        16
    %      |          27          |
    %      3    15     |     17   |
    %       \          |          7
    %          \       18       /
    %           6      |     8
    %               \  |  /
    %                  9
    %                
        if size(coord,2) < 3
            error('Error: three coordinates needed for the B27 element')
        else
            N = zeros(27,1);
            dNdxi = zeros(27,3);
            xi = coord(1); 
            eta = coord(2); 
            zeta = coord(3);
            c = 1;
            for i = 1:3
                [Ni,dNdxI] = L3at(zeta,i);
                for j = 1:3
                    [Nj,dNdxj] = L3at(eta,j);
                    for k = 1:3
                        [Nk,dNdxk] = L3at(xi,k);
                        N(c) = Ni*Nj*Nk;
                        dNdxi(c,1) = Ni*Nj*dNdxk;
                        dNdxi(c,2) = Ni*dNdxj*Nk;
                        dNdxi(c,3) = dNdxI*Nj*Nk;
                        c = c+1;
                    end
                end
            end       
        end

    otherwise 
    %%%%%%%%%%%%%% OTHER ELEMENT TYPES UNDEFINED %%%%%%%%%%%%%%%%%
        error(['Element ',type,' not yet supported'])
end

%% Convert to Voigt Notation

% Shape functions in Voigt Notation
Nv = getNv(N, dim);

%% Calculate shape functions for L3 elements

function [Ni,dNdxi] = L3at(xi,index)
    N = [(1-xi)*xi/(-2);
        1-xi^2;
        (1+xi)*xi/2];
    dNdx = [xi-0.5; -2*xi; xi+0.5];
    Ni = N(index);
    dNdxi = dNdx(index);
end

end