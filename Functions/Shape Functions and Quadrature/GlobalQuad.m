function Quad = GlobalQuad(nsd, elem_type, qo)
%GLOBALQUAD Quadrature weights, points, and shape functions and 
%derivatives at the quadrature points. 
%   Quad = GLOBALQUAD(nsd, elem_type, qo) is a structure array with the
% 	quadrature weights, points, shape functions, and shape function 
% 	derivatives evaluated at each quadrature point
% 
%   --------------------------------------------------------------------
% 	Input
%   --------------------------------------------------------------------
%   nsd:    	number of spatial dimensions of the problem.
%   elem_type: 	the toplogical class of finite element; it is in the 
%           	general form 'topology-#of nodes' ie a three node 
%           	triangle is T3 a four node quadralateral is Q4 a 4 node 
%           	tetrahedra is H4 a 27 node brick is B27 etc. Presently 
%           	defined are L2, L3, L4, T3, T4(cubic bubble), T6, Q4, 
% 				Q9, H4, H10, B8 and B27.  
%   qo:  		quadrature order
% 
%   --------------------------------------------------------------------
% 	Output
%   --------------------------------------------------------------------
% 	Quad: 	Structure array with fields,
%       	.nq:     Number of quadrature points	
%       	.W:      Vector of quadrature weights (size nq x 1)      
%       	.Q:      Vector of quadrature points	(size nq x nsd)
%       	.Nq:     Cell array (size nq x 1) with shape functions  
%       	         evaluated at each quadrature point
%       	.dNdxiq: Cell array (size nq x 1) with derivative of shape 
%       	         functions w.r.t. parent coordinates evaluated at 
%       	         each quadrature point
%       	.Nv:     Cell array (size nq x 1) with shape functions 
%       	         evaluated at each quadrature point in Voigt form

% Acknowledgments: Matin Parchei Esfahani

% quadrature weights and points
[Quad.W, Quad.Q] = quadrature(qo, 'GAUSS', nsd);

% number of quad points
Quad.nq = length(Quad.W); 

% cell array of shape functions at each quadrature point
Quad.Nq = cell(Quad.nq,1);

% cell array of derivative of shape functions at each quadrature point 
Quad.dNdxiq = cell(Quad.nq,1);
Quad.Nv = cell(Quad.nq,1);

for q = 1:Quad.nq
    % quadrature point in parent coordinate
    xi = Quad.Q(q,:);       
    % N: shape function evaluated at quadrature point
    % dNdxi: derivative of shape function wrt parent coordinate, 
    %          xi, evaluated at quadrature point
    [Quad.Nq{q}, Quad.dNdxiq{q}] = lagrange_basis(elem_type, xi, nsd);  
    Quad.Nv{q} = getNv(Quad.Nq{q}, nsd);
end

end