function Quad = GlobalQuad(nsd, elem_type, qo)

% Acknowledgments: Matin Parchei Esfahani

%% Quadrature for regular elements

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
        [Quad.Nq{q}, Quad.dNdxiq{q}] = lagrange_basis(elem_type, xi);  
        Quad.Nv{q} = getNv(Quad.Nq{q}, nsd);
    end

 end