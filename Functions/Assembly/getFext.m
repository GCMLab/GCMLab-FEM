function F = getFext(Mesh, BC, Quad)
%GETFEXT External forces
%   F = GETFEXT(Mesh, BC, Quad) is a column vector of external forces 
%   acting on each degree of freedom (size ndof x 1 in which ndof is the
%   number of degrees of freedom)
%   
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   Mesh:   Structure array with the following fields,
%           .nDOF:      Total number of DOFs
%           .ne:        Total number of elements in the mesh
%           .conn:      Array of element connectivity (size ne x nne)
%           .x:         Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)
%           .DOF:       Array of DOF indices (size nn x nsd)
%           .nDOFe:     Number of DOFs per element
%           .type:      the toplogical class of finite element; it is in 
%                       the general form 'topology-#of nodes' ie a three 
%                       node triangle is T3 a four node quadralateral is 
%                       Q4 a 4 node tetrahedra is H4 a 27 node brick is 
%                       B27 etc. Presently defined are L2, Q4, and Q9. 
%           .nsd:       Number of spatial dimensions
% 
%   BC:     Structure array with the following fields,
%           .b                          Anonymous function of distributed
%                                       body force (size 1 x nsd)
%           .traction_force_node        Column vector of nodes with 
%                                       prescribed tractions
%           .traction_force_value       Column vector of prescribed tractions
%                                       on nodes
%           .traction_force_dof         Column vector of degrees of freedom
%                                       with prescribed tractions
%           .traction_force_dof_value   Column vector of prescribed tractions
%                                       on DOF
% 
%   Quad:   Structure array with the following fields,
%           .W:         Vector of quadrature weights (size nq x 1)      
%           .Q:         Vector of quadrature points (size nq x nsd)
%           .nq:        Number of quadrature points 
%           .Nq:        Cell array (size nq x 1) with shape functions  
%                       evaluated at each quadrature point
%           .dNdxiq:    Cell array (size nq x 1) with derivative of shape 
%                       functions w.r.t. parent coordinates evaluated at 
%                       each quadrature point

% Acknowledgements: Chris Ladubec

%% Initialize values

	% initialize source (body) force vector
	F = zeros(Mesh.nDOF, 1);       

    if ~isempty(BC.traction_force_node)
        % counter to assess whether point load has been yet accounted for
        t_count = zeros(size(BC.traction_force_node));
        traction_value = BC.traction_force_value; 
        traction_force_node = BC.traction_force_node;
    else
        traction_force_node = [];
    end

%% return zeros if there is no applied force
	if isempty(BC.traction_force_node) ...
	    && strcmp(func2str(BC.b),'@(x)[]') ...
	    && isempty(BC.traction_force_dof)
	    return
	end

%% loop through all elements

for e = 1:Mesh.ne

    %% Element variables
	    % nodal ids of the element's nodes
	    enodes = Mesh.conn(e,:);    
	    % global coordinates of the element's nodes
	    xI = Mesh.x(enodes,:);  
	    % DOFs of element nodes
        dofE = Mesh.DOF(enodes,:);
        dofE = reshape(dofE',Mesh.nDOFe,[]);
        % number of degrees of freedom
        ndofE = length(dofE);

    %% Shape functions and derivatives in parent coordinates
        W = Quad.W;
        Q = Quad.Q;
        nq = Quad.nq;

    %% Calculate element body force

	    % initialize element body force vector
	    Fbe = zeros(ndofE, 1);

        % length of element. used to check that quadrature points and weights 
        % are correct. A = Sum Wi*Ji
        A = 0;

	    % if there is a body force run through the quadrature loop
	    if ~strcmp(func2str(BC.b),'@(x)[]')
            for q = 1:nq

                % Shape functions and derivatives in parent coordinates
                N = Quad.Nq{q};
                dNdxi = Quad.dNdxiq{q};

                % quadrature point in physical coordinates
                Xi = xI'*N;

                % Jacobian of the transformation between parent and global 
                % coordinates
                Je = dNdxi'*xI;
               
                % determinant of the Jacobian
                dJe = det(Je);

                % Applied body force
                Fbe = Fbe + W(q)*Nv*BC.b(Xi,t)*dJe;

                % quadrature debug tool
                A = A + W(q)*dJe; 
            end
	    end
    
    %% Calculate element traction force vector
	    % initialize element traction force vector
	    Fte = zeros(ndofE, 1);   

	    % loop through all tractions
        for i = 1:length(traction_force_node)  

	        % check whether traction is applied to the element
	        if ~isempty(find(enodes == BC.traction_force_node(i), 1)) ...
	                && t_count(i) == 0 

	            % local node number at traction location
	            t_node = find(enodes == BC.traction_force_node(i),1); 
	            % parent coordinates of traction location
	            xi = getXI(t_node,Mesh.type);
	            % Shape function at traction location
                [N,dNdxi] = lagrange_basis(Mesh.type, xi, Mesh.nsd);
                Nv = getNv(N, Mesh.nsd);

	            Fte = Fte + Nv*traction_value(i,:)';

	            t_count(i) = 1;
	        end
    	end

    %% Assemble element forces
    	F(dofE) = F(dofE) + Fbe + Fte;
end

%% Add forces prescribed at dofs
    F(BC.traction_force_dof) = F(BC.traction_force_dof) ...
                                    + BC.traction_force_dof_value;

end