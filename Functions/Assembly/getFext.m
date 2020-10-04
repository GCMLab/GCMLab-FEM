function F = getFext(Mesh, BC, Material, Quad, Control)
%GETF
%   F = getF(t)
%
%   --------------------------------------------------------------------
%   Created by Chris Ladubec
%
%   Modified by Endrina Rivas
%       endrina.rivas@uwaterloo.ca
%       Department of Civil Engineering
%       University of Waterloo
%       
%   Last Updated June 2016
%   --------------------------------------------------------------------

%% initialize values

	% initialize source (body) force vector
	F = zeros(Mesh.nDOF, 1);       

    if isfield(BC,'traction_force_node')
        % counter to assess whether point load has been yet accounted for
        t_count = zeros(size(BC.traction_force_node));
        traction_value = BC.traction_force_value; 
        traction_force_node = BC.traction_force_node;
    else
        traction_force_node = [];
    end

%% return zeros if there is no applied force
	if ( ~isfield(BC,'traction_force_node') || isempty(BC.traction_force_node) )...
	        && strcmp(func2str(BC.b),'@(x)[]') ...
	        && ~isfield(BC,'traction_force_dof')
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
                if ~strcmp(func2str(BC.b),'@(x)[]')
                    Fbe = Fbe + W(q)*Nv*BC.b(Xi,t)*dJe;
                end

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
                [N,dNdxi] = lagrange_basis(Mesh.type, xi);
                Nv = getNv(N, Mesh.nsd);

	            Fte = Fte + Nv*traction_value(i,:)';

	            t_count(i) = 1;
	        end
    	end

    %% Assemble element forces
    	F(dofE) = F(dofE) + Fbe + Fte;
end

%% Add forces prescribed at dofs
    if isfield(BC,'traction_force_dof')
        F(BC.traction_force_dof) = F(BC.traction_force_dof) ...
                                    + BC.traction_force_dof_value;
    end

end