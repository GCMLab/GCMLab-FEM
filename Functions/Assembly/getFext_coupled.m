function F = getFext_coupled(Mesh, BC, Quad, t, Quad_edge)
%GETFEXT_COUPLED External forces for coupled problem
%   F = GETFEXT_COUPLED(Mesh, BC, Quad) is a column vector of external forces
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
%   
%   Quad_edge: Structure array with same fields as Quad
%               Used computation of tractions with edge elements

% Acknowledgements: Chris Ladubec

%% Initialize values

% Check input arguments
if nargin<4
    t = 0; % Set default time if not specified
end

% initialize source (body) force vector
F = zeros(Mesh.nDOF, 1);

if ~isempty(BC.traction_force_node)
    % counter to assess whether point load has been yet accounted for
    t_count = zeros(size(BC.traction_force_node));
    traction_force_node = BC.traction_force_node;
    if isa(BC.traction_force_value,'function_handle')
        traction_value = BC.traction_force_value(t);
    else
        traction_value = BC.traction_force_value;
    end
else
    traction_force_node = [];
end

if ~isempty(BC.flux_node)
    % counter to assess whether point flux has been yet accounted for
    t_count = zeros(size(BC.flux_node));
    flux_node = BC.flux_node;
    if isa(BC.flux_value,'function_handle')
        flux_value = BC.flux_value(t);
    else
        flux_value = BC.flux_value;
    end
else
    flux_node = [];
end

%% return zeros if there is no applied force or flux
if isempty(BC.traction_force_node) ...
        && strcmp(func2str(BC.b),'@(x)[]') ...
        && isempty(BC.traction_force_dof) ...
        && isempty(BC.flux_node) ...
        && strcmp(func2str(BC.s),'@(x)[]') ...
        && isempty(BC.flux_dof) ...
        && ~isa(BC.c_N_t_f,'cell')
    return
end

%% loop through all elements
% number of degrees of freedom displacement field
ndofE_mech = Mesh.nne*Mesh.nsd;
% number of degrees of freedom temperature field
ndofE_sca = Mesh.nne;

for e = 1:Mesh.ne
    
    %% Element variables
    % nodal ids of the element's nodes
    enodes = Mesh.conn(e,:);
    % global coordinates of the element's nodes
    xI = Mesh.x(enodes,:);
    % displacement DOFs (mechanical problem)
    dofE_mech = Mesh.DOF_mech(enodes,:);
    dofE_mech = reshape(dofE_mech', Mesh.nDOFe_mech,[]);
    % temperature DOFs (scalar problem)
    dofE_sca = Mesh.DOF_sca(enodes,:);
    dofE_sca = reshape(dofE_sca', Mesh.nDOFe_sca,[]);
    
    %% Shape functions and derivatives in parent coordinates
    W = Quad.W;
    Q = Quad.Q;
    nq = Quad.nq;
    
    %% Calculate element body force
    
    % initialize element body force vector
    Fbe = zeros(ndofE_mech, 1);
    
    % initialize element flux source vector
    Fse = zeros(ndofE_sca, 1);
    
    % length of element. used to check that quadrature points and weights
    % are correct. A = Sum Wi*Ji
    A = 0;
    
    % if there is a body force run through the quadrature loop
    if ~strcmp(func2str(BC.b),'@(x,t)[]') && ~strcmp(func2str(BC.b),'@(x)[]')
        for q = 1:nq
            
            % Shape functions and derivatives in parent coordinates
            N = Quad.Nq{q};
            dNdxi = Quad.dNdxiq{q};
            Nv = getNv(N, Mesh.nsd);
            
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
    
    % if there is a flux source run through the quadrature loop
    if ~strcmp(func2str(BC.s),'@(x,t)[]') && ~strcmp(func2str(BC.s),'@(x)[]')
        for q = 1:nq
            
            % Shape functions and derivatives in parent coordinates
            N = Quad.Nq{q};
            dNdxi = Quad.dNdxiq{q};
            Nv = getNv(N, Mesh.nsd);
            
            % quadrature point in physical coordinates
            Xi = xI'*N;
            
            % Jacobian of the transformation between parent and global
            % coordinates
            Je = dNdxi'*xI;
            
            % determinant of the Jacobian
            dJe = det(Je);
            
            % Applied body force
            Fse = Fse + W(q)*N*BC.s(Xi,t)*dJe;
            
            % quadrature debug tool
            A = A + W(q)*dJe;
        end
    end
        
    %% Calculate element traction force vector
    % initialize element traction force vector
    % Fte = zeros(ndofE_mech, 1);
    
    % initialize element traction force vector
    % Ffe = zeros(ndofE_sca, 1);
    %
    %----  Saving code for possible future adaption to edge element
    %integration-----
    %
    % 	    % loop through all tractions
    %         for i = 1:length(traction_force_node)
    %
    % 	        % check whether traction is applied to the element
    % 	        if ~isempty(find(enodes == BC.traction_force_node(i), 1)) ...
    % 	                && t_count(i) == 0
    %
    % 	            % local node number at traction location
    % 	            t_node = find(enodes == BC.traction_force_node(i),1);
    % 	            % parent coordinates of traction location
    % 	            xi = getXI(t_node,Mesh.type);
    % 	            % Shape function at traction location
    %                 [N,dNdxi] = lagrange_basis(Mesh.type, xi, Mesh.nsd);
    %                 Nv = getNv(N, Mesh.nsd);
    %
    % 	            Fte = Fte + Nv*traction_value(i,:)';
    %
    % 	            t_count(i) = 1;
    % 	        end
    %     	end
    
    %% Assemble element forces
    F(dofE_mech) = F(dofE_mech) + Fbe;% + Fte;
    F(dofE_sca) = F(dofE_sca) + Fse;% + Ffe; 

end

%% Loop through boundary elements 

%--------------------------------------------------------------------------------------------------
% Note: this algorithm works only for straight boundary elements. If
% curvature edges were included the integration along the curvature line
% would have to be implemented
%--------------------------------------------------------------------------------------------------

if (isa(BC.c_N_t_f,'cell')) && ( strcmp(Mesh.type, 'Q4') || strcmp(Mesh.type, 'T3') )

    % Shape functions and derivatives in parent coordinates
    W = Quad_edge.W;
    Q = Quad_edge.Q;
    nq = Quad_edge.nq;

    % Auxiliary variables to account for the cases of a equilibrium or a
    % diffusion problem in a single code
        temp_1 = cell(1,2);
        temp_1{1} = [0,0];
        temp_1{2} = 1;
        temp_2 = cell(1,2);
        temp_2{1} = 1;
        temp_2{2} = zeros(2,2);
    
    % Loop through all sets
    for set_i = 1:length(Mesh.c_BC_N_t)
        
        t_f = BC.c_N_t_f{set_i};
        set_BC_N_t = Mesh.c_BC_N_t{set_i}; 
        set_BC_N_t_n_m = Mesh.c_BC_N_t_n_m{set_i}; % Normals
        set_BC_N_t_t_m = Mesh.c_BC_N_t_t_m{set_i}; % Tangents
        
        % Define if this set_i is applied on 
        BC_case = BC.c_N_t_flag_case{set_i}; 
        if (BC_case-1)
            DOF_set = Mesh.DOF_mech;
        else
            DOF_set = Mesh.DOF_sca;
        end

        % Loop through all edge elements
        for e_t = 1:length(Mesh.c_BC_N_t{set_i})

            % nodal ids of the element's nodes
            enodes = set_BC_N_t(e_t,:);    
            % global coordinates of the element's nodes
            xI = Mesh.x(enodes,:);  
            % local coordinates in the tangental direction
            l = ((xI(1,1)-xI(2,1))^2 + (xI(1,2)-xI(2,2))^2)^0.5;
            % local coordinates in the tangental direction
            sI = [0; l];
            % degrees of freedom
            dofE = DOF_set(enodes,:);
            dofE = reshape(dofE',BC_case*2,[]);
            % number of degrees of freedom
            ndofE = length(dofE);
            % Initialize force vector
            Fte_n = zeros(ndofE,1);
            
            %--------------------
            % Get normal vector
            n_e = set_BC_N_t_n_m(e_t,:);
            % Get transversal vector
            t_e = set_BC_N_t_t_m(e_t,:);
            
            % Loop through all quadrature points
            for q = 1:nq

                % Shape functions and derivatives in parent coordinates
                N = Quad_edge.Nq{q};
                dNdxi = Quad_edge.dNdxiq{q};

                Nv = getNv(N, Mesh.nDOFn);
                
                % quadrature point in physical coordinates in the
                % global coordinate system
                Xi = xI'*N;
                
                % quadrature point in physical coordinates in the
                % local system of the element 
                Si = sI'*N;

                % Jacobian of the transformation between parent and global 
                % coordinates
                Je = dNdxi'*sI;

                % determinant of the Jacobian
                dJe = det(Je);

                % Get tractions evaluated at gauss point
                t_f_g = [(t_f(Xi,t))'*n_e'; (t_f(Xi,t))'*t_e'].*...
                        (1 - BC.c_N_t_flag(set_i)) + ... % when tractions are given in  x and y directions
                        t_f(Xi,t).*BC.c_N_t_flag(set_i); % when tractions are given in n and s directions               
                % Account for equilibrium (vector) or diffusion (scalar) problem
                t_f_g = t_f_g(1:Mesh.nDOFn);

                % Applied traction force in normal and tangential
                % directions
                Fte_n = Fte_n + W(q)*Nv*t_f_g*dJe;

            end
            
            % Transform local force vector to global system and add to
            % total force vector
            % rotation(transformation) matrix
            r_m = temp_1{Mesh.nDOFn}*[n_e', t_e']*temp_1{Mesh.nDOFn}' + temp_2{Mesh.nDOFn};
            T = [r_m, zeros(Mesh.nDOFn,Mesh.nDOFn); zeros(Mesh.nDOFn,Mesh.nDOFn), r_m];
            Fte = T*Fte_n;

            % Assemble tractions
            F(dofE) = F(dofE) + Fte;

        end
    end
end

%% Apply pre-integrated boundary node point forces
if ~isempty(BC.traction_force_node)
    switch Mesh.nDOFn_mech
        case 1
            F(BC.traction_force_node) = F(BC.traction_force_node) + traction_value(:,1);
        case 2
            traction_xdofs = 2*BC.traction_force_node-1;
            traction_ydofs = 2*BC.traction_force_node;
            
            F(traction_xdofs) = F(traction_xdofs) + traction_value(:,1);
            F(traction_ydofs) = F(traction_ydofs) + traction_value(:,2);
        case 3
            traction_xdofs = 3*BC.traction_force_node-2;
            traction_ydofs = 3*BC.traction_force_node-1;
            traction_zdofs = 3*BC.traction_force_node;
            
            F(traction_xdofs) = F(traction_xdofs) + traction_value(:,1);
            F(traction_ydofs) = F(traction_ydofs) + traction_value(:,2);
            F(traction_zdofs) = F(traction_zdofs) + traction_value(:,3);
    end
end

%% Apply pre-integrated boundary node point fluxes
if ~isempty(BC.flux_node)
    flux_dofs = 3*BC.flux_node;
    
    F(flux_dofs) = F(flux_dofs) + flux_value(:,1);
end

%% Add forces prescribed at dofs
F(BC.traction_force_dof) = F(BC.traction_force_dof) ...
    + BC.traction_force_dof_value;

%% Add fluxes prescribed at dofs
F(BC.flux_dof) = F(BC.flux_dof) ...
    + BC.flux_dof_value;

end