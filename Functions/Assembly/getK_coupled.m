function K = getK_coupled(Mesh, Quad, Material)
%GETK Stiffness matrix 
%   K = GETK_COUPLED(Mesh, Quad, Material) is the global thermoelasticity
%   stiffness matrix for parameters defined in the structure arrays Mesh, 
%   Quad, and Material. The sparse matrix has size Mesh.nDOF x Mesh.nDOF
%   
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   Mesh:       Structure array with the following fields,
%               .ne:    Total number of elements in the mesh
%               .nne:   Vector of number of nodes per element (size nsd x 1)
%               .nsd:   Number of spatial dimensions
%               .conn:  Array of element connectivity (size ne x nne)
%               .x:     Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)
%               .DOF:   Array of DOF indices (size nn x nsd)
%               .nDOFe: Number of DOFs per element
%               .nDOF:  Total number of DOFs
%  
%   Quad:       Structure array with the following fields,
%               .W:      Vector of quadrature weights (size nq x 1)      
%               .nq:     Number of quadrature points 
%               .Nq:     Cell array (size nq x 1) with shape functions  
%                        evaluated at each quadrature point
%               .dNdxiq: Cell array (size nq x 1) with derivative of shape 
%                        functions w.r.t. parent coordinates evaluated at 
%                        each quadrature point
% 
%   Material:   Structure array with the following fields,
%               .t:         Material thickness

% Acknowledgements: Chris Ladubec

% initialize stiffness matrix
Kuuvec_size = Mesh.ne*(Mesh.nne * Mesh.nsd)^2;          % vector size
Kttvec_size = Mesh.ne*(Mesh.nne)^2;                     % vector size
Kutvec_size = Mesh.ne*Mesh.nDOFe_mech*Mesh.nDOFe_sca;   % vector size

rowu = zeros(Kuuvec_size, 1);                   % vector of row indices
colu = zeros(Kuuvec_size, 1);                   % vector of column indices

rowt = zeros(Kttvec_size, 1);                   % vector of row indices
colt = zeros(Kttvec_size, 1);                   % vector of column indices

rowut = zeros(Kutvec_size, 1);                   % vector of row indices
colut = zeros(Kutvec_size, 1);                   % vector of column indices

Kuuvec = zeros(Kuuvec_size, 1);                  % vectorized stiffness matrix
Kttvec = zeros(Kttvec_size, 1);                  % vectorized stiffness matrix
Kutvec = zeros(Kutvec_size, 1);                  % vectorized stiffness matrix

count_u = 1;                                  % DOF counter
count_t = 1;                                  % DOF counter
count_ut = 1;                                 % DOF counter

% for each element, compute element stiffness matrix and add to global
for e = 1:Mesh.ne

    %% Element variables
        % nodal ids of the element's nodes
        enodes = Mesh.conn(e,:);    
        % global coordinates of the element's nodes
        xI = Mesh.x(enodes,:);  
        % number of degrees of freedom displacement field
        ndofE_mech = Mesh.nne*Mesh.nsd;
        % number of degrees of freedom temperature field
        ndofE_sca = Mesh.nne;
        % displacement DOFs (mechanical problem)
        dofE_mech = Mesh.DOF_mech(enodes,:);
        dofE_mech = reshape(dofE_mech', Mesh.nDOFe_mech,[]);
        % temperature DOFs (scalar problem)
        dofE_sca = Mesh.DOF_sca(enodes,:);
        dofE_sca = reshape(dofE_sca', Mesh.nDOFe_sca,[]);
        
    %% Constitutive matrix
        nMat = Mesh.MatList(e); % element material type
        % mechanic problem
        Material.Dtype = Material.Dtype_mech;
        D_mech = getD(nMat, Material, Mesh);
        % thermal problem
        Material.Dtype = Material.Dtype_therm;
        D_sca = getD_TH1(nMat, Material, Mesh);
        %coupled matrix coefficient
        beta = Material.Prop(nMat).alpha*Material.Prop(nMat).E0/(1-2*Material.Prop(nMat).nu);
        
    %% Shape functions and derivatives in parent coordinates
        W = Quad.W;
        nq = Quad.nq;

    %% Assemble stiffness matrix
    
        % length of element. used to check that quadrature points and weights 
        % are correct. A = Sum Wi*Ji
        A = 0;

        % initialize element stiffness matrix
        Kuue = zeros(ndofE_mech, ndofE_mech); 
        Ktte = zeros(ndofE_sca, ndofE_sca);
        Kute = zeros(ndofE_mech, ndofE_sca);
        
        % loop through all quadrature points
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
            if dJe < 0
               error('Element %d has a negative Jacobian.', e)
            end

            % derivative of shape function in physical coordinates 
            % (tensor form)
            dNdxi = dNdxi';
            B = Je\dNdxi;

            % convert N matric to a Voigt form
            Nv = getNv(N, Mesh.nsd);
            
            % convert B matrix to Voigt form
            Bv = getBv(B', Mesh.nsd);
            
            % for 2D, volume integral includes the thickness
            switch Mesh.nsd 
                case 1
                    L = Material.t(Xi);
                    m = 1;
                case 2
                    L = Material.t(Xi);
                    m = [1;1;0]; % mapping vector
                case 3
                    L = 1;                
            end

            % Calculate local stiffness matrix
            Kuue = Kuue + W(q)*Bv*D_mech*Bv'*L*dJe;
            
            % Calculate local conductivity matrix
            Ktte = Ktte + W(q)*B'*D_sca*B*L*dJe;
            
            % Calculate local coupled matrix
            Kute = Kute + W(q)*Bv*m*beta*N'*L*dJe;
            
            % quadrature debug tool
            A = A + W(q)*dJe; 
        end
   
    %% Forming the vectorized stiffness matrix
    % elasticity stiffness matrix
        count_u = count_u + ndofE_mech^2;
        Kuue = reshape(Kuue, [ndofE_mech^2, 1]);
        rowmatrix_u = dofE_mech*ones(1,ndofE_mech);        
        rowe_u = reshape(rowmatrix_u,[ndofE_mech^2, 1]);
        cole_u = reshape(rowmatrix_u',[ndofE_mech^2, 1]);
        
        Kuuvec(count_u-ndofE_mech^2:count_u-1) = Kuue;
        rowu(count_u-ndofE_mech^2:count_u-1) = rowe_u;
        colu(count_u-ndofE_mech^2:count_u-1) = cole_u;
        
    % conductivity stiffness matrix
        count_t = count_t + ndofE_sca^2;
        Ktte = reshape(Ktte, [ndofE_sca^2, 1]);
        rowmatrix_t = dofE_sca*ones(1,ndofE_sca);        
        rowe_t = reshape(rowmatrix_t,[ndofE_sca^2, 1]);
        cole_t = reshape(rowmatrix_t',[ndofE_sca^2, 1]);
        
        Kttvec(count_t-ndofE_sca^2:count_t-1) = Ktte;
        rowt(count_t-ndofE_sca^2:count_t-1) = rowe_t;
        colt(count_t-ndofE_sca^2:count_t-1) = cole_t;
        
    % coupled stiffness matrix
        count_ut = count_ut + ndofE_mech*ndofE_sca;
        Kute = reshape(Kute, [ndofE_mech*ndofE_sca, 1]);       
        rowe_ut = reshape(dofE_mech*ones(1,ndofE_sca),[ndofE_mech*ndofE_sca, 1]);
        cole_ut = reshape(ones(ndofE_mech,1)*dofE_sca',[ndofE_mech*ndofE_sca, 1]);
        
        Kutvec(count_ut-ndofE_mech*ndofE_sca:count_ut-1) = Kute;
        rowut(count_ut-ndofE_mech*ndofE_sca:count_ut-1) = rowe_ut;
        colut(count_ut-ndofE_mech*ndofE_sca:count_ut-1) = cole_ut;
end

% sparse elasticity stiffness matrix
Kuu = sparse(rowu, colu, Kuuvec, Mesh.nDOF, Mesh.nDOF);

% sparse conductivity stiffness matrix
Ktt = sparse(rowt, colt, Kttvec, Mesh.nDOF, Mesh.nDOF);

% sparse coupling matrix 
Kut = sparse(rowut, colut, Kutvec, Mesh.nDOF, Mesh.nDOF);

% assemble global coupled matrix
K = Kuu + Ktt - Kut;

end