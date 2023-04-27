function [strain, stress] = getStrain_NLelastic(d, Mesh, Material, calc_type, Quad, dnm1)
%GETSTRAIN Evaluate stress and strain
%   [strain, stress] = GETSTRAIN_NLELASTIC(d, Mesh, Material) returns two 
%   matrices of strains and stresses computed at the center of each element. 
%   The matrices are of size dim x ne, in which dim = 1 for 1D elements, 3 
%   for 2D elements, and 6 for 3D elements.
%   
%   The non linear elastic material is defined by the law:
%   E = E0 + E1*I1^2
%
%   [strain, stress] = GETSTRAIN_NLELASTIC(d, Mesh, Material, 'none') does 
%   not compute the stresses or strains
%
%   [strain, stress] = GETSTRAIN_NLELASTIC(d, Mesh, Material, 'nodal') 
%   returns matrices of nodal-averaged strains (size dim x nn).
%
%   [strain, stress] = GETSTRAIN_NLELASTIC(d, Mesh, Material, 'center') 
%   returns matrices of strains computed at the center of each element
%   (size dim x ne).
%
%   [strain, stress] = GETSTRAIN_NLELASTIC(d, Mesh, Material, 'L2projection') 
%   returns matrices of L2-projected stresses at the nodes (size dim x nn).
%
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   d:          Vector of displacements at each DOF (size ndof x 1 in
%               which ndof is the number of degrees of freedom)
%
%   Mesh:       Structure array with the following fields,
%               .nsd:   Number of spatial dimensions
%               .ne:    Total number of elements in the mesh
%               .nn:    Total number of nodes
%               .nne:   Vector of number of nodes per element (size nsd x 1)
%               .type:  the toplogical class of finite element; it is in
%                       the general form 'topology-#of nodes' ie a three
%                       node triangle is T3 a four node quadralateral is
%                       Q4 a 4 node tetrahedra is H4 a 27 node brick is
%                       B27 etc. Presently defined are L2, Q4, and Q9.
%               .nDOFe: Number of DOFs per element
%               .conn:  Array of element connectivity (size ne x nne)
%               .x:     Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)
%               .DOF:   Array of DOF indices (size nn x nsd)
%
%   Material:   Structure array with the following fields,
%               .E:     Modulus of elasticity
%               .nu:    Poisson's ratio
%               .Dtype: 2D approximation ('PlaneStrain' or 'PlainStress')
%
% 	Quad: 	Structure array with following fields,
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

% initialize D matrix file pointer
[~,DMatrix_functn] = fileparts(Material.ConstitutiveLawFile);

if nargin < 4
    calc_type = 'center';
end

% Specify dimension of the strain/stress matrix
switch Mesh.nsd
    case 1
        dim = 1;
    case 2
        dim = 3;
    case 3
        dim = 6;
end

if strcmp(calc_type, 'none')
    strain = zeros(dim, Mesh.ne);
    stress = zeros(dim, Mesh.ne);
else
    % Specify type of strain/stress matrix/cell
    switch calc_type
        case 'nodal'
            strain = zeros(dim, Mesh.nn);
            stress = zeros(dim, Mesh.nn);
            count = zeros(dim, Mesh.nn);
        case 'center'
            strain = zeros(dim, Mesh.ne);
            stress = zeros(dim, Mesh.ne);
        case 'L2projection'
            vec_size = Mesh.ne*Mesh.nne;
            A = zeros(vec_size,1); % matrix of integrals of shape functions - vectorized
            row = zeros(vec_size,1); % vector of row indices
            col = zeros(vec_size,1); % vector of column indices
            count = 1;
            
            be = zeros(Mesh.nn, dim); % column vectors of strains
            bs = zeros(Mesh.nn, dim); % column vectors of stresses
    end
    
    % Loop through all elements
    for e = 1:Mesh.ne
        
        %% Element variables
        % nodal ids of the element's nodes
        enodes = Mesh.conn(e,:);
        % global coordinates of the element's nodes
        xI = Mesh.x(enodes,:);
        % DOFs of element nodes
        dofE = Mesh.DOF(enodes,:);
        dofE = reshape(dofE',Mesh.nDOFe,[]);
        % local displacement vector
        de = d(dofE);
        % element material type
        nMat = Mesh.MatList(e); 
        
        switch calc_type
            case 'center'
                xi = zeros(1,Mesh.nsd);
                [~, dNdxi] = lagrange_basis(Mesh.type, xi, Mesh.nsd);
                Je = dNdxi'*xI;
                B = Je\(dNdxi');
                Bv = getBv(B', Mesh.nsd);
                strain(:, e) = Bv'*de;
                stress(:, e) = D*strain(:, e);
            case 'nodal'
                % initialize strain element
                strain_e = zeros(dim, Mesh.nne);
                stress_e = zeros(dim, Mesh.nne);
                
                % loop through all nodes and calculate nodal strains/stresses
                for n = 1:Mesh.nne
                    
                    % node point in parent coordinate
                    xi = getXI(n, Mesh.type);
                    
                    % Shape function and derivatives in parent coordinates
                    [~, dNdxi] = lagrange_basis(Mesh.type, xi, Mesh.nsd);
                    
                    % Jacobian of the transformation between parent and global
                    % coordinates
                    Je = dNdxi'*xI;
                    
                    % derivative of shape function in physical coordinates
                    % (tensor form)
                    dNdxi = dNdxi';
                    B = Je\dNdxi;
                    
                    % convert B matrix to Voigt form
                    Bv = getBv(B', Mesh.nsd);
                    
                    % strain_e = [strainx_n1  strainx_n2...;
                    %strainy_n1  strainy_n2...;
                    %strainxy_n1 strainxy_n2...];
                    strain_e(:,n) = Bv'*de;
                    stress_e(:,n) = D*strain_e(:,n);
                end
                % Add to global strains
                strain(:,enodes) = strain(:,enodes) + strain_e;
                stress(:,enodes) = stress(:,enodes) + stress_e;
                count(:,enodes) = count(:,enodes) + 1;
            case 'L2projection'
                % initialize elemental matrices and vectors
                A_e = zeros(Mesh.nne, Mesh.nne);
                de_e = zeros(Mesh.nne,dim);    % column vector of strains exx
                ds_e = zeros(Mesh.nne,dim);    % column vector of stresses sxx
                
                % strain in element from previous step
                strain_e_old = strain_old(:,e);
                
                % loop through all quadrature points
                for q = 1:Quad.nq
                    
                    % Shape functions and derivatives in parent coordinates
                    N = Quad.Nq{q}';
                    dNdxi = Quad.dNdxiq{q};
                    
                    % Jacobian of the transformation between parent and global
                    % coordinates
                    Je = dNdxi'*xI;
                    
                    % determinant of the Jacobian
                    dJe = det(Je);
                    
                    % derivative of shape function in physical coordinates
                    % (tensor form)
                    dNdxi = dNdxi';
                    B = Je\dNdxi;
                    
                    % convert B matrix to Voigt form
                    Bv = getBv(B', Mesh.nsd);
                    
                    % calculate strain at quadrature point
                    strain_q = Bv'*de;
                    
                    % calculate constitutive matrix
                    D = feval(DMatrix_functn, nMat, Material, Mesh, strain_e_old);
                    
                    % calculate stress at quadrature point
                    stress_q = D*strain_q;
                    
                    % Element level integral of L2-projections
                    A_e = A_e + N'*N*Quad.W(q)*dJe;
                    
                    for i = 1:dim
                        de_e(:,i) = de_e(:,i) +  N'*Quad.W(q)*dJe*strain_q(i);
                        ds_e(:,i) = ds_e(:,i) +  N'*Quad.W(q)*dJe*stress_q(i);
                    end
                end
                
                % Add to global matrices and vectors
                % Form the vectorized A matrix
                nAe = Mesh.nne^2; % number of entries in A matrix
                count = count + nAe;
                A_e = reshape(A_e, [nAe,1]);
                rowmatrix = enodes'*ones(1,Mesh.nne);
                row_e = reshape(rowmatrix, [nAe,1]);
                col_e = reshape(rowmatrix',[nAe,1]);
                
                A(count-nAe:count-1) = A_e;
                row(count-nAe:count-1) = row_e;
                col(count-nAe:count-1) = col_e;
                
                % Add element vectors to global vectors
                be(enodes,:) = be(enodes,:) + de_e;
                bs(enodes,:) = bs(enodes,:) + ds_e;
        end
    end
    
    
    switch calc_type
        case 'nodal'
            % For nodal strains, divide by count to get the average
            strain = strain./count;
            stress = stress./count;
        case 'L2projection'
            % For L2 projection, solve the system of equations and assemble the
            % nodal matrix of stresses and strains
            A = sparse(row, col, A, Mesh.nn, Mesh.nn);
            eL2 = zeros(size(be));
            sL2 = zeros(size(bs));
            
            for i = 1:dim
                eL2(:,i) = A\be(:,i);
                sL2(:,i) = A\bs(:,i);
            end
            
            strain = eL2';
            stress = sL2';
    end
end

end