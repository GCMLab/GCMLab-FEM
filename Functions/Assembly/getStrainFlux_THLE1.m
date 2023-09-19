function [strain, stress, gradT, flux] = getStrainFlux_THLE1(d, Mesh, Material, calc_type, Quad,~)
%GETSTRAINFLUX_THLE1 Evaluate stress, strain, temperature gradient, and
%   flux
%   [strain, stress, gradT, flux] = GETSTRAINFLUX_THLE1(d, Mesh, Material)
%   returns four matrices of gradients and fluxes/stresses computed at the
%   center of each element. The matrices are of size dim x ne, in which
%   dim = 1 for 1D elements, 3 for 2D elements, and 6 for 3D elements.
%
%   [strain, stress, gradT, flux] = GETSTRAINFLUX_THLE1(d, Mesh, Material, 'none') does not
%   compute the stresses/fluxes or gradients
%
%   [strain, stress, gradT, flux] = GETSTRAINFLUX_THLE1(d, Mesh, Material, 'nodal') returns
%   matrices of nodal-averaged gradients (size dim x nn).
%
%   [strain, stress, gradT, flux] = GETSTRAINFLUX_THLE1(d, Mesh, Material, 'center') returns
%   matrices of gradients computed at the center of each element
%   (size dim x ne).
%
%   [strain, stress, gradT, flux] = GETSTRAINFLUX_THLE1(d, Mesh, Material, 'L2projection') returns
%   matrices of L2-projected stresses/fluxes at the nodes (size dim x nn).
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

if nargin < 4
    calc_type = 'center';
end

% Specify dimension of the strain/stress matrix
switch Mesh.nsd
    case 1
        dim_mech = 1;
    case 2
        dim_mech = 3;
    case 3
        dim_mech = 6;
end

% Specify dimension of the gradT/flux matrix
dim_sca = Mesh.nsd;

if strcmp(calc_type, 'none')
    strain = zeros(dim_mech, Mesh.ne);
    stress = zeros(dim_mech, Mesh.ne);
    gradT = zeros(dim_sca, Mesh.ne);
    flux = zeros(dim_sca, Mesh.ne);
else
    % Specify type of strain/stress matrix/cell
    switch calc_type
        case 'nodal'
            strain = zeros(dim_mech, Mesh.nn);
            stress = zeros(dim_mech, Mesh.nn);
            gradT = zeros(dim_sca, Mesh.nn);
            flux = zeros(dim_sca, Mesh.nn);
            count_mech = zeros(dim_mech, Mesh.nn);
            count_sca = zeros(dim_sca, Mesh.nn);
        case 'center'
            strain = zeros(dim_mech, Mesh.ne);
            stress = zeros(dim_mech, Mesh.ne);
            gradT = zeros(dim_sca, Mesh.ne);
            flux = zeros(dim_sca, Mesh.ne);
        case 'L2projection'
            vec_size = Mesh.ne*Mesh.nne;
            A_mech = zeros(vec_size,1); % matrix of integrals of shape functions - vectorized
            A_sca = zeros(vec_size,1); % matrix of integrals of shape functions - vectorized
            
            row = zeros(vec_size,1); % vector of row indices
            col = zeros(vec_size,1); % vector of column indices
            count = 1;
            
            b_strain = zeros(Mesh.nn, dim_mech); % column vectors of strains
            b_stress = zeros(Mesh.nn, dim_mech); % column vectors of stresses
            b_gradT = zeros(Mesh.nn, dim_sca); % column vectors of strains
            b_flux = zeros(Mesh.nn, dim_sca); % column vectors of stresses
    end
    
    % Loop through all elements
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
        % local displacement vector
        de_mech = d(dofE_mech);
        de_sca = d(dofE_sca);
        %% Constitutive matrices
        nMat = Mesh.MatList(e); % element material type
        % mechanic problem
        Material.Dtype = Material.Dtype_mech;
        D_mech = getD(nMat, Material, Mesh);
        % thermal problem
        Material.Dtype = Material.Dtype_therm;
        D_sca = getD_TH1(nMat, Material, Mesh);
        
        switch calc_type
            case 'center'
                xi = zeros(1,Mesh.nsd);
                [~, dNdxi] = lagrange_basis(Mesh.type, xi, Mesh.nsd);
                Je = dNdxi'*xI;
                B = Je\(dNdxi');
                Bv = getBv(B', Mesh.nsd);
                strain(:, e) = Bv'*de_mech;
                stress(:, e) = D_mech*strain(:, e);
                gradT(:, e) = B*de_sca;
                flux(:, e) = D_sca*gradT;
            case 'nodal'
                % initialize strain element
                strain_e = zeros(dim_mech, Mesh.nne);
                stress_e = zeros(dim_mech, Mesh.nne);
                
                % initialize flux element
                gradT_e = zeros(dim_sca, Mesh.nne);
                flux_e = zeros(dim_sca, Mesh.nne);
                
                % loop through all nodes
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
                    
                    strain_e(:,n) = Bv'*de_mech;
                    stress_e(:,n) = D_mech*strain_e(:,n);
                    gradT_e(:,n) = B*de_sca;
                    flux_e(:,n) = D_sca*gradT_e;
                end
                % Add to global gradients
                strain(:,enodes) = strain(:,enodes) + strain_e;
                stress(:,enodes) = stress(:,enodes) + stress_e;
                count_mech(:,enodes) = count_mech(:,enodes) + 1;
                
                flux(:,enodes) = flux(:,enodes) + flux_e;
                gradT(:,enodes) = gradT(:,enodes) + gradT_e;
                count_sca(:,enodes) = count_sca(:,enodes) + 1;
            case 'L2projection'
                % initialize elemental matrices and vectors
                A_mech_e = zeros(Mesh.nne, Mesh.nne);
                strain_e = zeros(Mesh.nne,dim_mech);    % column vector of strains exx
                stress_e = zeros(Mesh.nne,dim_mech);    % column vector of stresses sxx
                
                A_sca_e = zeros(Mesh.nne, Mesh.nne);
                gradT_e = zeros(Mesh.nne,dim_sca);    % column vector of gradients
                flux_e = zeros(Mesh.nne,dim_sca);    % column vector of fluxes
                
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
                    
                    % calculate stress/flux and gradients at quadrature point
                    strain_q = Bv'*de_mech;
                    stress_q = D_mech*strain_q;
                    
                    gradT_q = B*de_sca;
                    flux_q = D_sca*gradT_q;
                    
                    % Element level integral of L2-projections
                    A_mech_e = A_mech_e + N'*N*Quad.W(q)*dJe;
                    A_sca_e = A_sca_e + N'*N*Quad.W(q)*dJe;
                    
                    for i = 1:dim_mech
                        strain_e(:,i) = strain_e(:,i) +  N'*Quad.W(q)*dJe*strain_q(i);
                        stress_e(:,i) = stress_e(:,i) +  N'*Quad.W(q)*dJe*stress_q(i);
                    end
                    
                    for i = 1:dim_sca
                        gradT_e(:,i) = gradT_e(:,i) +  N'*Quad.W(q)*dJe*gradT_q(i);
                        flux_e(:,i) = flux_e(:,i) +  N'*Quad.W(q)*dJe*flux_q(i);
                    end
                end
                
                % Add to global matrices and vectors
                % Form the vectorized A matrix
                nA_e = Mesh.nne^2; % number of entries in A matrix
                count = count + nA_e;
                A_mech_e = reshape(A_mech_e, [nA_e,1]);
                A_sca_e = reshape(A_sca_e, [nA_e,1]);
                
                rowmatrix = enodes'*ones(1,Mesh.nne);
                row_e = reshape(rowmatrix, [nA_e,1]);
                col_e = reshape(rowmatrix',[nA_e,1]);
                
                A_mech(count-nA_e:count-1) = A_mech_e;
                A_sca(count-nA_e:count-1) = A_sca_e;
                
                row(count-nA_e:count-1) = row_e;
                col(count-nA_e:count-1) = col_e;
                
                % Add element vectors to global vectors
                b_strain(enodes,:) = b_strain(enodes,:) + strain_e;
                b_stress(enodes,:) = b_stress(enodes,:) + stress_e;
                
                b_gradT(enodes,:) = b_gradT(enodes,:) + gradT_e;
                b_flux(enodes,:) = b_flux(enodes,:) + flux_e;                
        end
    end
    
    switch calc_type
        case 'nodal'
            % For nodal gradients, divide by count to get the average
            strain = strain./count_mech;
            stress = stress./count_mech;
            gradT = gradT./count_sca;
            flux = flux./count_sca;
        case 'L2projection'
            % For L2 projection, solve the system of equations and assemble the
            % nodal matrix of stresses/fluxes and gradients
            A_mech = sparse(row, col, A_mech, Mesh.nn, Mesh.nn);
            strainL2 = zeros(size(b_strain));
            stressL2 = zeros(size(b_stress));
            
            A_sca = sparse(row, col, A_sca, Mesh.nn, Mesh.nn);
            gradTL2 = zeros(size(b_gradT));
            fluxL2 = zeros(size(b_flux));
            
            for i = 1:dim_mech
                strainL2(:,i) = A_mech\b_strain(:,i);
                stressL2(:,i) = A_mech\b_stress(:,i);
            end
            
            for i = 1:dim_sca
                gradTL2(:,i) = A_sca\b_gradT(:,i);
                fluxL2(:,i) = A_sca\b_flux(:,i);
            end
            
            strain = strainL2';
            stress = stressL2';
            
            gradT = gradTL2';
            flux = fluxL2';
    end
    
end

end