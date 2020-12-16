function [strain, stress] = getStrain(d, Mesh, Material, calc_type, Quad)
%GETSTRAIN Evaluate stress and strain
%   strain = GETSTRAIN(d, Mesh, Material) is a cell array of  
%   nodal strains in each element of the mesh. The cell array is of size 
%   dim x ne, in which dim = 1 for 1D elements, 3 for 2D elements, and 6 
%   for 3D elements.
%
%   [strain, stress] = GETSTRAIN(d, Mesh, Material) also returns a cell 
%   array of nodal stresses in each element of the mesh (size dim x ne).
% 
%   [strain, stress] = GETSTRAIN(d, Mesh, Material, 'nodal') returns 
%   matrices of nodal-averaged strains (size dim x nn).
% 
%   [strain, stress] = GETSTRAIN(d, Mesh, Material, 'center') returns 
%   matrices of strains computed at the center of each element 
%   (size dim x ne).
%
%   [strain, stress] = GETSTRAIN(d, Mesh, Material, 'L2projection') returns 
%   matrices of L2-projected stresses at the nodes (size dim x nn).
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
    calc_type = 'none';
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

% Specify type of strain/stress matrix/cell
switch calc_type
    case 'none'
        strain = cell(dim,Mesh.ne);
        stress = cell(dim,Mesh.ne);
    case 'nodal'
        strain = zeros(dim, Mesh.nn);
        stress = zeros(dim, Mesh.nn);
        count = zeros(dim, Mesh.nn);
    case 'center'
        strain = zeros(dim, Mesh.ne);
        stress = zeros(dim, Mesh.ne);
    case 'L2projection'
        if Mesh.nsd ~= 2
            error('L2 projection of stresses currently only supported for 2 dimensions')
        end
        vec_size = Mesh.ne*Mesh.nne;
        A = zeros(vec_size,1); % matrix of integrals of shape functions - vectorized
        row = zeros(vec_size,1); % vector of row indices
        col = zeros(vec_size,1); % vector of column indices
        count = 1;
        
        dexx = zeros(Mesh.nn,1);    % column vector of strains exx
        deyy = zeros(Mesh.nn,1);    % column vector of strains eyy
        dexy = zeros(Mesh.nn,1);    % column vector of strains exy
        dsxx = zeros(Mesh.nn,1);    % column vector of stresses sxx
        dsyy = zeros(Mesh.nn,1);    % column vector of stresses syy
        dsxy = zeros(Mesh.nn,1);    % column vector of stresses sxy
end

% Loop through all elements
if ~strcmp(calc_type,'none')
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

        switch calc_type
            case 'center'
                xi = zeros(1,Mesh.nsd);
                [N, dNdxi] = lagrange_basis(Mesh.type, xi, Mesh.nsd);
                Xi = xI'*N;
                Je = dNdxi'*xI;
                B = Je\(dNdxi');
                Bv = getBv(B', Mesh.nsd);
                D = getD(Material.E(Xi), Material.nu(Xi), Mesh.nsd, Material.Dtype);
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
                    [N, dNdxi] = lagrange_basis(Mesh.type, xi, Mesh.nsd);

                    % quadrature point in physical coordinates
                    Xi = xI'*N;

                    % Jacobian of the transformation between parent and global 
                    % coordinates
                    Je = dNdxi'*xI;

                    % derivative of shape function in physical coordinates 
                    % (tensor form)
                    dNdxi = dNdxi';
                    B = Je\dNdxi;

                    % convert B matrix to Voigt form
                    Bv = getBv(B', Mesh.nsd);

                    D = getD(Material.E(Xi), Material.nu(Xi), Mesh.nsd, Material.Dtype);    

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
                dexx_e = zeros(Mesh.nne,1);    % column vector of strains exx
                deyy_e = zeros(Mesh.nne,1);    % column vector of strains eyy
                dexy_e = zeros(Mesh.nne,1);    % column vector of strains exy
                dsxx_e = zeros(Mesh.nne,1);    % column vector of stresses sxx
                dsyy_e = zeros(Mesh.nne,1);    % column vector of stresses syy
                dsxy_e = zeros(Mesh.nne,1);    % column vector of stresses sxy
                
                % loop through all quadrature points
                for q = 1:Quad.nq     

                    % Shape functions and derivatives in parent coordinates
                    N = Quad.Nq{q}';
                    dNdxi = Quad.dNdxiq{q};

                    % quadrature point in physical coordinates
                    Xi = xI'*N';

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

                    D = getD(Material.E(Xi), Material.nu(Xi), Mesh.nsd, Material.Dtype);    

                    % calculate stress and strain at quadrature point
                    strain_q = Bv'*de;
                    stress_q = D*strain_q;
                        
                    % Element level integral of L2-projections
                    A_e = A_e + N'*N*Quad.W(q)*dJe;
                    dexx_e = dexx_e + N'*Quad.W(q)*dJe*strain_q(1);
                    deyy_e = deyy_e + N'*Quad.W(q)*dJe*strain_q(2);
                    dexy_e = dexy_e + N'*Quad.W(q)*dJe*strain_q(3);
                    dsxx_e = dsxx_e + N'*Quad.W(q)*dJe*stress_q(1);
                    dsyy_e = dsyy_e + N'*Quad.W(q)*dJe*stress_q(2);
                    dsxy_e = dsxy_e + N'*Quad.W(q)*dJe*stress_q(3);
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
                dexx(enodes) = dexx(enodes) + dexx_e;
                deyy(enodes) = deyy(enodes) + deyy_e;
                dexy(enodes) = dexy(enodes) + dexy_e;
                dsxx(enodes) = dsxx(enodes) + dsxx_e;
                dsyy(enodes) = dsyy(enodes) + dsyy_e;
                dsxy(enodes) = dsxy(enodes) + dsxy_e;            
                
        end
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
        exxL2 = A\dexx;
        eyyL2 = A\deyy;
        exyL2 = A\dexy;
        
        strain = [exxL2';
                  eyyL2';
                  exyL2'];
              
        sxxL2 = A\dsxx;
        syyL2 = A\dsyy;
        sxyL2 = A\dsxy;
        
        stress = [sxxL2';
                  syyL2';
                  sxyL2'];
end

end