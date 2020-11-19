function K = getK(Mesh, Quad, Material)
% Acknowledgements: Chris Ladubec, Endrina Rivas

% initialize stiffness matrix
vec_size = Mesh.ne*(Mesh.nne * Mesh.nsd)^2; % vector size (solid dofs)
row = zeros(vec_size, 1);                   % vector of row indices
col = zeros(vec_size, 1);                   % vector of column indices
Kvec = zeros(vec_size, 1);                  % vectorized stiffness matrix

count = 1;                                  % DOF counter

% for each element, compute element stiffness matrix and add to global
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
        ndofE = numel(dofE);
        
    %% Shape functions and derivatives in parent coordinates
        W = Quad.W;
        nq = Quad.nq;

    %% Assemble stiffness matrix
    
        % length of element. used to check that quadrature points and weights 
        % are correct. A = Sum Wi*Ji
        A = 0;

        % initialize element stiffness matrix
        Ke = zeros(ndofE, ndofE); 

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

            % derivative of shape function in physical coordinates 
            % (tensor form)
            dNdxi = dNdxi';
            B = Je\dNdxi;

            % convert B matrix to Voigt form
            Bv = getBv(B', Mesh.nsd);

            D = getD(Xi, Mesh.nsd, Material);    
            
            % for 2D, volume integral includes the thickness
            switch Mesh.nsd 
                case 2
                    L = Material.t(Xi);
                case 3
                    L = 1;                
            end

            % Calculate local stiffness matrix
            Ke = Ke + W(q)*Bv*D*Bv'*L*dJe;

            % quadrature debug tool
            A = A + W(q)*dJe; 
        end
   
    %% Forming the vectorized stiffness matrix
        count = count + ndofE^2;
        Ke = reshape(Ke, [ndofE^2, 1]);
        rowmatrix = dofE*ones(1,ndofE);        
        rowe = reshape(rowmatrix,[ndofE^2, 1]);
        cole = reshape(rowmatrix',[ndofE^2, 1]);
        
        Kvec(count-ndofE^2:count-1) = Ke;
        row(count-ndofE^2:count-1) = rowe;
        col(count-ndofE^2:count-1) = cole;
end

% sparse stiffness matrix
K = sparse(row, col, Kvec, Mesh.nDOF, Mesh.nDOF);

end