function K = getK(Mesh, Quad, Material)
%GETK 
%   K = getK()
%
%   ----------------------------------------------------------------------
%   Created by Endrina Rivas
%       endrina.rivas@uwaterloo.ca
%       Department of Civil Engineering
%       University of Waterloo
%   
%   Last updated: June 2016
%   ----------------------------------------------------------------------

% initialize stiffness matrix
K = sparse(Mesh.nDOF, Mesh.nDOF); 

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
            B = dNdxi/Je;

            % convert B matrix to Voigt form
            Bv = getBv(B, Mesh.nsd);

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
    
    %% Assemble element matrices
        K(dofE,dofE) = K(dofE,dofE) + Ke;
end

end