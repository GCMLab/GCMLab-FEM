function Bv = getBv(B, dim)
%GETBV Shape function derivative matrix in Voigt form
%   Bv = GETBV(B, dim) transforms the matrix of shape function 
%   derivatives, B, to Voigt notation, returning a matrix with n*dim rows 
%   in which n is the number of element nodes and dim is the number of 
%   spatial dimensions). For dim = 1, the resulting matrix has one column;
%   for  dim = 2, the resulting matrix has 3 columns; and for dim = 3, 
%   the matrix has 6 columns.
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   B:      vector of shape function derivatives (size n x dim where n 
%           is the number of element nodes)
%   dim:    number of spatial dimensions
% 
%   --------------------------------------------------------------------
%   Example
%   --------------------------------------------------------------------
%   B = |N1,1 N2,1 ... Nn,1|^T
%       |N1,2 N2,2 ... Nn,2|
%   dim = 2
%   Bv = |N1,1 0     N2,1 0     ...  Nn,1 0   |^T
%        |0    N1,2  0    N2,2  ...  0    Nn,2|
%        |N1,2 N1,1  N2,2 N2,1  ...  Nn,2 Nn,1|

%   Acknowledgements: Matin Parchei Esfahani, Jack Chessa

% Number of degrees of freedom
n = size(B,1);

switch dim
case 1
	Bv = B;
case 2
	Bv = zeros(dim*n,3);

    Bv(1:dim:dim*n-1,1) = B(:,1);
    Bv(2:dim:dim*n,2)   = B(:,2);

    Bv(1:dim:dim*n-1,3) = B(:,2);
    Bv(2:dim:dim*n,3)   = B(:,1);

case 3
	Bv = zeros(dim*n,6);

    Bv(1:dim:dim*n-2,1) = B(:,1);
    Bv(2:dim:dim*n-1,2) = B(:,2);
    Bv(3:dim:dim*n,3)   = B(:,3);

    Bv(2:dim:dim*n-1,4) = B(:,3);
    Bv(3:dim:dim*n,4)   = B(:,2);

    Bv(3:dim:dim*n,5)   = B(:,1);
    Bv(1:dim:dim*n-2,5) = B(:,3);

    Bv(1:dim:dim*n-2,6) = B(:,2);
    Bv(2:dim:dim*n-1,6) = B(:,1);

end
	

end