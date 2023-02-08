function Nv = getNv(N, dim)
%GETNV Shape function matrix in Voigt form
% 	Nv = GETNV(N, dim) transforms the vector of shape functions, N,
% 	into Voigt notation, returning a matrix of size n x dim where n is 
% 	the number of element nodes and dim is the number of spatial 
% 	dimensions. 
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
% 	N:		vector of shape functions (size n x 1 where n is the number 
% 			of element nodes)
% 	dim:	number of spatial dimensions
% 
%   --------------------------------------------------------------------
% 	Example
%   --------------------------------------------------------------------
% 	N = |N1 N2 ... Nn|^T
%	dim = 2
%	Nv = |N1 0   N2 0   ...  Nn 0  |^T
%      	 |0  N1  0  N2  ...  0  Nn |

% 	Acknowledgements: Matin Parchei Esfahani

n = size(N,1);

I = eye(dim);
Nv = [];
for i = 1:n
    Nv = [Nv;I*N(i)];
end

end