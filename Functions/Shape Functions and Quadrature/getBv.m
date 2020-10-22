function Bv = getBv(B, dim)
%GETB

%   Acknowledgements: Matin Parchei Esfahani, Jack Chessa

% B = |N1,1 N2,1 ... Nn,1|^T
%     |N1,2 N2,2 ... Nn,2|
%
% Bv = |N1,1 0     N2,1 0     ...  Nn,1 0   |^T
%      |0    N1,2  0    N2,2  ...  0    Nn,2|
%      |N1,2 N1,1  N2,2 N2,1  ...  Nn,2 Nn,1|

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