function Nv = getNv(N,dim)
%GETNv

%   ----------------------------------------------------------------------
%   Created by Endrina Rivas
%       endrina.rivas@uwaterloo.ca
%       Department of Civil Engineering
%       University of Waterloo
%       November 2015
%
%   Acknowledgements: Matin Parchei Esfahani
%   ----------------------------------------------------------------------
% N = |N1 N2 ... Nn|^T
%
% Nv = |N1 0   N2 0   ...  Nn 0  |^T
%      |0  N1  0  N2  ...  0  Nn |

n = size(N,1);

% Alternate:
% if isempty(N)
%     Nv = [];
%     return
% end
% Ntemp = reshape([N, zeros(n,dim-1)]',n*dim,1);
% for j = 1:dim
%     Nv(:,j) = [zeros(j-1,1);Ntemp(1:end-j+1)];
% end

I = eye(dim);
Nv = [];
for i = 1:n
    Nv = [Nv;I*N(i)];
end

end