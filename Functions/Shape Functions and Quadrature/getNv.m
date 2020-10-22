function Nv = getNv(N,dim)
%GETNv

% Acknowledgements: Matin Parchei Esfahani
% N = |N1 N2 ... Nn|^T
%
% Nv = |N1 0   N2 0   ...  Nn 0  |^T
%      |0  N1  0  N2  ...  0  Nn |

n = size(N,1);

I = eye(dim);
Nv = [];
for i = 1:n
    Nv = [Nv;I*N(i)];
end

end