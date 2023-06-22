function D = getD_TH1(nMat, Material, ~)
%GETD_TH1 - Themal conductivity tensor
%   D = GETD_TH1(nMat, Material, Mesh) is the conductivity tensor.
%       Supports model options ISO - Isotropic or ORTHO - Orthotropic
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   nMat:       Material type
%   Mesh:       Structure array with the following fields, may be updated
%               with new fields
%               .nsd:   Number of spatial dimensions
%   Material:   Structure array with the following fields, may be updated
%               with new fields
%               .k
%   --------------------------------------------------------------------

k1 = Material.Prop(nMat).k1;
if isfield(Material.Prop(nMat),'k2')
    k2 = Material.Prop(nMat).k2;
end

switch Material.Dtype
    case 'ISO'
        D = k1;
    case 'ORTHO'
        D = [k1,    0;
             0,     k2];
end

end