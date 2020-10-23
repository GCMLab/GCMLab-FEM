function Mesh = NodeDOFs(Mesh)

    Mesh.nDOFe = Mesh.nne*Mesh.nsd;         % number of DOF per element
    Mesh.nDOF = Mesh.nn*Mesh.nsd;           % total number of DOF
    Mesh.DOF = zeros(Mesh.nn,Mesh.nsd); 

    for sd = 1:Mesh.nsd
       Mesh.DOF(:,sd) = (sd : Mesh.nsd : (Mesh.nDOF-(Mesh.nsd-sd)))';
    end

end