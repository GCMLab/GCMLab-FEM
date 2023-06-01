function [Material, stiffnessmatrixfile_name, stressstrainfile_name] = setMaterialModel(Material)
%SETMATERIALMODEL - creates pointers to material model files

    switch Material.Model
        case 'LE1' % Linear elasticity
            Material.ConstitutiveLawFile = 'getD';
            Material.StiffnessMatrixFile = 'getK_LE1';
            Material.StressStrainFile = 'getStrain';
            Material.ProblemType      = 1;
        case 'ST1' % Stiffening Model via 1st Invariant of strain
            Material.ConstitutiveLawFile = 'getD_ST1';
            Material.StiffnessMatrixFile = 'getK_ST1';
            Material.StressStrainFile = 'getStrain_ST1';
            Material.ProblemType      = 1;
        case 'TH1'
            Material.ConstitutiveLawFile = 'getD_TH1';
            Material.StiffnessMatrixFile = 'getK_TH1';
            Material.StressStrainFile = 'getFlux_TH1';
            Material.ProblemType      = 2;
    end

    % Create function pointers
        % Tangent Matrix
        [~,stiffnessmatrixfile_name] = fileparts(Material.StiffnessMatrixFile);

        % Stress/strain 
        [~,stressstrainfile_name] = fileparts(Material.StressStrainFile);
        
    % D matrix is called inside assembly functions
end

