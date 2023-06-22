function [Material, stiffnessmatrixfile_name, stressstrainfile_name] = setMaterialModel(Material,Control)
%SETMATERIALMODEL - creates pointers to material model files

    switch Material.Model
        case 'LE1' % Linear elasticity
            Material.ConstitutiveLawFile = 'getD';
            Material.StiffnessMatrixFile = 'getK_LE1';
            Material.StressStrainFile = 'getStrain';
            Control.transient = 0; 
        case 'LET1' % Linear Elastic Transient
            Material.ConstitutiveLawFile = 'getD';
            Material.StiffnessMatrixFile = 'getK_LET1';
            Material.StressStrainFile = 'getStrain';
            Control.transient = 1; 
        case 'ST1' % Stiffening Model via 1st Invariant of strain
            Material.ConstitutiveLawFile = 'getD_ST1';
            Material.StiffnessMatrixFile = 'getK_ST1';
            Material.StressStrainFile = 'getStrain_ST1';
            Control.transient = 0; 
        case 'TR2' % Transient model with stiffening model via 1st invariant of strain
            Material.ConstitutiveLawFile = 'getD_ST1';
            Material.StiffnessMatrixFile = 'getK_TR2';
            Material.StressStrainFile = 'getStrain_ST1';
            Control.transient = 1; 
        case 'VE1' % Linear Viscoelastic Kelvin-Voigt Model
            Material.ConstitutiveLawFile = 'getD';
            Material.StiffnessMatrixFile = 'getK_VE1';
            Material.StressStrainFile = 'getStrain';
            Material.DampingFile = 'getC_VE1'; 
            Control.transient = 1; 
    end

    % Create function pointers
        % Tangent Matrix
        [~,stiffnessmatrixfile_name] = fileparts(Material.StiffnessMatrixFile);

        % Stress/strain 
        [~,stressstrainfile_name] = fileparts(Material.StressStrainFile);
        
    % D matrix is called inside assembly functions
end

