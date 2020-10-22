function Material = Material_shale()
%% Material Properties - Shale

	% Young's modulus [Pa]
	Material.E = @(x) 7.5e9;  

	% Constitutive law: 'PlaneStrain' or 'PlaneStress' 
	Material.Dtype = 'PlaneStrain'; 

	% Thickness [m] (set as default to 1)
	Material.t = @(x) 1;

	% Poisson's ratio (set as default to 0)
	Material.nu = @(x) 0.25;

end