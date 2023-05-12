%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defines a specular reflector in a homogeneous medium
% Author: Gayathri Malamal(121814001@smail.iitpkd.ac.in;
%                          malamalgayathri@gmail.com)
% Date: 04-05-2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function medParam = Edematous_Lung(Nx,Ny)
%==========================================================================
% Define the medium properties
% =========================================================================
%Skin layer
cSkin = 1700;     % Speed of sound [m/s]
rhoSkin = 1150;   % Density [kg/m^3]
alphaSkin = 0.6;  % Attenuation [db/cm/MHz]

% Soft tissue
cSoftTissue = 1540;         % Speed of sound [m/s]
rhoSoftTissue = 1100;       % Density [kg/m^3]
alphaSoftTissue = 0.5;      % Attenuation [db/cm/MHz]

%Air
cAir = 330;     % Speed of sound [m/s]
rhoAir = 100;   % Density [kg/m^3]
alphaAir= 14;   % Attenuation [db/cm/MHz]

%Pleura
cPleura = cSoftTissue;          % Speed of sound [m/s]
rhoPleura = rhoSoftTissue;      % Density [kg/m^3]
alphaPleura = alphaSoftTissue;  % Attenuation [db/cm/MHz]

%Fluid (Water at 20 degree)
cFluid =  1430;         % Speed of sound [m/s]
rhoFluid = 1000;        % Density [kg/m^3]
alphaFluid  = 0.0022;   % Attenuation [db/cm/MHz] 

%=========================================================================
%Creating background maps of random scatterers
%=========================================================================
rng(1); % Seed initialization to ensure same scatterer distribution in every dataset
scat_dist = randn(Nx,Ny); % Random distribution of scatterers

background_map_mean = 1;
background_map_std = 0.05;
background_map = background_map_mean + background_map_std *scat_dist;

background_map_std_skin = 0.005;
background_map_skin = background_map_mean + background_map_std_skin*scat_dist;

% =========================================================================
% Defining edematous lung
% =========================================================================

medParam.sound_speed_map = cSoftTissue* ones(Nx, Ny) .* background_map;
medParam.density_map = rhoSoftTissue * ones(Nx, Ny) .* background_map;
medParam.alpha_coeff = alphaSoftTissue* ones(Nx, Ny);

%Define the skin layer
skin_layer=makeLine(Nx, Ny,  [1 1], [1 Ny]);
for i=1:30
    skin_layer=skin_layer+makeLine(Nx, Ny,  [1+i 1], [1+i Ny]);
end

slab                    = zeros(Nx, Ny);
slab(348:348+10, :)        = 1; %Pleura
slab(348+11:end, :)        = 2; %Air

rad_cavity = 40;
fluid_region = makeDisc(Nx, Ny, 427+abs(rad_cavity-40)-30, 200, rad_cavity);

for i = 2*rad_cavity:83:Nx-800
    fluid_region =  fluid_region + makeDisc(Nx, Ny, 427+abs(rad_cavity-40)+i-30, 200, rad_cavity);
end

for k = 0:83:Nx-800
    for i = 600:400: Ny
        fluid_region = fluid_region +  makeDisc(Nx, Ny, 427+abs(rad_cavity-40)+k-30, i, rad_cavity);
    end
end

medParam.sound_speed_map(skin_layer==1) = cSkin.*background_map_skin(skin_layer == 1);
medParam.density_map(skin_layer==1) = rhoSkin.*background_map_skin(skin_layer == 1);
medParam.alpha_coeff(skin_layer == 1) = alphaSkin;

medParam.sound_speed_map(slab == 1) = cPleura;
medParam.sound_speed_map(slab == 2) = cAir;

medParam.density_map(slab == 1) = rhoPleura;
medParam.density_map(slab == 2) = rhoAir;

medParam.alpha_coeff(slab == 1) = alphaPleura;
medParam.alpha_coeff(slab == 2) = alphaAir;

medParam.sound_speed_map(fluid_region == 1) = cFluid;
medParam.density_map(slab == 1) = rhoPleura;
medParam.density_map(slab == 2) = rhoAir;
medParam.density_map(fluid_region == 1) = rhoFluid;

medParam.alpha_coeff(slab == 1) = alphaPleura;
medParam.alpha_coeff(slab == 2) = alphaAir;
medParam.alpha_coeff(fluid_region == 1) = alphaFluid;

end
