%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defines a specular reflector in a homogeneous medium
% Author: Gayathri Malamal(121814001@smail.iitpkd.ac.in;
%                          malamalgayathri@gmail.com)
% Date: 04-05-2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function medParam = Consolidation_Lung(Nx,Ny)
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

%Fluid+Air
cAir = 600;          % Speed of sound [m/s]
rhoAir = 400;        % Density [kg/m^3]
alphaAir= 2;         % Attenuation [db/cm/MHz]

%Pleura
cPleura = cSoftTissue;          % Speed of sound [m/s]
rhoPleura = rhoSoftTissue;      % Density [kg/m^3]
alphaPleura = alphaSoftTissue;  % Attenuation [db/cm/MHz]

%Scarred Pleura
cScarredPleura =  1200;         % Speed of sound [m/s]
rhoScarredPleura = 800;         % Density [kg/m^3]
alphaScarredPleura  = 2;        % Attenuation [db/cm/MHz] 

cFluid      = 900;          % Speed of sound [m/s]
rhoFluid    = 600;        % Density [kg/m^3]
alphaFluid  = 2;      % Attenuation [db/cm/MHz] 
%=========================================================================
%Creating background maps of random scatterers
%=========================================================================
rng(1); % Seed initialization to ensure same scatterer distribution in every dataset
scat_dist = randn(Nx,Ny); % Random distribution of scatterers

background_map_mean = 1;
background_map_std = 0.05;
background_map = background_map_mean + background_map_std*scat_dist;

background_map_std_skin = 0.005;
background_map_skin = background_map_mean + background_map_std_skin*scat_dist;

background_map_mean_air = 2;
background_map_std_air = 0.09;
background_map_air = background_map_mean_air + background_map_std_air*scat_dist;

% =========================================================================
% Defining consolidated lung
% =========================================================================

medParam.sound_speed_map = cSoftTissue* ones(Nx, Ny) .* background_map;
medParam.density_map = rhoSoftTissue* ones(Nx, Ny) .* background_map;
medParam.alpha_coeff = alphaSoftTissue* ones(Nx, Ny);

%Define the skin layer
skin_layer = makeLine(Nx, Ny,  [1 1], [1 Ny]);
for i = 1:30
    skin_layer = skin_layer+makeLine(Nx, Ny,  [1+i 1], [1+i Ny]);
end

slab                             = zeros(Nx, Ny);
slab(348:348+40, :)              = 1; % Thickened Pleura
slab(348+41: 348+50, 1:5:end)    = 3; % Scarred Pleura
slab(348+51: 348+60, 1:8:end)    = 4; % Air
slab(348+41:348+60, [160:316, 660:810])      = 1; %Pleura
slab(348+61:end,:)                           = 2; %Fluid+Air(Below pleura)

consolidation = zeros(Nx, Ny);
consolidation = consolidation +  makeDisc(Nx, Ny, 348+90-5, 200, 40);
consolidation = consolidation +  makeDisc(Nx, Ny, 360+90-5, 278, 40);
consolidation = consolidation +  makeDisc(Nx, Ny, 360+155-5, 200, 40);
consolidation = consolidation +  makeDisc(Nx, Ny, 360+139-5, 340, 40);

consolidation = consolidation +  makeDisc(Nx, Ny, 348+90-5, 200+500, 40);
consolidation = consolidation +  makeDisc(Nx, Ny, 360+100-5, 278+500, 40);
consolidation = consolidation +  makeDisc(Nx, Ny, 360+165-5, 200+500, 40);
consolidation = consolidation +  makeDisc(Nx, Ny, 360+139-5, 340+500, 40);

medParam.sound_speed_map(skin_layer==1) = cSkin.*background_map_skin(skin_layer == 1);
medParam.density_map(skin_layer==1) = rhoSkin.*background_map_skin(skin_layer == 1);
medParam.alpha_coeff(skin_layer == 1) = alphaSkin;

medParam.sound_speed_map(slab == 1) = cPleura;
medParam.sound_speed_map(slab == 2) = cAir.*background_map_air(slab == 2);
medParam.sound_speed_map(slab == 3) = cScarredPleura.*background_map(slab == 3);
medParam.sound_speed_map(slab == 4) = cFluid.*background_map(slab == 4);
medParam.sound_speed_map(consolidation == 1) = cPleura;

medParam.density_map(slab == 1) = rhoPleura;
medParam.density_map(slab == 2) = rhoAir.*background_map_air(slab == 2);
medParam.density_map(slab == 3) = rhoScarredPleura.*background_map(slab == 3);
medParam.density_map(slab == 4) = rhoFluid.*background_map(slab == 4);
medParam.density_map(consolidation == 1) = rhoPleura;

medParam.alpha_coeff(slab == 1) = alphaPleura;
medParam.alpha_coeff(slab == 2) = alphaAir;
medParam.alpha_coeff(slab == 3) = alphaScarredPleura;
medParam.alpha_coeff(slab == 4) = alphaFluid;
medParam.alpha_coeff(consolidation == 1)  = alphaPleura;

end

