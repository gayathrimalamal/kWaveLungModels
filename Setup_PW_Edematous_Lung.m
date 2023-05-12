%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abnormal lung model  - Presence of intermittent B-lines
% Script for simulating an edematous lung model with k-Wave toolbox for PW
% transmissions. 
% Author: Gayathri Malamal (121814001@smail.iitpkd.ac.in;
%                           malamalgayathri@gmail.com)
% Date: 12-05-2023
% Cite as: G. Malamal and M. R. Panicker, "Pixel Intensity Vector Field:
% An Inside Out Approach of Looking at Ultrasound Reflections from the 
% Lung at High Frame Rates," 2021 43rd Annual International Conference of
% the IEEE Engineering in Medicine & Biology Society (EMBC), Mexico, 2021,
% pp. 2708-2711, doi: 10.1109/EMBC46164.2021.9629896
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;

addpath(genpath(pwd));
filepath = pwd;
dataset_name = 'Edema_Lung';

GPU             = false;   % default running on CPU
if (gpuDeviceCount)
    GPU         = true;   % set to true to simulate using GPU, to false to simulate in CPU
end

USE_CPU_C_CODE  = false;  % set to true to simulate using CPU cores in the absence of CUDA
SimPlot         = false; % set to true to view wave propagation. Will be active only if USE_CPU_CORES == false
RecordMovie     = false; % set to true to record the movie of simulation, else set to false
DO_BEAMFORMING  = true;  % set to true to perform beamforming with the generated dataset, else set to false

%% Define L11-5V transducer
Trans.name = 'L11-5v';                % Transducer name
Trans.tx_N_elements = 128;            % Number of transducer elements
Trans.tx_element_width = 270e-6;      % Width of a single transducer element [m]
Trans.tx_kerf = 30e-6;                % Kerf of the transducer [m]
Trans.tx_pitch = 300e-6;              % Pitch of the transducer [m]
Trans.tx_elevation_height = 5e-3;     % Elevation height of the transducer [m]
Trans.tx_elevation_focus = 18e-3;     % Elevation focus of the transducer [m]
Trans.tx_Fc = 4.8e6;                  % Transmit frequencies [Hz]
Trans.tx_Fmax = 1.2*Trans.tx_Fc;      % Maximum possible transducer center frequency [Hz]
Trans.tx_c0 = 1540;                   % SoS of the transducer [m/s]

Trans.tx_width = Trans.tx_N_elements * Trans.tx_pitch;      % Total width of the transducer [m]
Trans.angles = linspace(-18,18,7);     % Steering angle of the transmit signal between [-18, +18] [degrees]
Trans.na = length(Trans.angles);              % Number of plane wave steering angles
probe_geometry = linspace(-Trans.tx_width/2, Trans.tx_width/2, Trans.tx_N_elements)';

%% Define the k-Wave grid

%Here x is the depth direction, y is the sensor (or lateral) direction 

% set the size of the perfectly matched layer (PML)
PML_X_size = 10;
PML_Y_size = 10;

% Number of grid points choosen based on 2 points per wavelength and assumed transmit frequency of 3.8 MHz
c_min = 330;
lambda = c_min/Trans.tx_Fc;

% Calculate the spacing between the grid points
dx = lambda/2;                % [m]
dy = dx;                      % [m]

Nx = 1200 - 2 * PML_X_size;   % [grid points]
Ny = 1200 - 2 * PML_Y_size;   % [grid points]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% create the time array
t_end = (Nx * dx) * 2.2 / Trans.tx_c0;   % [s]
CFL = 0.05;
kgrid.t_array = makeTime(kgrid, Trans.tx_c0, CFL, t_end);

% create the time array
t_end = (Nx * dx) * 2.2 / Trans.tx_c0;   % [s]
CFL = 0.05;
kgrid.t_array = makeTime(kgrid, Trans.tx_c0, CFL, t_end);

%% Define the medium parameters
med = Edematous_Lung(Nx,Ny);

medium.sound_speed = med.sound_speed_map;
medium.density = med.density_map;
medium.alpha_coeff = med.alpha_coeff;
medium.alpha_power = 2;
medium.BonA = 6;

figure,imagesc(medium.sound_speed(:,:,1)); % Plotting the speed of sound map
title('Abnormal lung - SoS map');
colorbar;
figure,imagesc(medium.density (:,:,1)); % Plotting the density map
title('Abnormal lung - Density map');
colorbar;

%% Define the input signal

source_strength     = 4e6;      % [Pa], Maximum is 4e6
source_cycles       = 3;        % number of tone burst cycles
source_freq         = Trans.tx_Fc;
source.p_mask   = zeros(Nx, Ny);
source_width    = round(Trans.tx_element_width/dy);
source_kerf     = round(Trans.tx_kerf/dy);
element_height  = round(Trans.tx_elevation_height/dx);
rotation        = 0;

% Create empty kWaveArray
karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);

% Creating each transducer element by combining grid points
for Nc = 1: Trans.tx_N_elements
    karray.addRectElement([kgrid.x_vec(1), probe_geometry(Nc, 1)], dy, Trans.tx_element_width, rotation);
end

karray.setArrayPosition([0, 0], 0);
source.p_mask = karray.getArrayBinaryMask(kgrid);

display_mask=medium.sound_speed/norm(medium.sound_speed);
display_mask=histeq(display_mask);
display_mask(display_mask<0.7)=0;
display_mask(display_mask>0.7)=1;

sensor.mask = zeros(Nx, Ny);
sensor.mask = source.p_mask;

%% Run the simulation

acq = 'PW';
for txIdx = 1: length(Trans.angles)
    
    disp(txIdx);
    
    % Calculate the steering delay
    if(Trans.angles(txIdx)>0)
        Trans.delay(:,txIdx) = (probe_geometry(Trans.tx_N_elements,1)-probe_geometry(1:Trans.tx_N_elements,1))*sind(Trans.angles(txIdx))/(Trans.tx_c0*kgrid.dt);
    else
        Trans.delay(:,txIdx) = (probe_geometry(1,1)-probe_geometry(1:Trans.tx_N_elements,1))*sind(((Trans.angles(txIdx))))/(Trans.tx_c0*kgrid.dt);
    end
    
    % Note that the source strength is not scaled by acoustic impedance as a pressure source is used (Z=P/V)
    source_sig = source_strength*toneBurst(1/kgrid.dt, source_freq, source_cycles,'SignalOffset',Trans.delay(:,txIdx));
    source.p = karray.getDistributedSourceSignal(kgrid, source_sig);
    
    if(GPU)
        DATA_CAST       = 'gpuArray-single';
        input_args = {'PMLSize', [PML_X_size, PML_Y_size], 'PlotPML', false, ...
            'PMLInside', false, 'PlotScale', [-1, 1]*source_strength, ...
            'DisplayMask', display_mask, 'DataCast', DATA_CAST, 'DataRecast', true,...
            'RecordMovie', RecordMovie, 'MovieProfile', 'MPEG-4'};
        if(USE_CPU_C_CODE == true)
            sensor_data_grid_points = kspaceFirstOrder2DC(kgrid, medium, source, sensor, input_args{:});
        elseif(SimPlot == true)
            sensor_data_grid_points = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
        else
            sensor_data_grid_points = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});
        end
    else
        DATA_CAST       = 'single';
        if(SimPlot == false)
            input_args = {'PMLSize', [PML_X_size, PML_Y_size], 'PMLInside', false, 'DataCast', DATA_CAST,  'DataRecast', true,'PlotSim', false};
        else
            input_args = {'PMLSize', [PML_X_size, PML_Y_size], 'PlotPML', false, ...
            'PMLInside', false, 'PlotScale', [-1, 1]*source_strength, ...
            'DisplayMask', display_mask, 'DataCast', DATA_CAST, 'DataRecast', true,...
            'RecordMovie', RecordMovie, 'MovieProfile', 'MPEG-4'};
        end
        sensor_data_grid_points = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    end
    sensor_data(:,:,txIdx) = karray.combineSensorData(kgrid, sensor_data_grid_points); % Combine the grid points to get in terms of physical transducer elements
    
end

RcvData = permute(sensor_data, [2 1 3]);

% Saving the dataset
dataset = dataExtractFunc(dataset_name,acq,Trans,probe_geometry,RcvData,kgrid,medium);

%%
if(DO_BEAMFORMING)
    
    % DAS PW Beamforming
    Fnum = 1.5; % Receive F# for beamforming

    acq = 'PW';

    rawData=double(dataset.rawData);
    Fs = dataset.Fs*1e6;
    timeVector = (0:(size(rawData,1)-1)).'/(Fs);
    
    z_axis = 0.5*1540*timeVector;
    [x_grid,z_grid] = meshgrid(dataset.probe_geometry,z_axis);
    
    x_lim = [min(x_grid(:)) max(x_grid(:))]*1e3;
    z_lim = [min(z_grid(:)) max(z_grid(:))]*1e3;
    
    ang = dataset.angles;

    for txIdx = 1: length(ang)
        beamformedDataPW (:,:,txIdx) = DAS_RF(acq,squeeze(rawData(:, :, txIdx)), timeVector, x_grid, z_grid, dataset.probe_geometry, -ang(txIdx), dataset.Trans.tx_c0*ones(size(x_grid)), Fnum);
    end
    
    beamformedDataDAS = sum(beamformedDataPW(:,:,:), 3);
    envelopeDAS = abs(hilbert(beamformedDataDAS));
    beamformedDataDASImage = (envelopeDAS(:,:)./max(max(envelopeDAS(:,:))));
    figure,imagesc(dataset.probe_geometry.*1000,z_axis.*1000,20*log10(beamformedDataDASImage));
    colormap(gray);
    colorbar;
    vrange = [-60 0];
    caxis(vrange);
    xlabel('x [mm]');
    ylabel('z [mm]');
    title('Edema-Lung');
    set(gca,'fontsize',20);
    axis([x_lim z_lim]);
    
    fsavepath = strcat(filepath,filesep,'kWaveImages',filesep,dataset_name);
    if ~exist(fsavepath, 'dir')
        mkdir(fsavepath);
    end
    
    saveas(gcf,[fsavepath,filesep,dataset_name,'_',acq,'_',dataset.timeStamp,'.fig']);
    
end

