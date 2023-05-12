%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --Script for saving the RF data and beamforming parameters
% --Authors: Gayathri M and Mahesh Raveendranatha Panicker
% --Inputs:
%         dataset_name - Name for the dataset
%         acq - Type of transmit (PW or STA)
%         Trans - Trans structure that contain the transmit parameters
%         probe_geometry - Geometry of the probe (x) 
%         RcvData - Sensor Data generated in kWave (time samples x channels)
%         kgrid - kgrid structure
%         medium - kWave medium structure
% --Date: 04-05-2023
% --Center for Computational Imaging, Indian Institute of Technology Palakkad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  dataset = dataExtractFunc(dataset_name,acq,Trans,probe_geometry,RcvData,kgrid,medium)

% =========================================================================
% Creating Save Directory
% =========================================================================
filepath = pwd;
date_Str=datestr(now, 'dd_mmm_yyyy_HH_MM_SS_AM');
fsaveDir = strcat(filepath,filesep,'kWaveData',filesep,dataset_name);

if ~exist(fsaveDir, 'dir')
    mkdir(fsaveDir);
end

% =========================================================================
% Downsampling the k-Wave data
% =========================================================================
% The data in k_wave is heavily oversampled.
% To reduce the sampling rate to 5 time the transmit frequency.
Fs = (1/kgrid.dt); %Sampling frequency in MHz
freqDown = round(Fs/(5*Trans.tx_Fc));

% =========================================================================
% Extracting raw RF data: Synthetic Transmit Aperture
% =========================================================================
if(strcmp(acq, 'STA'))
    active_aperture = Trans.tx_active_SA ;
    for txIdx = 1:length(active_aperture)
        rawData_resample(:,:,txIdx) = resample(RcvData(:,:,txIdx),1,freqDown);
    end
end

% =========================================================================
% Extracting raw RF data: Plane Wave Transmission
% =========================================================================
if(strcmp(acq, 'PW'))
    angles = Trans.angles;
    dataset.angles = angles.*pi/180; %angles in radian
    
    for txIdx = 1:length(angles)
        init_delay = (Trans.delay(:, txIdx).*kgrid.dt);
        maxDelay = (max(init_delay)/2)*(1/kgrid.dt);
        
        if(maxDelay==0)
            maxDelay=1;
        end
        
        rawDataTemp = squeeze(RcvData(round(maxDelay):end,:,txIdx));
        rawData(1:size(rawDataTemp,1),:,txIdx) = rawDataTemp;
    end
    for txIdx = 1:length(angles)
        rawData_resample(:,:,txIdx) = resample(rawData(:,:,txIdx),1,freqDown);
    end
end

% =========================================================================
% Saving variables to dataset structure
% =========================================================================
dataset.probe_geometry = probe_geometry;
dataset.Trans = Trans; %Tra
dataset.c0 = Trans.tx_c0; %Speed of sound
dataset.Fc = Trans.tx_Fc/1e6; %Transmit center frequency in MHz
dataset.Fs = (Fs*1e-6/freqDown); %Sampling frequency in MHz
dataset.medium = medium; %Medium Paramaters
dataset.timeStamp = date_Str; %Time stamp
dataset.rawData = single(rawData_resample);

save([fsaveDir,filesep,'dataset','_',dataset_name, '_',acq,'_',date_Str,'.mat'], 'dataset','-v7.3');
