%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Implements the DAS RF beamforming with pixel wise addition
%-- This is inspired from DAS in USTB                                      
%-- UltraSound ToolBox (https://www.ustb.no/)  
%-- Inputs:
%         acq - Type of transmit (PW or STA)
%         rawData - RF data for a single transmit
%         timeVector - Total time
%         x_grid - 'x' co-ordinates for every pixel
%         z_grid - 'z' co-ordinates for every pixel
%         probe_geometry - Geometry of the probe (x)
%         idx - PW angle or Transmit Center in STA
%         c_map - SoS map
%         Fnum - F-number for DAS beamforming
%-- Authors: Gayathri Malamal and Mahesh Raveendranatha Panicker
%-- Date: 04-05-2023
%-- Center for Computational Imaging, Indian Institute of Technology Palakkad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beamformedData = DAS_RF(acq,rawData,timeVector,x_grid,z_grid,probe_geometry,idx,c_map, Fnum)

N_c = size(rawData,2);
rows = size(x_grid,1);
columns = size(x_grid,2);
beamformedData = zeros(rows,columns);
delay_compensation = zeros(rows*columns,N_c);

%==========================================================================
% Calculate the transmit distance
% =========================================================================
if(strcmp(acq, 'PW'))
    transmit_distance = z_grid(:)*cos(idx)+x_grid(:)*sin(idx);
elseif(strcmp(acq, 'STA'))
    transmit_distance = sqrt((x_grid(:)-idx).^2+z_grid(:).^2); 
end

%==========================================================================
% Depth based apodization
%==========================================================================

rx_f_number = Fnum;
rx_aperture = z_grid(:)/rx_f_number;
aperture = rx_aperture*ones(1,N_c);
rx_aperture_distance = abs(x_grid(:)*ones(1,N_c)-ones(rows*columns,1)*probe_geometry(:,1).'); 

% Hanning Apodization
receive_apodization = double(rx_aperture_distance<=aperture/2).*(0.5 + 0.5*cos(2*pi*rx_aperture_distance./aperture)); 

%==========================================================================
% Delay compensation and beamforming
%==========================================================================

for nrx=1:N_c
    
    receive_distance = sqrt((probe_geometry(nrx,1)-x_grid(:)).^2+z_grid(:).^2);
    delay = (transmit_distance+receive_distance)./c_map(:);
       
    delay_compensation(:,nrx) = interp1(timeVector,rawData(:,nrx),delay,'spline',0);
end

beamformedData(:) = sum(receive_apodization.*delay_compensation,2);
beamformedData(isnan(beamformedData)) = 0;
beamformedData = reshape(beamformedData,size(x_grid));