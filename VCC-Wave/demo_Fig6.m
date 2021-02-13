

% =========================================================================
% Figure 6: imaging of high bandwidth
% =========================================================================
% Add path
addpath('utils');
addpath('ESPIRiT_code');
addpath('coilCompression_code');
addpath('ismrm');

% Load data
load('data\HighBW\DATA.mat');
load('data\HighBW\DATA_wave.mat');
load('data\HighBW\PsfY.mat');

% SENSE and Wave
useVCC = 0;
Recon;

% Reload data
load('data\HighBW\DATA.mat');
load('data\HighBW\DATA_wave.mat');
load('data\HighBW\PsfY.mat');

% VCC and VCC-Wave
useVCC = 1;
Recon;

% Load full-sampled data as reference
load('data\HighBW\DATA.mat'); % reload fully sampled data
load('data\HighBW\mask_roi.mat'); % image ROI mask

% Calculate error maps and RMSEs
Process;

% Imshow
figure(6); 
imshow([img_it_sense,img_it_wave;resVCCESPIRiT,resVCCWaveESPIRiT],[]);

figure(7);
imshow([img_it_sense_err,img_it_wave_err;resVCCESPIRiT_err,resVCCWaveESPIRiT_err],[]); % error maps
colorbar; % and then choose jet

% see RMSEs
img_it_sense_rmse, img_it_wave_rmse,
resVCCESPIRiT_rmse, resVCCWaveESPIRiT_rmse






