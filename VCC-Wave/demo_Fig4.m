

% =========================================================================
% Figure 4: imaging of high resolution
% =========================================================================
% Add path
addpath('utils');
addpath('ESPIRiT_code');
addpath('coilCompression_code');
addpath('ismrm');

% Load data
load('data\HighRes\DATA.mat');
load('data\HighRes\DATA_wave.mat');
load('data\HighRes\PsfY.mat');

% SENSE and Wave
useVCC = 0;
Recon;

% Reload data
load('data\HighRes\DATA.mat');
load('data\HighRes\DATA_wave.mat');
load('data\HighRes\PsfY.mat');

% VCC and VCC-Wave
useVCC = 1;
Recon;

% Load full-sampled data as reference
load('data\HighRes\DATA.mat'); % reload fully sampled data
load('data\HighRes\mask_roi.mat'); % image ROI mask

% Calculate error maps and RMSEs
Process;

% Imshow
figure(4); 
imshow([img_it_sense,img_it_wave;resVCCESPIRiT,resVCCWaveESPIRiT],[]);

figure(5);
imshow([img_it_sense_err,img_it_wave_err;resVCCESPIRiT_err,resVCCWaveESPIRiT_err],[]); % error maps
colorbar; % and then choose jet






