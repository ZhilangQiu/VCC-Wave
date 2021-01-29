


addpath('utils');
addpath('ESPIRiT_code');
addpath('coilCompression_code');
addpath('ismrm');

% =========================================================================
% Figure 4: imaging of high bandwidth
% =========================================================================
load('data\HighRes\DATA.mat');
load('data\HighRes\DATA_wave.mat');
load('data\HighRes\PsfY.mat');

% SENSE and Wave
useVCC = 0;
Recon;

% VCC and VCC-Wave
useVCC = 1;
Recon;

% calculate error maps and RMSEs
load('data\HighRes\DATA.mat'); % reload fully sampled data
load('data\HighRes\mask_roi.mat'); % image ROI mask
Process;
















