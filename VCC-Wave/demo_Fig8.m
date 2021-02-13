

% =========================================================================
% Figure 8: brain imaging containing blood vessel signals
% =========================================================================
% Add path
addpath('utils');
addpath('ESPIRiT_code');
addpath('coilCompression_code');
addpath('ismrm');

% Load data
load('data\VesselSig\DATA.mat');
load('data\VesselSig\DATA_wave.mat');
load('data\VesselSig\PsfY.mat');

% SENSE and Wave
useVCC = 0;
Recon;

% Reload data
load('data\VesselSig\DATA.mat');
load('data\VesselSig\DATA_wave.mat');
load('data\VesselSig\PsfY.mat');

% VCC and VCC-Wave
useVCC = 1;
Recon;


% Modified ESPIRiT reconstruction for VCC-Wave
A = abs(resVCCWaveESPIRiT(:,:,1));
A = abs(ifft2c(fft2c(A) .* (hanning(size(A,1))*hanning(size(A,2))')));
A = (A - min(A(:)))./(max(A(:))-min(A(:)));
B = 1 - exp(-A*25);

weights_VW = weights;
weights_VW(:,:,1) = B;

tic;
WaveESP = WaveESPIRiT(maps, weights_VW, [], PsfY);
[reskVCCWaveESPIRiT_M, resVCCWaveESPIRiT_M] = cgESPIRiT(DATA_wave, WaveESP, nIterCG, 0.0, zeros(size(DATA_wave)));
toc;


% Modified ESPIRiT reconstruction for VCC
A = abs(resVCCESPIRiT(:,:,1));
A = abs(ifft2c(fft2c(A) .* (hanning(size(A,1))*hanning(size(A,2))')));
A = (A - min(A(:)))./(max(A(:))-min(A(:)));
B = 1 - exp(-A*25);

weights_W = weights;
weights_W(:,:,1) = B;

tic;
ESP = ESPIRiT(maps, weights_W);
[reskVCCESPIRiT_M, resVCCESPIRiT_M] = cgESPIRiT(DATA, ESP, nIterCG, 0.0, zeros(size(DATA)));
toc;


% Load full-sampled data as reference
load('data\VesselSig\DATA.mat'); % reload fully sampled data
load('data\VesselSig\mask_roi.mat'); % image ROI mask

% Calculate error maps and RMSEs
Process_Fig8;

% Imshow
figure(8); 
imshow([img_it_sense,img_it_wave;resVCCESPIRiT,resVCCWaveESPIRiT;resVCCESPIRiT_M,resVCCWaveESPIRiT_M],[0,3.0]);

figure(9);
imshow([img_it_sense_err,img_it_wave_err;resVCCESPIRiT_err,resVCCWaveESPIRiT_err;resVCCESPIRiT_M_err,resVCCWaveESPIRiT_M_err],[]);
colorbar; % and then choose jet

% see RMSEs
img_it_sense_rmse, img_it_wave_rmse,
resVCCESPIRiT_rmse, resVCCWaveESPIRiT_rmse,
resVCCESPIRiT_M_rmse, resVCCWaveESPIRiT_M_rmse










