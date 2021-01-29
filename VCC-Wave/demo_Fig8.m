




addpath('utils');
addpath('ESPIRiT_code');
addpath('coilCompression_code');
addpath('ismrm');

% =========================================================================
% Figure 4: imaging of high bandwidth
% =========================================================================
load('data\VesselSig\DATA.mat');
load('data\VesselSig\DATA_wave.mat');
load('data\VesselSig\PsfY.mat');

% SENSE and Wave
useVCC = 0;
Recon;

% VCC and VCC-Wave
useVCC = 1;
Recon;


% modified ESPIRiT reconstruction for VCC-Wave
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


% modified ESPIRiT reconstruction for VCC
A = abs(resVCCESPIRiT(:,:,1));
A = abs(ifft2c(fft2c(A) .* (hanning(size(A,1))*hanning(size(A,2))')));
A = (A - min(A(:)))./(max(A(:))-min(A(:)));
B = 1 - exp(-A*25);

weights_W = weights;
weights_W(:,:,1) = B;

tic;
ESP = ESPIRiT(maps, weights_W);
[reskVCCESPIRiT_M, resVCCESPIRiT_M] = cgESPIRiT(DATA2, ESP, nIterCG, 0.0, DATA2*0);
toc;


% calculate error maps and RMSEs
load('data\VesselSig\DATA.mat'); % reload fully sampled data
load('data\VesselSig\mask_roi.mat'); % image ROI mask
Process_Fig8;






