

% =========================================================================
% ESPIRiT Reconstruction 
% =========================================================================
% Weight the eigenvectors with soft-senses eigen-values
weights = W(:,:,end-num_maps+1:end);
weights = (weights - eigThresh_2)./(1-eigThresh_2).* (W(:,:,end-num_maps+1:end) > eigThresh_2);
weights = -cos(pi*weights)/2 + 1/2;
nIterCG = 100;

ESP = ESPIRiT(maps, weights);
WaveESP = WaveESPIRiT(maps, weights, [], PsfY);

tic;
if useVCC == 0
    [reskESPIRiT, resESPIRiT] = cgESPIRiT(DATA2, ESP, nIterCG, 0.0, DATA2*0);
    [reskWaveESPIRiT, resWaveESPIRiT] = cgESPIRiT(DATA_wave2, WaveESP, nIterCG, 0.0, DATA_wave2*0);
    % [reskWaveESPIRiT_v2, resWaveESPIRiT_v2] = cgESPIRiT_v2(DATA_wave2, WaveESP, 80, 0.0, DATA_wave2*0);
elseif useVCC == 1
    [reskVCCESPIRiT, resVCCESPIRiT] = cgESPIRiT(DATA2, ESP, nIterCG, 0.0, DATA2*0);
    [reskVCCWaveESPIRiT, resVCCWaveESPIRiT] = cgESPIRiT(DATA_wave2, WaveESP, nIterCG, 0.0, DATA_wave2*0);
end
toc;


if 1
A0 = abs(resVCCWaveESPIRiT(:,:,1));
% A = abs(ifft2c(fft2c(A0) .* (hamming(size(A0,1))*hamming(size(A0,2))')));
% A = abs(ifft2c(fft2c(A0) .* (tukeywin(size(A0,1),0.85)*tukeywin(size(A0,2),0.85)')));
A = abs(ifft2c(fft2c(A0) .* (hanning(size(A0,1))*hanning(size(A0,2))')));

% % -------------------------------------
% cd('cs_sense');
% XFM = Wavelet('Daubechies',4,2);
% TV = TVOP;
% cd('../');
% 
% cd('cs_sense');
% reg = abs(resVCCWaveESPIRiT(:,:,1));
% [fe,pe] = size(reg);
% img_reg = zpad(reg, [256,256]);
% img_xfm = XFM * img_reg;
% img_xfm = SoftThresh(img_xfm, 6e4);
% img_reg = XFM' * img_xfm;
% img_reg = crop(img_reg, [fe,pe]);
% reg = img_reg;
% cd('../');
% 
% A = reg;
% % -------------------------------------

A = (A - min(A(:)))./(max(A(:))-min(A(:)));
B = 1 - exp(-A*25);
weights(:,:,1) = B;

tic;
WaveESP = WaveESPIRiT(maps, weights, [], PsfY);
[reskVCCWaveESPIRiT_M, resVCCWaveESPIRiT_M] = cgESPIRiT(DATA_wave2, WaveESP, nIterCG, 0.0, DATA_wave2*0);
toc;
end


% QZL
A0 = abs(resVCCESPIRiT(:,:,1));
A = abs(ifft2c(fft2c(A0) .* (hanning(size(A0,1))*hanning(size(A0,2))')));

A = (A - min(A(:)))./(max(A(:))-min(A(:)));
B = 1 - exp(-A*25);
weights(:,:,1) = B;

ESP = ESPIRiT(maps, weights);
[reskVCCESPIRiT_M, resVCCESPIRiT_M] = cgESPIRiT(DATA2, ESP, nIterCG, 0.0, DATA2*0);
% END


if 0
% ESPIRiT CG reconstruction with 1 map
SNS = ESPIRiT(maps(:,:,:,end));
WaveSNS = WaveESPIRiT(maps(:,:,:,end), weights(:,:,end), [], PsfY);

if useVCC == 0
    [reskSENSE, resSENSE] = cgESPIRiT(DATA2, SNS, nIterCG, 0.0, DATA2*0);
    [reskWaveSENSE, resWaveSENSE] = cgESPIRiT(DATA_wave2, WaveSNS, nIterCG, 0.0, DATA_wave2*0);
elseif useVCC == 1
    [reskVCCSENSE, resVCCSENSE] = cgESPIRiT(DATA2, SNS, nIterCG, 0.0, DATA2*0);
    [reskVCCWaveSENSE, resVCCWaveSENSE] = cgESPIRiT(DATA_wave2, WaveSNS, nIterCG, 0.0, DATA_wave2*0);
end

% reg = abs(resVCCESPIRiT);
% reg = abs(ifft2c(fft2c(reg) .* repmat(hamming(size(reg,1))*hamming(size(reg,2))', [1,1,2])));
% [reskVCCESPIRiT_reg, resVCCESPIRiT_reg] = cgESPIRiT_QZL(DATA2, ESP, nIterCG, 0.0, reg, DATA2*0);
% 
% reg = abs(resVCCWaveESPIRiT);
% reg = abs(ifft2c(fft2c(reg) .* repmat(hamming(size(reg,1))*hamming(size(reg,2))', [1,1,2])));
% [reskVCCWaveESPIRiT_reg, resVCCWaveESPIRiT_reg] = cgESPIRiT_QZL(DATA_wave2, WaveESP, nIterCG, 0.0, reg, DATA_wave2*0);
end






