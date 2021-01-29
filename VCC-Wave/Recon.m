

% =========================================================================
% Subsampling	
% =========================================================================
[Nx,Ny,Nc] = size(DATA);
Nx_os = size(PsfY,1);

acs_len = 48; % ACS lines
CalibSize = [acs_len, acs_len];
kCalib = crop(DATA,[CalibSize,Nc]);

Ry = 6; % acceleration factor

mask = zeros(Nx,Ny);
mask(:,floor(Ny/2)+1:-Ry:1) = 1;
mask(:,floor(Ny/2)+1:Ry:end) = 1;

mask_wave = zeros(Nx_os,Ny);
mask_wave(:,floor(Ny/2)+1:-Ry:1) = 1;
mask_wave(:,floor(Ny/2)+1:Ry:end) = 1;

DATA = DATA .* repmat(mask,[1,1,Nc]);
DATA_wave = DATA_wave .* repmat(mask_wave,[1,1,Nc]); 


% =========================================================================
% Coil Compression
% =========================================================================
% nCHA_cc = 16; % compressed coil number
% [sccmtx] = calcSCCMtx(kCalib);
% sccmtx_cc = sccmtx(:,1:nCHA_cc);
% 
% kCalib = CC(kCalib,sccmtx_cc);
% DATA = CC(DATA,sccmtx_cc);
% DATA_wave = CC(DATA_wave,sccmtx_cc);
% Nc = nCHA_cc;


% =========================================================================
% Extend along coil dimension
% =========================================================================
PsfY = repmat(PsfY,[1,1,Nc]);
samp_mat = repmat(mask,[1,1,Nc]);
samp_mat_wave = repmat(mask_wave,[1,1,Nc]);


% =========================================================================
% Virtual Conjugate Coil	
% =========================================================================
if useVCC == 1
    v_kCalib = conj(flip(flip(kCalib,1),2)); 
    if mod(size(v_kCalib,1),2) == 0
        v_kCalib = circshift(v_kCalib,[1,0]);
    end
    if mod(size(v_kCalib,2),2) == 0
        v_kCalib = circshift(v_kCalib,[0,1]);
    end
    kCalib = cat(3, kCalib, v_kCalib);
    
    v_DATA = conj(flip(flip(DATA,1),2)); 
    if mod(size(v_DATA,1),2) == 0
        v_DATA = circshift(v_DATA,[1,0]);
    end
    if mod(size(v_DATA,2),2) == 0
        v_DATA = circshift(v_DATA,[0,1]);
    end
    DATA = cat(3,DATA,v_DATA);
    
    v_DATA_wave = conj(flip(flip(DATA_wave,1),2)); 
    if mod(size(v_DATA_wave,1),2) == 0
        v_DATA_wave = circshift(v_DATA_wave,[1,0]);
    end
    if mod(size(v_DATA_wave,2),2) == 0
        v_DATA_wave = circshift(v_DATA_wave,[0,1]);
    end
    DATA_wave = cat(3,DATA_wave,v_DATA_wave);
    
    v_samp_mat = conj(flip(flip(samp_mat,1),2)); 
    if mod(size(v_samp_mat,1),2) == 0
        v_samp_mat = circshift(v_samp_mat,[1,0]);
    end
    if mod(size(v_samp_mat,2),2) == 0
        v_samp_mat = circshift(v_samp_mat,[0,1]);
    end
    samp_mat = cat(3,samp_mat,v_samp_mat);
    
    v_samp_mat_wave = conj(flip(flip(samp_mat_wave,1),2)); 
    if mod(size(v_samp_mat_wave,1),2) == 0
        v_samp_mat_wave = circshift(v_samp_mat_wave,[1,0]);
    end
    if mod(size(v_samp_mat_wave,2),2) == 0
        v_samp_mat_wave = circshift(v_samp_mat_wave,[0,1]);
    end
    samp_mat_wave = cat(3,samp_mat_wave,v_samp_mat_wave);

    PsfY = cat(3,PsfY,conj(circshift(flip(PsfY,1),[1,0]))); % please note
    Nc = Nc * 2; % coils doubled
end


% =========================================================================
% Sensitivity Esimatioon
% =========================================================================
kSize = [8,8];

if useVCC == 0
    eigThresh_1 = 0.02;
    eigThresh_2 = 0.95;
else
    eigThresh_1 = 0.015;
    eigThresh_2 = 0.75;
end

tic;
[K,S] = dat2Kernel(kCalib, kSize);
idx = find(S >= S(1)*eigThresh_1, 1, 'last');
[M,W] = kernelEig(K(:,:,:,1:idx), [Nx,Ny]);
toc;


% =========================================================================
% Reconstruction 
% =========================================================================
if useVCC == 0 % SENSE recon
    ecalib_map = M(:,:,:,end) .* repmat(W(:,:,end) > eigThresh_2, [1,1,Nc]);
    
    [img_it_sense] = ismrm_cartesian_iterative_SENSE_QZL(DATA, samp_mat, ecalib_map);
    [img_it_wave] = ismrm_cartesian_iterative_Wave(DATA_wave, samp_mat_wave, ecalib_map, PsfY);
elseif useVCC == 1 % ESPIRiT recon with 2 sets of maps
    num_maps = 2;
    maps = M(:,:,:,end-num_maps+1:end);

    weights = W(:,:,end-num_maps+1:end);
    weights = (weights - eigThresh_2)./(1-eigThresh_2).* (W(:,:,end-num_maps+1:end) > eigThresh_2);
    weights = -cos(pi*weights)/2 + 1/2;
    nIterCG = 100;
    
    ESP = ESPIRiT(maps, weights);
    WaveESP = WaveESPIRiT(maps, weights, [], PsfY);
    
    [reskVCCESPIRiT, resVCCESPIRiT] = cgESPIRiT(DATA, ESP, nIterCG, 0.0, zeros(size(DATA)));
    [reskVCCWaveESPIRiT, resVCCWaveESPIRiT] = cgESPIRiT(DATA_wave, WaveESP, nIterCG, 0.0, zeros(size(DATA_wave)));
end


















