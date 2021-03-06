


img_ref = sos(ifft2c(DATA)) .* mask_roi;
img_it_sense = abs(img_it_sense) .* mask_roi;
img_it_wave = abs(img_it_wave) .* mask_roi;
resVCCESPIRiT = sos(resVCCESPIRiT) .* mask_roi;
resVCCWaveESPIRiT = sos(resVCCWaveESPIRiT) .* mask_roi;


% intensity normalization
IMG = img_ref;            N3_v2;  img_ref = img;
IMG = img_it_sense;       N3_v2;  img_it_sense = img;
IMG = img_it_wave;        N3_v2;  img_it_wave = img;
IMG = resVCCESPIRiT;      N3_v2;  resVCCESPIRiT = img;
IMG = resVCCWaveESPIRiT;  N3_v2;  resVCCWaveESPIRiT = img;


img_ref = img_ref ./ mean(img_ref(mask_roi == 1));
img_it_sense = img_it_sense ./ mean(img_it_sense(mask_roi == 1));
img_it_wave = img_it_wave ./ mean(img_it_wave(mask_roi == 1));
resVCCESPIRiT = resVCCESPIRiT ./ mean(resVCCESPIRiT(mask_roi == 1));
resVCCWaveESPIRiT = resVCCWaveESPIRiT ./ mean(resVCCWaveESPIRiT(mask_roi == 1));


img_it_sense_err = abs(img_it_sense - img_ref);
img_it_wave_err = abs(img_it_wave - img_ref);
resVCCESPIRiT_err = abs(resVCCESPIRiT - img_ref);
resVCCWaveESPIRiT_err = abs(resVCCWaveESPIRiT - img_ref);


img_it_sense_rmse = sqrt(sum(sum(abs(img_it_sense_err).^2,1),2)) / sqrt(sum(sum(abs(img_ref).^2,1),2));
img_it_wave_rmse = sqrt(sum(sum(abs(img_it_wave_err).^2,1),2)) / sqrt(sum(sum(abs(img_ref).^2,1),2));
resVCCESPIRiT_rmse = sqrt(sum(sum(abs(resVCCESPIRiT_err).^2,1),2)) / sqrt(sum(sum(abs(img_ref).^2,1),2));
resVCCWaveESPIRiT_rmse = sqrt(sum(sum(abs(resVCCWaveESPIRiT_err).^2,1),2)) / sqrt(sum(sum(abs(img_ref).^2,1),2));






