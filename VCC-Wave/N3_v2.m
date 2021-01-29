%截至2014-04-03，2D已经基本完成，resample/shrink=4效果可以
%下一步，开始3D的试验
% clc; clear;
% imtool close all;
% M, N: Size of source image IMG
% m, n: Size of down sampling image img
% v   : Log(img)
% V   : Probability density of v, Hist of v
% fX  : The FFT result of X

% Gaussian Distribution
% f(x) = 1/(σ(2π)^0.5) * exp(-(x-u)^2/2σ^2)
% σ: 标准差, u:期望值(Mu in English)
% σ=1时，称为标准正态分布
% FWHM: 2*sqrt(2ln2)*σ = 2.35482*σ

% [IMG,info] = loaddicom;
IMG = double(IMG);

[M, N] = size(IMG);

src_max = max(IMG(:));
level = graythresh(IMG/src_max);
BW = im2bw(IMG/src_max,level*0.0);
BW = imfill(BW,'holes');
% BW = ones(M,N);
%降采样4倍不错，越小校正越强
%作者使用的好像是3
resample = 4;
% pixelspacing = info.PixelSpacing;
% RoPeMax = double(max(M,N))*pixelspacing(1);
% resample = floor(1000/RoPeMax);

img = IMG(1:resample:M, 1:resample:N);
v_uc = log(img+eps);
[m,n] = size(v_uc);

bw = BW(1:resample:M, 1:resample:N);
BW = reshape(BW, 1, M*N);
bw = reshape(bw, 1, m*n);

nbins = 200;
FWHM  = 0.15;
Z     = 2;
Z2    = Z.^2;

padded_size = floor(2^(ceil(log(nbins)/log(2))+1)+0.5);
offset = (padded_size-nbins)/2;
F  = zeros(1,padded_size);
V_padded = zeros(1,padded_size);

fs = zeros(m,n);
fe = zeros(m,n);

tic;
err = 1;
loop = 1;
maxLoop = 500;
while loop<maxLoop && err>0.001
    v = v_uc-fs;
    v(v<0) = 0;
    v = reshape(v, 1, m*n);
    
%     V = hist(v(:), nbins);
%     vmax = max(v(:));
%     vmin = min(v(:));
%     slope  = (vmax-vmin)/(nbins-1);
    
    vtemp = v(find(bw));
    vmax = max(vtemp);
    vmin = min(vtemp);
    slope  = (vmax-vmin)/(nbins-1);
    V = hist(vtemp, nbins);
    
    %gaussian
    fwhm = FWHM/slope;
    g_factor = 4.0*log(2.0)/(fwhm*fwhm);
    g_scale  = 2.0*(log(2.0)/pi)^0.5/fwhm;

    F(1) = g_scale;
    x = 1:(padded_size-1)/2;
    F(x+1) = g_scale*exp(-(x).^2*g_factor);
    F(padded_size-x+1) = F(x+1);
    F(padded_size/2+1) = g_scale*exp(-padded_size.^2*g_factor/4);
    blur = fft(F);
    
%     x = -(padded_size/2-1):(padded_size/2);
%     F = g_scale*exp(-x.^2*g_factor);
%     blur = fft(fftshift(F));
    %end gaussian
    
    %weiner, Z2 is noise
    G = conj(blur)./(conj(blur).*blur + Z2);
    %end weiner
    
    V_padded(offset+1:offset+nbins) = V;
    fV_padded = fft(V_padded);
    U = real(ifft(fV_padded.*G));
    U(U<0) = 0;
    
    moment = U .* (vmin + ((0:padded_size-1)-offset)*slope);
    
    numerator   = real(ifft(fft(moment).*blur));
    denominator = real(ifft(fft(U).*blur));
    Y_padded = zeros(1,padded_size);
    Y_padded(denominator~=0) = numerator(denominator~=0) ./ denominator(denominator~=0);
    
    Y = Y_padded(offset+1:offset+nbins);
    u = zeros(size(vtemp));
    for ii=1:length(vtemp)
        rIdx = (vtemp(ii)-vmin)/slope + 1;
        iIdx = floor(rIdx);
        if iIdx < nbins
            u(ii) = Y(iIdx) + (Y(iIdx+1)-Y(iIdx))*(rIdx-iIdx);
        else
            u(ii) = Y(nbins);
        end
    end

    fe(find(bw)) = v_uc(find(bw)) - u;  % v or v_uc? ITK中使用v, 效果不好
    th = -max(fe(:));
    if th<0
        fe(fe<th) = min(th,0);
    end

    fs_last = fs;
    fs = smoothn(fe,5);
    
    r = exp(fs_last-fs);
    if mean(r(:))~=0
        err = abs(std(r(:))/mean(r(:)));
    end
    loop = loop+1;
end

fs = imresize(fs, [M, N]);
fs = exp(fs);
vfs = fs(find(BW));
minfs = min(vfs);
fs = fs ./ minfs;   %这个处理需要再确认

fs  = reshape(fs, 1, M*N);
img = reshape(IMG, 1, M*N);
for ii=1:M*N
    if BW(ii)
        img(ii) = IMG(ii)./fs(ii);
    elseif fs(ii)>1
        img(ii) = 0.5*IMG(ii)./fs(ii);  %这一步需要再确认
    end
end

fs  = reshape(fs, M, N);
img = reshape(img, M, N);
BW  = reshape(BW, M, N);
toc;


