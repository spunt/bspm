function [denoisedMap] = waveidtb_waveletDenoise(img , mask , ParameterInfo)

% Updated August 2011 (Wavelet filter kernels hardcoded for those who do
% not have the wavelet toolbox provided by MATLAB.

[x y z w] = size( img );
%% Deterending,reshaping and masking;
mixedMat = detrend(reshape(img,[],w),0);
clear img;

% Reshape in to an image
spatialData = reshape(mixedMat,[x y z w]);
clear mixedMat;

%% Get the wavelet filters (From MATLAB's WAVELET TOOLBOX)

if strcmp(ParameterInfo.waveType,'sym2')
    
    Lo_D = [-0.129409522550921,0.224143868041857,0.836516303737469,0.482962913144690];
    Hi_D = [-0.482962913144690,0.836516303737469,-0.224143868041857,-0.129409522550921];
    Lo_R = [0.482962913144690,0.836516303737469,0.224143868041857,-0.129409522550921];
    Hi_R = [-0.129409522550921,-0.224143868041857,0.836516303737469,-0.482962913144690];

elseif strcmp(ParameterInfo.waveType,'db2')
    
    Lo_D =  [-0.129409522550921,0.224143868041857,0.836516303737469,0.482962913144690];
    Hi_D =  [-0.482962913144690,0.836516303737469,-0.224143868041857,-0.129409522550921];
    Lo_R = [0.482962913144690,0.836516303737469,0.224143868041857,-0.129409522550921];
    Hi_R = [-0.129409522550921,-0.224143868041857,0.836516303737469,-0.482962913144690];

elseif strcmp(ParameterInfo.waveType,'coif2')
    
    Lo_D = [-0.000720549445364512,-0.00182320887070299,0.00561143481939450,0.0236801719463341,-0.0594344186464569,-0.0764885990783064,0.417005184421693,0.812723635445542,0.386110066821162,-0.0673725547219630,-0.0414649367817592,0.0163873364635221];
    Hi_D = [-0.0163873364635221,-0.0414649367817592,0.0673725547219630,0.386110066821162,-0.812723635445542,0.417005184421693,0.0764885990783064,-0.0594344186464569,-0.0236801719463341,0.00561143481939450,0.00182320887070299,-0.000720549445364512];
    Lo_R = [0.0163873364635221,-0.0414649367817592,-0.0673725547219630,0.386110066821162,0.812723635445542,0.417005184421693,-0.0764885990783064,-0.0594344186464569,0.0236801719463341,0.00561143481939450,-0.00182320887070299,-0.000720549445364512];
    Hi_R = [-0.000720549445364512,0.00182320887070299,0.00561143481939450,-0.0236801719463341,-0.0594344186464569,0.0764885990783064,0.417005184421693,-0.812723635445542,0.386110066821162,0.0673725547219630,-0.0414649367817592,-0.0163873364635221];
    
elseif strcmp(ParameterInfo.waveType,'bior1.3')
    
    Lo_D = [-0.0883883476483185,0.0883883476483185,0.707106781186548,0.707106781186548,0.0883883476483185,-0.0883883476483185];
    Hi_D = [0,0,-0.707106781186548,0.707106781186548,0,0];
    Lo_R = [0,0,0.707106781186548,0.707106781186548,0,0];
    Hi_R = [-0.0883883476483185,-0.0883883476483185,0.707106781186548,-0.707106781186548,0.0883883476483185,0.0883883476483185];
end 

% DOES NOT WORK IF WAVELETS TOOLBOX IS NOT INSTALLED (SEE CODE ABOVE)
% [Lo_D,Hi_D] = wfilters(ParameterInfo.waveType);
% [Lo_R,Hi_R] = wfilters(ParameterInfo.waveType,'r');

lev = ParameterInfo.levels + 1;
L = ParameterInfo.levels;

%% 3-D denoising module
inRecon = cell(1,1);
dat = spatialData.*mask;
[ lll , llh , lhl , lhh , hll , hlh , hhl , hhh ] = waveidtb_dwt3D( dat, lev, [Lo_D'./sqrt(2),Hi_D'./sqrt(2)]);

[sighhl_hhh, siglhl_lhh, sighll_hlh , siglll_llh ] = waveidtb_getStd( lhl , lhh , hll , hlh , hhl , hhh , lll , llh );
%     [sighhl_hhh, siglhl_lhh, sighll_hlh , siglll_llh ] = getStd( hhh , hhl , llh );
K = 2;

y_est_hhh = hhh;
y_est_hhl = hhl;
y_est_lhh = lhh;
y_est_lhl = lhl;
y_est_hlh = hlh;
y_est_hll = hll;
y_est_llh = llh;

for j = size(lll,2)-1:-1:1; % Loop through levels (starting from 2 in this case)
    
    [y_est_hhh{j}] = waveidtb_getEstimates3D( hhh{j} , y_est_hhh{j+1} , sighhl_hhh(j) , K);fprintf('.');
    [y_est_hhl{j}] = waveidtb_getEstimates3D( hhl{j} , y_est_hhl{j+1} , sighhl_hhh(j) , K);fprintf('.');
    
    [y_est_lhh{j}] = waveidtb_getEstimates3D( lhh{j} , y_est_lhh{j+1} , siglhl_lhh(j) , K);fprintf('.');
    [y_est_lhl{j}] = waveidtb_getEstimates3D( lhl{j} , y_est_lhl{j+1} , siglhl_lhh(j) , K);fprintf('.');
    
    [y_est_hlh{j}] = waveidtb_getEstimates3D( hlh{j} , y_est_hlh{j+1} , sighll_hlh(j) , K);fprintf('.');
    [y_est_hll{j}] = waveidtb_getEstimates3D( hll{j} , y_est_hll{j+1} , sighll_hlh(j) , K);fprintf('.');
    
    %                     [y_est_lll_xy{j}(:,:,z)] = lll{j}(:,:,z);
    [y_est_llh{j}] = waveidtb_getEstimates3D( llh{j} , y_est_hll{j+1} , siglll_llh(j) , K);fprintf('.');
    
    % Prepare data for reconstruction
    inRecon{j}(1) = { y_est_llh{j} };
    inRecon{j}(2) = { y_est_lhl{j} };
    inRecon{j}(3) = { y_est_lhh{j} };
    inRecon{j}(4) = { y_est_hll{j} };
    inRecon{j}(5) = { y_est_hlh{j} };
    inRecon{j}(6) = { y_est_hhl{j} };
    inRecon{j}(7) = { y_est_hhh{j} };
    inRecon{j}(8) = {lll(j)};
end;

% 3-D Wavelet Reconstruction
reconData = waveidtb_idwt3D(inRecon, L, [Lo_R'./sqrt(2),Hi_R'./sqrt(2)]);
%         output3D(i,:,:) = reshape(recLev3 , x*y, size(recLev3,3));
denoisedMap = circshift(abs(reconData)*1.8, [lev-1 lev-1 0 ]).*(mask>0); % circshift to compensate for sampling during filtering.
fprintf('Done.\n');
