function bad_volumes = bramila_artrep(cfg)
% bramila_artrep.m
% This function computes EPI+motion data diagnostics as implemented in ArtRepair Toolbox v2.5.
% !! Data is not actually modified/fixed (run the toolbox if needed) !!
% input:
%   cfg = one-subject cfg structure
% output:
%   bad_volumes = list of bad volumes, empty for clean data, -1 if failed to run
%   corresponding diagnostics figure is saved into "bramila" sub-folder (may not work in command-line mode)

% default quality thresholds
MvmtThresh = 0.5;   % max allowed motion [mm]
PercentThresh=1.5;  % max signal change tolerance in percentages
ZThresh = 3.0;      % max signal change tolerance in Z-scores

fprintf('Computing ArtRepair diagnostics\n');

% load data
if isfield(cfg,'vol') && ~isempty(cfg.vol)
    img=double(cfg.vol);
elseif isfield(cfg,'infile') && ~isempty(cfg.infile)
    nii=load_nii(cfg.infile);
    img=double(nii.img);
else
    error('No input EPI data found!');
end

% if we have a mask and if we have 4D data, then apply the mask
hasmask=0;
if (isfield(cfg,'analysis_mask') && ~isempty(cfg.analysis_mask) && size(img,4)>0)
    mask=double(cfg.analysis_mask);
    sz=size(img);
    if(~any(size(mask) ==sz(1:3)))
        error(['The specified mask has a different size than the fMRI data. Quitting.'])
    end
    for t=1:size(img,4)
        img(:,:,:,t)=mask.*img(:,:,:,t);
    end
    hasmask=1;
elseif (isfield(cfg,'mask') && ~isempty(cfg.mask) && size(img,4)>0)
    mask=double(cfg.mask);
    sz=size(img);
    if(~any(size(mask) ==sz(1:3)))
        error(['The specified mask has a different size than the fMRI data. Quitting.'])
    end
    for t=1:size(img,4)
        img(:,:,:,t)=mask.*img(:,:,:,t);
    end
    hasmask=1;
end

if hasmask==1
    artrep_mask=mask;
else
    artrep_mask = [];
end

motion_params = load(cfg.motionparam);

prepro_suite='fsl-fs';
if ~isfield(cfg,'prepro_suite')
    cfg.prepro_suite = prepro_suite;
end

try
    bad_volumes = art_global_bramila(img,motion_params,artrep_mask,2,PercentThresh,ZThresh,MvmtThresh,cfg);
    fprintf('Total %i bad volumes found (%4.1f%%)\n',length(bad_volumes),100*length(bad_volumes)/size(img,4));
catch
    bad_volumes = -1;
    warning('ArtRepair diagnostics failed to run')
end

end


function bad_volumes = art_global_bramila(Images,RealignmentFile,Mask,RepairType,PercentThresh,ZThresh,MvThresh,cfg)
% MODIFIED art_global.m file (no need for SPM) 
% 3.7.2014 JanneK
%
% FORMAT art_global                         (v.2.5)
%
%     Art_global allows visual inspection of average intensity and
% scan to scan motion of fMRI data, and offers methods to repair outliers
% in the data. Outliers are scans whose global intensity is very different
% from the mean, or whose scan-to-scan motion is large. Thresholds that
% define the outliers are shown, and they can be adjusted by the user.
% Margins are defined around the outliers to assure that the repaired
% data will satisfy the "slowly-varying" background assumption of GLM
% models. Outliers can also be defined by user manual edits. When all
% repair parameters are set, the user writes new repaired files by using
% the Repair button in the GUI. Repairs can be done by interpolation
% between the nearest non-repaired scans (RECOMMNEDED),
% or despike interpolation using the
% immediate before and after scans, or inserting the mean image in place
% of the image to be repaired. Repairs will change the scans marked by
% red vertical lines. Scans marked by green vertical lines are unchanged,
% but will be added to the deweighting text file.
%     Repaired images have the prefix "v" and are written to the same
% folder as the input images. The input images are preserved in place.
% Unchanged images are copied to a "v" named file, so SPM can run with
% the "v" images. The program writes a file art_repaired.txt with a list
% of all repaired images, and art_deweighted.txt with a list of all images
% to be deweighted during estimation. Deweighting is recommended for both
% repaired images and margin images. Repairs bias contrasts lower, but
% informally, the effect seems small when fewer than 5-10% of scans are
% repaired. For deeper repairs, deweighting must be implemented in a
% batch script, e.g. as in the art_redo function.
%     A set of default thresholds is suggested ( 1.3% variation in global
% intensity, 0.5 mm/TR scan-to-scan motion.) Informally, these values
% are OK for small to moderate repairs. The thresholds can be reduced for
% good data (motion -> 0.3) , and should be raised for severely noisy data
% (motion -> 1.0). The values of the intensity and motion thresholds are
% linked and will change together. As a default, all images with total
% movement > 3 mm will also be marked for repair. For special situations,
% the outlier edit feature can mark additional scans. When motion
% regressors will be used, suggest setting the motion threshold to 0.5
% and not applying the margins.
%
% For MULTIPLE SESSIONS, we suggest realigning each session separately,
% and repairing each session separately. This approach agrees more
% with SPM standard practice, and clinical subjects often move between
% sessions. This approach differs from previous versions of this program.
%
% For batch scripts, use
% !!
% !! FORMAT art_global(Images, RealignmentFile, HeadMaskType, RepairType)
% !!
%    Images  = Full path name of images to be repaired.
%       For multiple sessions, all Images are in one array, e.g.SPM.xY.P
%    RealignmentFile = Full path name of realignment file
%       For multiple sessions, cell array with one name per session.
%    HeadMaskType  = 1 for SPM mask, = 4 for Automask
%    RepairType = 1 for ArtifactRepair alone (0.5 movement and add margin).
%               = 2 for Movement Adjusted images  (0.5 movement, no margin)
%               = 0 No repairs are done, bad scans are found.
%                   Listed in art_suspects.txt for motion adjustment.
%    Hardcoded actions:
%       Does repair, does not force repair of first scan.
% ----------------------------------------------------------------------
% May 2011 Modified by Eerik Puska
%    Takes thresholds as input arguments, no GUI.

% v2.5, May 2009 pkm  - adds SPM8, RepairType=0.

% v2.4, Mar 2009 pkm
%    Only one session allowed. Realign and repair each session separately.
%    Increment global intensity at 0.05%, instead of 0.1 in z-score.
%    Movement adaptive threshold now affects global threshold.
% v2.3, July 2008 PKM
%    Allows automatic adaptive threshold for movement
%    Compatible SPM5 and SPM2
%    Lowered default thresholds, to catch more deep breath artifacts.
% v2.2, July 2007 PKM
%    Also marks scan before a big movement for mv_out repair.
%    Allows two kinds of RepairType in batch, with and without margin.
%    Fixed logic for user masked mean
%    Adds user GUI option to repair all scans with movement > 3 mm.
%    small fix to subplots
% v2.1, Jan. 2007. PKM.
%    Prints copy of art_global figure in batch mode.
%    Allows multiple sessions in batch mode.
%       Images is full path name of all sessions, e.g. SPM.xY.P
%       Realignment file is one cell per session with full rp file name.
%    For SPM5, finds size of VY.dim.  ( Compatible with SPM2)
%
% Paul Mazaika, September 2006.
% This version replaces v.1 by Paul Mazaika and Sue Whitfield 2005,
%  originally derived from artdetect4.m, by Sue Whitfield,
%  Jeff Cooper, and Max Gray in 2002.

% -----------------------
% Initialize, begin loop
% -----------------------
bad_volumes=[];
pfig = [];
% Identify spm version
% spmv = spm('Ver'); spm_ver = 'spm2';
% if (strcmp(spmv,'SPM5') | strcmp(spmv,'SPM8b') | strcmp(spmv,'SPM8') )
%     spm_ver = 'spm5'; end

% ------------------------
% Default values for outliers
% ------------------------
% When std is very small, set a minimum threshold based on expected physiological
% noise. Scanners have about 0.1% error when fully in spec.
% Gray matter physiology has ~ 1% range, ~0.5% RMS variation from mean.
% For 500 samples, expect a 3-sigma case, so values over 1.5% are
% suspicious as non-physiological noise. Data within that range are not
% outliers. Set the default minimum percent variation to be suspicious...
Percent_thresh = PercentThresh; % default 1.3
%  Alternatively, deviations over 2*std are outliers, if std is not very small.
z_thresh = ZThresh;  % default 2.5
% Large intravolume motion may cause image reconstruction
% errors, and fast motion may cause spin history effects.
% Guess at allowable motion within a TR. For good subjects,
% would like this value to be as low as 0.3. For clinical subjects,
% set this threshold higher.
mv_thresh = MvThresh;  % default 0.5

if isempty(Mask)
    HeadMaskType = 4;
else
    HeadMaskType = 3;
end
% ------------------------
% Collect files
% ------------------------

% num_sess = 1;   % Only one session;
num_sess = 1;%size(RealignmentFile,1);
global_type_flag = HeadMaskType;
realignfile = 1;
P{1} = Images;
M{1} = RealignmentFile;
repair1_flag = 0;   % Only repair scan 1 when necessary
%repair1_flag = 1;   % Force repair scan 1  (Reisslab)
GoRepair = 1;       % Automatic Repair
if nargin == 7      % To stay backward compatible with v2.1
    if RepairType == 1
        GoRepair = 1;
    elseif RepairType == 2
        mv_thresh = MvThresh;
        GoRepair = 2;
    elseif RepairType == 0  % no repairs done, bad scans found.
        mv_thresh = MvThresh;
        GoRepair = 4;
    end
end

if global_type_flag==3
    maskY = Mask;%spm_read_vols(spm_vol(mask));
    %maskXYZmm = maskXYZmm(:,find(maskY==max(max(max(maskY)))));
    maskcount = sum(sum(sum(maskY)));  %  Number of voxels in mask.
    voxelcount = prod(size(maskY));    %  Number of voxels in 3D volume.
end

if global_type_flag == 4   %  Automask option
    fprintf('No mask given, computing new EPI mask\n')
    kk=size(Images);
    T=kk(4);
    Automask =  ones(kk(1:3));
    for t=1:T
        temp=squeeze(Images(:,:,:,t));
        Automask=Automask.*(temp>0.1*quantile(temp(:),.98));
    end
    maskcount = sum(sum(sum(Automask)));  %  Number of voxels in mask.
    voxelcount = prod(size(Automask));    %  Number of voxels in 3D volume.
end

mv_data = M{1};

nscans = size(Images,4);

g      = zeros(nscans,1);

fprintf('Calculating globals\n')

if global_type_flag == 3 % user masked mean
    %voxelcount = prod(size(Y));
    %vinv = inv(VY(1).mat);
    %[dummy, idx_to_mask] = intersect(XYZmm', maskXYZmm', 'rows');
    %maskcount = length(idx_to_mask);
    for i = 1:nscans
        Y = squeeze(Images(:,:,:,i));
        Y = Y.*maskY;
        %Y(idx_to_mask) = [];
        g(i) = mean(mean(mean(Y)))*voxelcount/maskcount;
        %g(i) = mean(Y(idx_to_mask));
    end
else   %  global_type_flag == 4  %  auto mask
    for i = 1:nscans
        Y = squeeze(Images(:,:,:,i));
        Y = Y.*Automask;
        if realignfile == 0
            output = art_centroid(Y);
            centroiddata(i,1:3) = output(2:4);
            g(i) = output(1)*voxelcount/maskcount;
        else     % realignfile == 1
            g(i) = mean(mean(mean(Y)))*voxelcount/maskcount;
        end
    end
    % If computing approximate translation alignment on the fly...
    %   centroid was computed in voxels
    %   voxel size is VY(1).mat(1,1), (2,2), (3,3).
    %   calculate distance from mean as our realignment estimate
    %   set rotation parameters to zero.
    if realignfile == 0    % change to error values instead of means.
        centroidmean = mean(centroiddata,1);
        for i = 1:nscans
            mv0data(i,:) = - centroiddata(i,:) + centroidmean;
        end
        % THIS MAY FLIP L-R  (x translation)
        mv_data(1:nscans,1) = mv0data(1:nscans,1)*VY(1).mat(1,1);
        mv_data(1:nscans,2) = mv0data(1:nscans,2)*VY(1).mat(2,2);
        mv_data(1:nscans,3) = mv0data(1:nscans,3)*VY(1).mat(3,3);
        mv_data(1:nscans,4:6) = 0;
    end
end

% in FSL, rotation and translation order is flipped
if(strcmp(cfg.prepro_suite,'fsl-fs'))
    mv_data=mv_data(:,[4,5,6,1,2,3]);
end
% Convert rotation movement to degrees
mv_data(:,4:6)= mv_data(:,4:6)*180/pi;

if global_type_flag==3
    fprintf('%g voxels in given mask\n', maskcount)
end
if global_type_flag==4
    fprintf('%g voxels in auto-generated mask\n', maskcount)
end

% ------------------------
% Compute default out indices by z-score, or by Percent-level is std is small.
% ------------------------
%  Consider values > Percent_thresh as outliers (instead of z_thresh*gsigma) if std is small.
gsigma = std(g);
gmean = mean(g);
pctmap = 100*gsigma/gmean;
mincount = Percent_thresh*gmean/100;
z_thresh = min( z_thresh, mincount/gsigma );
z_thresh = 0.1*round(z_thresh*10); % Round to nearest 0.1 Z-score value
zscoreA = (g - mean(g))./std(g);  % in case Matlab zscore is not available
glout_idx = (find(abs(zscoreA) > z_thresh))';

% ------------------------
% Compute default out indices from rapid movement
% ------------------------
%   % Rotation measure assumes voxel is 65 mm from origin of rotation.
if realignfile == 1 | realignfile == 0
    delta = zeros(nscans,1);  % Mean square displacement in two scans
    for i = 2:nscans
        delta(i,1) = (mv_data(i-1,1) - mv_data(i,1))^2 +...
            (mv_data(i-1,2) - mv_data(i,2))^2 +...
            (mv_data(i-1,3) - mv_data(i,3))^2 +...
            1.28*(mv_data(i-1,4) - mv_data(i,4))^2 +...
            1.28*(mv_data(i-1,5) - mv_data(i,5))^2 +...
            1.28*(mv_data(i-1,6) - mv_data(i,6))^2;
        delta(i,1) = sqrt(delta(i,1));
    end
end

% Also name the scans before the big motions (v2.2 fix)
deltaw = zeros(nscans,1);
for i = 1:nscans-1
    deltaw(i) = max(delta(i), delta(i+1));
end
delta(1:nscans-1,1) = deltaw(1:nscans-1,1);

% Adapt the threshold  (v2.3 fix)
if RepairType == 2 || GoRepair == 4
    delsort = sort(delta);
    if delsort(round(0.75*nscans)) > mv_thresh
        mv_thresh = min(1.0,delsort(round(0.75*nscans)));
        words = ['Automatic adjustment of movement threshold to ' num2str(mv_thresh)];
        disp(words)
        Percent_thresh = mv_thresh + PercentThresh - MvThresh;    % Modified
    end
end

mvout_idx = find(delta > mv_thresh)';

% Total repair list
out_idx = unique([mvout_idx glout_idx]);
if repair1_flag == 1
    out_idx = unique([ 1 out_idx]);
end
% Initial deweight list before margins
outdw_idx = out_idx;
% Initial clip list without removing large displacements
clipout_idx = out_idx;

bad_volumes = out_idx;
% -----------------------
% Draw initial figure
% -----------------------

try
    figure('Units', 'normalized', 'Position', [0.2 0.2 0.5 0.8],'Visible','off');
catch
    warning('Failed to create (invisible) figure, skipping...');
    return;
end
%figure
rng = max(g) - min(g);   % was range(g);
pfig = gcf;
% Don't show figure in batch runs
%if (nargin > 0); set(pfig,'Visible','off'); end

subplot(4,1,1);
plot(g);
%xlabel(['artifact index list [' int2str(out_idx') ']'], 'FontSize', 8, 'Color','r');
%ylabel(['Range = ' num2str(rng)], 'FontSize', 8);
ylabel('Global mean signal');
%xlabel('Vol');
title('Red vertical lines are to repair, Green vertical lines are to deweight');

% Add vertical exclusion lines to the global intensity plot
axes_lim = get(gca, 'YLim');
axes_height = [axes_lim(1) axes_lim(2)];
for i = 1:length(outdw_idx)   % Scans to be Deweighted
    line((outdw_idx(i)*ones(1, 2)), axes_height, 'Color', 'g');
end
if GoRepair == 2
    for i = 1:length(outdw_idx)   % Scans to be Deweighted
        line((outdw_idx(i)*ones(1, 2)), axes_height, 'Color', 'r');
    end
end
xlabel('Scan')

subplot(4,1,2);
%thresh_axes = gca;
%set(gca, 'Tag', 'threshaxes');
zscoreA = (g - mean(g))./std(g);  % in case Matlab zscore is not available
plot(abs(zscoreA));
ylabel('Std away from mean');
%xlabel('Vol');

thresh_x = 1:nscans;
thresh_y = z_thresh*ones(1,nscans);
line(thresh_x, thresh_y, 'Color', 'r');

%  Mark global intensity outlier images with vertical lines
axes_lim = get(gca, 'YLim');
axes_height = [axes_lim(1) axes_lim(2)];
for i = 1:length(glout_idx)
    line((glout_idx(i)*ones(1, 2)), axes_height, 'Color', 'r');
end
xlabel('Scan')

if realignfile == 1
    subplot(4,1,3);
    xa = [ 1:nscans];
    plot(xa,mv_data(:,1),'b-',xa,mv_data(:,2),'g-',xa,mv_data(:,3),'r-',...
        xa,mv_data(:,4),'r--',xa,mv_data(:,5),'b--',xa,mv_data(:,6),'g--');
    %plot(,'--');
    ylabel('Re-alignment param.');
    title('Translation [mm] solid lines, Rotation [deg] dashed lines');
    legend('x', 'y', 'z','pitch','roll','yaw','Location','EastOutside');
    h = gca;
    set(h,'Ygrid','on');
elseif realignfile == 0
    subplot(4,1,3);
    plot(mv0data(:,1:3));
    ylabel('Alignment (voxels)');
    %xlabel('Scans. VERY APPROXIMATE EARLY-LOOK translation in voxels.');
    legend('x', 'y', 'z',0);
    h = gca;
    set(h,'Ygrid','on');
end
xlabel('Scan')

subplot(4,1,4);   % Rapid movement plot
plot(delta);
ylabel('Motion [mm/TR]');
xlabel('Scan to scan movement [~mm], 65mm from origin assumed');
title('Fast motion (motion derivative)');
h = gca;
set(h,'Ygrid','on');

thresh_x = 1:nscans;
thresh_y = mv_thresh*ones(1,nscans);
line(thresh_x, thresh_y, 'Color', 'r');

% Mark all movement outliers with vertical lines
subplot(4,1,4)
axes_lim = get(gca, 'YLim');
axes_height = [axes_lim(1) axes_lim(2)];
for i = 1:length(mvout_idx)
    line((mvout_idx(i)*ones(1,2)), axes_height, 'Color', 'r');
end

if GoRepair == 1 || GoRepair == 2
    %figname = ['artglobal', subname, '.jpg'];
    figname = 'ArtRepair_report.png';
    try
        if ~exist([cfg.outpath,cfg.separator,'bramila'],'dir')
            mkdir([cfg.outpath,cfg.separator,'bramila']);
        end
        filepath = [cfg.outpath,cfg.separator,'bramila',cfg.separator,cfg.fileID,'_',figname];
        saveas(pfig,filepath);
    catch
        warning('Could not save figure!');
    end
    % Return to directory in use before writing the jpg.
    
    %art_repairvol(P);
end

end
