function icatb_gica_html_report(param_file)

global GICA_PARAM_FILE;

icatb_defaults;
global EXPERIMENTAL_TR;
global PARAMETER_INFO_MAT_FILE;
global CONVERT_Z;
global IMAGE_VALUES;
global THRESHOLD_VALUE;
global TIMECOURSE_POSTPROCESS;

if (~isempty(GICA_PARAM_FILE))
    param_file = GICA_PARAM_FILE;
end

if (~exist('param_file', 'var'))
    filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];
    param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', filterP);
end

outputDir = fileparts(param_file);
if (isempty(outputDir))
    outputDir = pwd;
end

load(param_file);

if (sesInfo.which_analysis == 2)
    
    %% ICASSO Plots
    icassoResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_icasso_results.mat']);
    try
        if (sesInfo.write_analysis_steps_in_dirs)
            icassoResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_ica_files', filesep, 'icasso_results.mat']);
        end
    catch
    end
    load(icassoResultsFile, 'sR');
    try
        icassoShow(sR, 'L', sesInfo.numComp, 'colorlimit', [.8 .9]);
    catch
    end
    
end

%% Mean Components
% Mean across all subjects and sessions is computed for each component
% a) Timecourse - Mean timecourse is converted to z-scores.
% b) Spectra - Timecourses spectra is computed for each data-set and
% averaged across sessions. Mean and standard error of mean is shown in the
% figure.
% c) Montage - Axial slices -40:4:72 are shown.
% d) Ortho slices - Ortho plot is shown for the peak voxel and coordinates
% are reported.
structFile = fullfile(fileparts(which('gift.m')), 'icatb_templates', 'ch2bet.nii');

compFileNaming = sesInfo.icaOutputFiles(1).ses(1).name;
currentFile = deblank(compFileNaming(1, :));

if ~exist(currentFile, 'file')
    [zipFileName, files_in_zip] = icatb_getViewingSet_zip(currentFile, [], 'real', sesInfo.zipContents);
    if (~isempty(zipFileName))
        icatb_unzip(regexprep(zipFileName, ['.*\', filesep], ''), fullfile(outputDir, fileparts(currentFile)));
    end
end

compFiles = icatb_rename_4d_file(fullfile(outputDir, compFileNaming));
icaTimecourse = icatb_loadICATimeCourse(compFiles);

postProcessFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_postprocess_results.mat']);
if (~exist(postProcessFile, 'file'))
    icatb_postprocess_timecourses(param_file);
end

load(postProcessFile);

if (sesInfo.numOfSess > 1)
    spectra_tc_all = reshape(mean(spectra_tc_all, 2), sesInfo.numOfSub, length(freq), sesInfo.numComp);
else
    spectra_tc_all = reshape(squeeze(spectra_tc_all), sesInfo.numOfSub, length(freq), sesInfo.numComp);
end

clear tc;

freq_limits = [0.1, 0.15];
try
    freq_limits = TIMECOURSE_POSTPROCESS.spectra.freq_limits;
catch
end

cmap = getCmap(IMAGE_VALUES);

for nF = 1:size(compFiles, 1)
    
    cn = deblank(compFiles(nF, :));
    
    [dd, fN, extn] = fileparts(cn);
    
    gH = icatb_getGraphics([fN, extn], 'graphics', 'imviewer', 'on');
    set(gH, 'resize', 'on');
    
    xOffSet = 0.05;
    yOffSet = 0.05;
    
    
    width2 = 0.5;
    height2 = 0.5;
    sh = axes('parent', gH, 'units', 'normalized', 'position', [0.01, yOffSet, width2, height2]);
    icatb_image_viewer(cn, 'structfile', structFile, 'image_values', IMAGE_VALUES, 'convert_to_zscores', CONVERT_Z, 'threshold', str2num(THRESHOLD_VALUE), ...
        'slices_in_mm', (-40:4:72), 'anatomical_view', 'axial', 'axesh', sh, 'colorbar', 0, 'labels', ' ');
    
    
    width = 0.4;
    height = 0.4;
    
    sh = axes('parent', gH, 'units', 'normalized', 'position', [width2 + 0.03, height2/2-(height/2), width, height]);
    plotStackedOrtho(cn, 'structfile', structFile, 'image_values', IMAGE_VALUES, 'convert_to_zscores', CONVERT_Z, 'threshold', str2num(THRESHOLD_VALUE), 'set_to_max_voxel', 1, ...
        'get_interp_data', 1, 'cmap', cmap, 'axesh', sh, 'colorbar', false, 'labels', 'Peak Coordinates (mm)', 'colorbar', true, 'colorbar_label', true);
    
    yPos = 0.65;
    width = 0.4;
    height = 0.25;
    sh = axes('parent', gH, 'units', 'normalized', 'position', [xOffSet, yPos, width, height]);
    tc =  icatb_detrend(icaTimecourse(:, nF));
    icatb_plotTimecourse('data', tc/std(tc), 'parent', sh, 'color', 'm', 'title', ['Component ', icatb_returnFileIndex(nF)], 'xlabelstr', 'Scans', 'ylabelstr', 'Z-scores');
    
    sh = axes('parent', gH, 'units', 'normalized', 'position', [xOffSet+width+0.1, yPos, width, height]);
    clear tc;
    tc.data = squeeze(spectra_tc_all(:, :, nF));
    dynamicrange = zeros(1, size(tc.data, 1));
    fALFF = dynamicrange;
    for nS = 1:length(dynamicrange)
        [dynamicrange(nS), fALFF(nS)] = icatb_get_spec_stats(tc.data(nS, :), freq, freq_limits);
    end
    tc.xAxis = freq;
    tc.isSpectra = 1;
    tc.xlabelStr = 'Frequency (Hz)';
    tc.ylabelStr = 'Power';
    tc.titleStr = [sprintf('Dynamic range: %0.3f, Power_L_F/Power_H_F: %0.3f', mean(dynamicrange), mean(fALFF));];
    icatb_plot_spectra(sh, tc);
    
    clear tc;
    
end



%% FNC correlations
% Functional network connectivity correlations are computed for each
% data-set and averaged across sessions.
if (sesInfo.numOfSess > 1)
    fnc_corrs_all = reshape(mean(fnc_corrs_all, 2), sesInfo.numOfSub, sesInfo.numComp, sesInfo.numComp);
else
    fnc_corrs_all = reshape(squeeze(fnc_corrs_all), sesInfo.numOfSub, sesInfo.numComp, sesInfo.numComp);
end

if (sesInfo.numOfSub > 1)
    fnc_corrs_all = squeeze(mean(fnc_corrs_all));
else
    fnc_corrs_all = squeeze(fnc_corrs_all);
end
fnc_corrs_all = icatb_z_to_r(fnc_corrs_all);
comps = (1:sesInfo.numComp)';
CLIM = max(abs(fnc_corrs_all(:)));
gH = figure('color', 'w');
set(gH, 'resize', 'on');
axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
icatb_plot_FNC(fnc_corrs_all, [-CLIM, CLIM], cellstr(num2str(comps)), (1:length(comps)), gH, [], axesH);
colormap(jet(64));
title('Average FNC Correlations', 'parent', axesH);



if (exist('files_in_zip', 'var') && ~isempty(files_in_zip))
    icatb_cleanupFiles(files_in_zip, outputDir);
end


function [sliceXY, sliceXZ, sliceYZ] = returnSlices(data, voxelcoords)

sliceXY = rot90(reshape(data(:, :, voxelcoords(end)), size(data, 1), size(data, 2)));
sliceXZ = rot90(reshape(data(:, voxelcoords(2), :), size(data, 1),size(data, 3)));
sliceYZ = rot90(reshape(data(voxelcoords(1), :, :), size(data, 2), size(data, 3)));

function data = stackData(slices, minVal)

[m1, n1] = size(slices{1});
[m2, n2] = size(slices{2});
[m3, n3] = size(slices{3});

maxSizeX = max([m1, m2, m3]);
maxSizeY = max([n1, n2, n3]);

data = minVal*ones(maxSizeX, [n1 + n2 + n3]);

e = 0;
for nS = 1:length(slices)
    tmp = slices{nS};
    s = e + 1;
    e = e + size(tmp, 2);
    inda = ceil(maxSizeX/2) - ceil(size(tmp, 1)/2) + 1;
    indb = inda + size(tmp, 1) - 1;
    data(inda:indb, s:e) = tmp;
end



function plotStackedOrtho(file_name, varargin)

icatb_defaults;
global FONT_COLOR;

useColorbar = 1;
threshold = 1;
convert_to_zscores = 'no';
image_values = 'positive';
labels = '';
cmap = hot(64);
structFile = fullfile (fileparts(which('gift.m')), 'icatb_templates', 'ch2bet.nii');

for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'structfile'))
        structFile = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'image_values'))
        image_values = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'cmap'))
        cmap = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'threshold'))
        threshold = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'convert_to_zscores'))
        convert_to_zscores = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'axesh'))
        sh = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'colorbar'))
        useColorbar = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'labels'))
        labels = varargin{n + 1};
    end
end

hD = icatb_orth_views(file_name, 'structfile', structFile, 'image_values', image_values, 'convert_to_zscores', convert_to_zscores, 'threshold', threshold, 'set_to_max_voxel', 1, ...
    'get_interp_data', 1, 'cmap', cmap);
if (~exist('sh', 'var'))
    sh = gca;
end
[sliceXY, sliceXZ, sliceYZ] = returnSlices(hD.data, hD.maxVoxelPos);
tmp = stackData({sliceYZ, sliceXZ, sliceXY}, 0);
imagesc(tmp, [1, 200]);
colormap(hD.cmap);
axis(sh, 'image');
set(sh, 'Ytick', []);
set(sh, 'Xtick', [])
realCoords = (hD.maxVoxelPos - hD.voxelOrigin).*hD.VOX;

str = charC (labels, ['(', num2str(realCoords(1)), ',', num2str(realCoords(2)), ',', num2str(realCoords(3)), ')']);
title(str, 'parent', sh, 'horizontalalignment', 'center', 'fontweight', 'bold');
if (useColorbar)
    ch = colorbar('location', 'southoutside');
    set(ch, 'units', 'normalized');
    chPos = get(sh, 'position');
    chPos(1) = chPos(1) + 0.5*chPos(3) - 0.2;
    chPos(2) = chPos(2) - 0.02;
    chPos(3) = 0.4;
    chPos(4) = 0.025;
    set(ch, 'position', chPos);
    set(ch, 'xLim', [1, 100]);
    xTicks = get(ch, 'xTick');
    set(ch, 'yTick', []);
    good_inds = abs(hD.data) > eps;
    minICAVal = hD.minICAIM;
    maxICAVal = hD.maxICAIM;
    set(ch, 'xTick', [xTicks(1), xTicks(end)]);
    set(ch, 'xTicklabel', cellstr(char(num2str(minICAVal, '%0.1f'), num2str(maxICAVal, '%0.1f'))));
    set(ch, 'XColor', FONT_COLOR, 'YColor', FONT_COLOR);
end

function cmap = getCmap(image_values)

load icatb_colors coldhot_sensitive;

returnValue = strmatch(lower(image_values), {'positive and negative', 'positive', 'absolute value', 'negative'}, 'exact');

if (returnValue == 1)
    cmap = coldhot_sensitive(1:4:end, :);
elseif (returnValue == 4)
    cmap = coldhot_sensitive(1:128, :);
    cmap = cmap(1:2:end, :);
else
    cmap = coldhot_sensitive(129:end, :);
    cmap = cmap(1:2:end, :);
end

function c = charC(a, b)

len = max([length(a), length(b)]);

c = repmat(' ', 2, len);

s = ceil(len/2 - length(a)/2) + 1;
c(1, s:s+length(a) - 1) = a;
s = ceil(len/2 - length(b)/2) + 1;
c(2, s:s + length(b) - 1) = b;