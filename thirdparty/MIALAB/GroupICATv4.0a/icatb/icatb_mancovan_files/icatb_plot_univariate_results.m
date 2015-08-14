function figs = icatb_plot_univariate_results(mancovanInfo, timeNo)
%% Plot Univariate results
%

icatb_defaults;
global UI_FONTNAME;
global UI_FS;

if (~exist('timeNo', 'var'))
    timeNo = 0;
end

outputDir = mancovanInfo.outputDir;
comps = mancovanInfo.comps;
results = mancovanInfo.outputFiles;

if (isempty(mancovanInfo.covariatesToPlot))
    error('Covariates of interest are not selected');
end

covariatesToPlot = mancovanInfo.covariatesToPlot;
covariatesToPlot = cellstr(covariatesToPlot);

figs = [];

for nF = 1:length(results)
    featureName = results(nF).feature_name;
    if (strcmpi(featureName, 'spatial maps'))
        figTitle = 'Univariate Results (Spatial maps)';
        nCols = 4;
        structFile = mancovanInfo.structFile;
        gH = icatb_getGraphics(figTitle, 'graphics', 'univariate_results_spatial_maps', 'off');
        disp('Writing out composite images of covariates of interest ...');
        img_files = write_composite_sig_effects(mancovanInfo.HInfo(1), covariatesToPlot, results(nF).filesInfo.result_files, mancovanInfo.prefix, outputDir, ...
            mancovanInfo.userInput.p_threshold, mancovanInfo.threshdesc, timeNo);
        countPlot = 0;
        for nCov = 1:length(covariatesToPlot)
            countPlot = countPlot + 1;
            axesH(1) = subplot(length(covariatesToPlot), nCols, countPlot);
            countPlot = countPlot + 1;
            axesH(2) = subplot(length(covariatesToPlot), nCols, countPlot);
            countPlot = countPlot + 1;
            axesH(3) = subplot(length(covariatesToPlot), nCols, countPlot);
            tmp = icatb_loadData(img_files{nCov});
            if (any(abs(tmp(:)) > eps))
                tmpTitle = ['Significant Effects Of ', covariatesToPlot{nCov}, ' (p < ', num2str(mancovanInfo.userInput.p_threshold), ')'];
            else
                tmpTitle = ['No Significant Effects Of ', covariatesToPlot{nCov}, ' (p < ', num2str(mancovanInfo.userInput.p_threshold), ')'];
            end
            clear tmp;
            
            if (timeNo > 0)
                tmpTitle = [tmpTitle, '(Time', num2str(timeNo), ')'];
            end
            
            hD = icatb_orth_views(img_files{nCov}, 'structfile', structFile, 'image_values', mancovanInfo.image_values, 'convert_to_zscores', 'no', 'labels', ...
                tmpTitle, 'get_interp_data', 1, 'set_to_max_voxel', 1);
            cmap1 = hD.cmap;
            colormap(cmap1);
            hD.coords = hD.maxVoxelPos;
            hD.currentFigure = gH;
            ch1 = drawSlices(axesH, hD);
            set(gH, 'visible', 'on');
            countPlot = countPlot + 1;
            
            sh = subplot(length(covariatesToPlot), nCols, countPlot);
            S = gather_univariate_effects(mancovanInfo, nF, covariatesToPlot{nCov}, timeNo);
            S.colorbarLabel = 'Fraction Of Component';
            
            
            vars2 = {S, sh, cmap1(1:ceil(size(cmap1, 1)/2), :)};
            ch2 = plot_univariate_effects(vars2{:});
            
            YLim = get(ch2, 'YLim');
            YLim(2) = YLim(1) + 0.5*diff(YLim);
            set(ch2, 'YLim', YLim);
            
            
            YTicksCB = get(ch2, 'YTick');
            
            if (length(YTicksCB) > 3)
                YTicksCB = [YTicksCB(1), median(YTicksCB), YTicksCB(end)];
                set(ch2, 'YTick', YTicksCB);
                set(ch2, 'YTickLabel', {'1', '0', '1'});
            end
            
            childAxesH = get(findobj(gH, 'Type', 'Axes'), 'children');
            %if (iscell(childAxesH))
            %    childAxesH = cell2mat(childAxesH);
            %end
            for nChildH = 1:length(childAxesH)
                if (iscell(childAxesH))
                    set(childAxesH{nChildH}, 'hittest', 'off');
                else
                    set(childAxesH(nChildH), 'hittest', 'off');
                end
            end
            set([axesH, ch1], 'ButtonDownFcn', {@openOrthoViews, hD});
            set([sh, ch2], 'ButtonDownFcn', {@openUnivEffects, vars2, ch2});
            
            clear hD vars2;
            
        end
        
        
    elseif (strcmpi(featureName, 'timecourses spectra'))
        figTitle = 'Univariate Results (Spectra)';
        gH = icatb_getGraphics(figTitle, 'graphics', 'univariate_results_sepctra', 'on');
        load icatb_colors coldhot;
        load(fullfile(mancovanInfo.outputDir, mancovanInfo.outputFiles(nF).filesInfo.result_files{1}), 'freq');
        countPlot = 0;
        nCols = 2;
        for nCov = 1:length(covariatesToPlot)
            S = gather_univariate_stats(mancovanInfo, nF, covariatesToPlot{nCov}, timeNo);
            countPlot = countPlot + 1;
            sh = subplot(length(covariatesToPlot), nCols, countPlot);
            axesH = sh;
            cmap1 = coldhot;
            S.c = ['Significant Effects Of ', covariatesToPlot{nCov}, ' (p < ', num2str(mancovanInfo.userInput.p_threshold), ')'];
            
            if (timeNo > 0)
                S.c = [S.c, '(Time', num2str(timeNo), ')'];
            end
            
            
            vars1 = {S, freq, mancovanInfo.comps, mancovanInfo.userInput.p_threshold, mancovanInfo.threshdesc, cmap1, 'Component', 'Frequency (Hz)', sh};
            ch1 = plotUnivStats (vars1{:});
            
            countPlot = countPlot + 1;
            sh = subplot(length(covariatesToPlot), nCols, countPlot);
            S = gather_univariate_effects(mancovanInfo, nF, covariatesToPlot{nCov}, timeNo);
            S.colorbarLabel = 'Fraction Of Spectrum';
            vars2 = {S, sh, cmap1};
            ch2 = plot_univariate_effects(vars2{:});
            
            YTicksCB = get(ch2, 'YTick');
            
            if (length(YTicksCB) > 3)
                YTicksCB = [YTicksCB(1), median(YTicksCB), YTicksCB(end)];
                set(ch2, 'YTick', YTicksCB);
                set(ch2, 'YTickLabel', {'1', '0', '1'});
            end
            
            childAxesH = get(findobj(gH, 'Type', 'Axes'), 'children');
            for nChildH = 1:length(childAxesH)
                if (iscell(childAxesH))
                    set(childAxesH{nChildH}, 'hittest', 'off');
                else
                    set(childAxesH(nChildH), 'hittest', 'off');
                end
            end
            % if (iscell(childAxesH))
            %    childAxesH = cell2mat(childAxesH);
            %end
            % set(childAxesH, 'hittest', 'off');
            set([axesH, ch1], 'ButtonDownFcn', {@openUnivStats, vars1});
            set([sh, ch2], 'ButtonDownFcn', {@openUnivEffects, vars2, ch2});
            
            clear vars1 vars2;
            
            axis(sh);
            
        end
        
    else
        figTitle = 'Univariate Results (FNC Correlations)';
        gH = icatb_getGraphics(figTitle, 'graphics', 'univariate_results_fnc', 'on');
        load icatb_colors coldhot;
        countPlot = 0;
        nCols = 1;
        for nCov = 1:length(covariatesToPlot)
            S = gather_univariate_stats(mancovanInfo, nF, covariatesToPlot{nCov}, timeNo);
            countPlot = countPlot + 1;
            sh = subplot(length(covariatesToPlot), nCols, countPlot);
            axesH = sh;
            cmap1 = coldhot;
            S.c = ['Significant Effects Of ', covariatesToPlot{nCov}, ' (p < ', num2str(mancovanInfo.userInput.p_threshold), ')'];
            if (timeNo > 0)
                S.c = [S.c, '(Time', num2str(timeNo), ')'];
            end
            M = S.logp ;
            M = icatb_vec2mat(M, 1);
            S.logp = M;
            vars = {S, mancovanInfo.comps, mancovanInfo.comps, mancovanInfo.userInput.p_threshold, mancovanInfo.threshdesc, cmap1, 'Component', 'Component', sh, 1};
            ch1 = plotUnivStats (vars{:});
            
            childAxesH = get(findobj(gH, 'Type', 'Axes'), 'children');
            for nChildH = 1:length(childAxesH)
                if (iscell(childAxesH))
                    set(childAxesH{nChildH}, 'hittest', 'off');
                else
                    set(childAxesH(nChildH), 'hittest', 'off');
                end
            end
            %             if (iscell(childAxesH))
            %                 childAxesH = cell2mat(childAxesH);
            %             end
            % set(childAxesH, 'hittest', 'off');
            set([axesH, ch1], 'ButtonDownFcn', {@openUnivStats, vars});
            
            clear vars;
            axis(sh, 'square');
            
            [levelNames, levelMeans, titleStr] = plotContextMenu(mancovanInfo,  covariatesToPlot{nCov}, nF, timeNo);
            if (~isempty(levelMeans))
                contextMenuH = uicontextmenu;
                item1 = uimenu(contextMenuH, 'Label', 'Level Means');
                set(sh, 'uicontextmenu', contextMenuH);
                lData.featureName = featureName;
                lData.levelNames = levelNames;
                lData.levelMeans = levelMeans;
                lData.title = titleStr;
                set(item1, 'userdata', lData);
                set(item1, 'callback', {@plotLevelMeans, mancovanInfo});
            end
            
        end
        
    end
    
    
    axisH = axes('Parent', gH, 'position', [0 0 1 1], 'visible', 'off');
    xPos = 0.5; yPos = 0.97;
    titleColor = 'c';
    text(xPos, yPos, figTitle, 'color', titleColor, 'fontweight', 'bold', 'fontsize', UI_FS + 2, 'HorizontalAlignment', 'center', 'FontName', UI_FONTNAME, 'parent', axisH);
    
    figs(length(figs) + 1).H = gH;
    
    
end


function S = gather_univariate_stats(mancovanInfo, index, term, timeNo)
%% Gather univariate stats
%

result_files = mancovanInfo.outputFiles(index).filesInfo.result_files;
outputDir = mancovanInfo.outputDir;
thresh = mancovanInfo.userInput.p_threshold;
featureName = mancovanInfo.outputFiles(index).feature_name;

for ii = 1:length(result_files)
    if (timeNo > 0)
        load(fullfile(outputDir, result_files{ii}), 'time');
        UNI = time.UNI{timeNo};
    else
        load(fullfile(outputDir, result_files{ii}), 'UNI');
    end
    %mIND = strmatch(term, UNI.tests, 'exact');
    [mIND, con_no] = getTermAndConNum(term, UNI);
    if (strcmpi(featureName, 'timecourses spectra'))
        load(fullfile(outputDir, result_files{ii}), 'freq');
        xdim = length(freq);
    else
        load(fullfile(outputDir, result_files{ii}), 'fnc_corrs');
        xdim = size(fnc_corrs, 2);
        clear fnc_corrs;
    end
    
    if (ii == 1)
        STATS = zeros(length(result_files), xdim);
        STATS2 = STATS;
    end
    
    if ~isempty(mIND)
        ps = UNI.p{mIND}(con_no, :);
        ts = UNI.t{mIND}(con_no, :);
        if (strcmpi(mancovanInfo.threshdesc, 'fdr'))
            p_masked = icatb_fdr(ps, thresh);
        else
            p_masked = thresh;
        end
        good_inds = ps < p_masked;
        STATS(ii, good_inds) = -log10(ps(good_inds) + eps).*sign(ts(good_inds));
        STATS2(ii, good_inds) = ps(good_inds);
        conname = UNI.stats{mIND}.Contrast{con_no};
    end
end


if (~exist('conname', 'var'))
    conname = term;
end

S.logp = STATS;
S.p = STATS2;
S.x = term;
S.c = conname;
S.y = cellstr(num2str(mancovanInfo.comps(:)))';

function S = gather_univariate_effects(mancovanInfo, index, term, timeNo)
%% Gather univariate effects
%

result_files = mancovanInfo.outputFiles(index).filesInfo.result_files;
outputDir = mancovanInfo.outputDir;
thresh = mancovanInfo.userInput.p_threshold;
featureName = mancovanInfo.outputFiles(index).feature_name;

for ii = 1:length(result_files)
    
    if (timeNo > 0)
        load(fullfile(outputDir, result_files{ii}), 'time');
        UNI = time.UNI{timeNo};
    else
        load(fullfile(outputDir, result_files{ii}), 'UNI');
    end
    
    %mIND = strmatch(term, UNI.tests, 'exact');
    [mIND, con_no] = getTermAndConNum(term, UNI);
    
    if ~isempty(mIND)
        ps = UNI.p{mIND}(con_no, :) + eps;
        ts = UNI.t{mIND}(con_no, :);
        betas = getBetaWeights(UNI.stats{mIND}, con_no);
        %betas = UNI.stats{mIND}.B(con_no + 1, :);
        if (strcmpi(mancovanInfo.threshdesc, 'fdr'))
            [p_masked, p] = icatb_fdr(ps, thresh);
        else
            p_masked = thresh;
        end
        a = ps < p_masked;
        if length(find(a == 1))<1
            D = [];
        else
            switch lower(featureName)
                case 'spatial maps'
                    load(fullfile(outputDir, result_files{ii}), 'mask_ind');
                    B = zeros(mancovanInfo.HInfo.dim(1:3));
                    B(mask_ind) = a;
                    D = get_clusters(B, 0);
                    
                    for jj = 1:length(D)
                        D{jj}.ICAmask = D{jj}.mask;
                        [TF,LOC] = ismember(D{jj}.mask, mask_ind);
                        D{jj}.mask = LOC;
                        D{jj}.fn = D{jj}.n/length(a);
                        
                    end
                    
                case 'timecourses spectra'
                    D = get_spec_clusters(a, 0);
                case 'fnc correlations'
                    D = [];
                otherwise
                    fprintf('data type unknown')
                    S = [];
                    return;
            end
            
            for jj = 1:length(D)
                if sign(mean(betas(D{jj}.mask))) == 1
                    %girl cluster
                    D{jj}.sign = 1;
                    
                else
                    D{jj}.sign = -1;
                end
                D{jj}.beta = mean(betas(D{jj}.mask));
                D{jj}.allbeta = betas(D{jj}.mask);
                D{jj}.allt = ts(D{jj}.mask);
                D{jj}.allp = -log10(ps(D{jj}.mask)).*sign(ts(D{jj}.mask));
            end
            conname = UNI.stats{mIND}.Contrast{con_no};
        end
    else
        D = [];
    end
    D = combine_clusters(D);
    S.B{ii} = D;
    
    clear UNI mask_ind;
    
end

if (~exist('conname', 'var'))
    conname = term;
end

S.x = term;
S.c = conname;
S.y = cellstr(num2str(mancovanInfo.comps(:)))';


function D = get_clusters(B, minclustersize)

C = zeros(size(B)); C(find(B > 0)) = 1;
[L,NUM] = icatb_spm_bwlabel(C,26);
D = [];
cnum = 0;
for ii = 1:NUM
    IND = find(L == ii);
    numvox = length(IND);
    fprintf('Cluster %d, size = %d\n',ii,numvox )
    if numvox >= minclustersize;
        cnum = cnum+1;
        D{cnum}.cluster = zeros(size(B));
        D{cnum}.cluster(IND) = 1;
        D{cnum}.mask = IND;
        D{cnum}.n = numvox;
    else
        
    end
end

function D = get_spec_clusters(a, minclustersize)

bfile = which('bwlabel.m');
if (isempty(bfile))
    [L, nnn] = icatb_spm_bwlabel(double(a), 6);
else
    L = bwlabel(a, 4);
end
NUM = length(unique(L))-1;
cnum = 0;
D = [];
for ii = 1:NUM
    IND = find(L == ii);
    numvox = length(IND);
    fprintf('Cluster %d, size = %d\n',ii,numvox )
    if numvox > minclustersize;
        cnum = cnum+1;
        D{cnum}.cluster = zeros(size(a));
        D{cnum}.cluster(IND) = 1;
        D{cnum}.mask = IND;
        D{cnum}.n = numvox;
        D{cnum}.fn = numvox/length(a);
    else
        
    end
end


function [sliceXY, sliceXZ, sliceYZ] = returnSlices(data, voxelcoords)
%% Slices (XY, XZ, YZ)

sliceXY = reshape(data(:, :, voxelcoords(end)), size(data, 1), size(data, 2));
sliceXZ = reshape(data(:, voxelcoords(2), :), size(data, 1),size(data, 3));
sliceYZ = reshape(data(voxelcoords(1), :, :), size(data, 2), size(data, 3));


function ch = drawSlices(axesH, handles_data)
%% Draw Slices
%

icatb_defaults;
global FONT_COLOR;

composite_data = handles_data.data;
CLIM = handles_data.CLIM;
minInterval = handles_data.minInterval;
maxInterval = handles_data.maxInterval;
coords = handles_data.coords;
structVol = handles_data.structVol;
labels = handles_data.labels;
minICAIM = handles_data.minICAIM;
maxICAIM = handles_data.maxICAIM;
gH = handles_data.currentFigure;

[sliceXY, sliceXZ, sliceYZ] = returnSlices(composite_data, coords);

sh = axesH(1);
realCoords = (coords - handles_data.voxelOrigin).*handles_data.VOX;
%realCoords = icatb_voxel_to_real(structVol, coords);
set(sh, 'units', 'normalized');
plotImage(sh, sliceYZ, CLIM);
xlabel(['X = ', num2str(realCoords(1)), ' mm'], 'parent', sh);
firstPos = get(sh, 'position');

sh = axesH(2);
set(sh, 'units', 'normalized');
plotImage(sh, sliceXZ, CLIM);
xlabel(['Y = ', num2str(realCoords(2)), ' mm'], 'parent', sh);

middleAxes = sh;

sh = axesH(3);
set(sh, 'units', 'normalized');
plotImage(sh, sliceXY, CLIM);
xlabel(['Z = ', num2str(realCoords(3)), ' mm'], 'parent', sh);
lastPos = get(sh, 'position');

% Plot label
title(labels{1}, 'parent', middleAxes);
width = lastPos(1) + lastPos(3) - firstPos(1);
%         xPos = -width/2 - 0.05/2;
%         yPos = 0.05/2;
% text(xPos, yPos, labels{ncoordss}, 'parent', sh, 'units', 'normalized');

% Plot horizontal colorbar
colorBarWidth = 0.3;
colorBarHeight = 0.016;
pos = [firstPos(1) + width/2 - colorBarWidth/2, lastPos(2)-0.06, colorBarWidth, colorBarHeight];
ch = colorbar('horiz');
set(ch, 'units', 'normalized');
set(ch, 'position', pos);
%ch = colorbar(pos);
set(ch, 'xLim', [minInterval, maxInterval]);
xTicks = get(ch, 'xTick');
set(ch, 'yTick', []);
set(ch, 'xTick', [xTicks(1), xTicks(end)]);
set(ch, 'xTicklabel', cellstr(num2str([minICAIM;maxICAIM])));
set(ch, 'XColor', FONT_COLOR, 'YColor', FONT_COLOR);
xlabel(ch, '-sign(t) log_1_0(p-value)', 'Interpreter', 'tex');
set(ch, 'tag', 'colorbar');
set(axesH, 'XColor', get(gH, 'color'), 'YColor', get(gH, 'color'));

xlabelH = get(axesH, 'XLabel');
% if (iscell(xlabelH))
%     xlabelH = cell2mat(xlabelH);
% end
for nChildH = 1:length(xlabelH)
    if (iscell(xlabelH))
        set(xlabelH{nChildH}, 'color', FONT_COLOR);
    else
        set(xlabelH(nChildH), 'color', FONT_COLOR);
    end
end
%set(xlabelH, 'color', FONT_COLOR);

set(gH, 'visible', 'on');


function CB = plot_univariate_effects(S, axesH, cmap)
% first use:
%S  = gather_univariate_effects(type, term, contrast)

if (~exist('axesH', 'var'))
    figure;
    axesH = gca;
end

if (~isfield(S, 'colorbarLabel'))
    colorbarLabel = 'Data Fraction';
else
    colorbarLabel = S.colorbarLabel;
end


s = [1 -1];

BETAS = nan*ones(2, length(S.B));
Ns = nan*ones(2, length(S.B));
s= [1, -1];
for jj= 1:length(S.B)
    if ~isempty(S.B{jj})
        for ii = 1:length(s)
            [mind] = find(S.B{jj}.sign == s(ii));
            if ~isempty(mind)
                BETAS(ii,jj) = weighted_mean( S.B{jj}.beta(mind), S.B{jj}.n(mind));
                Ns(ii,jj) = sum(S.B{jj}.fn(mind));
            end
        end
    end
end

ydesc = 'Average Beta';
Y{1} = BETAS;
Y{2} = Ns;
CB = make_bars(BETAS, Ns, S.y, ydesc, axesH, colorbarLabel, cmap);

if ( length(find(cellfun('isempty', S.B) == 1)) ~= length(S.B))
    title(['Effect Sizes Of ', S.x], 'parent', axesH);
else
    title(['No Significant Effect Sizes Of ', S.x], 'parent', axesH);
end


function CB = make_bars(Y, YC, xdesc, ydesc, axesH, colorbarLabel, cmap)

icatb_defaults;
global FONT_COLOR;

% B1 = barh(1:length(xdesc),Y(1,:)) ; hold on
% set_bar_color(B1, YC(1,:), 'r')
% B2 = barh(1:length(xdesc),Y(2,:)) ;
% set_bar_color(B2, YC(2,:), 'b')

for ii = 1:length(xdesc)
    B1 = barh(ii,Y(1,ii), 'parent', axesH) ; hold on
    CM1 = set_bar_color(B1, YC(1,ii), 'r', cmap);
    B2 = barh(ii,Y(2,ii), 'parent', axesH) ;
    CM2 = set_bar_color(B2, YC(2,ii), 'b', cmap);
end
%colormap(CM2);

CB = colorbar;

ylabel(CB, colorbarLabel, 'parent', CB);

set(axesH, 'YTick',1:length(xdesc),'YTickLabel',xdesc, 'TickDir', 'out')

axis(axesH, 'tight');
ax = axis;
axis ij;
mx = max(abs(ax(1:2)));
axis([-1.1*mx 1.1*mx 0.5 length(xdesc)+.5 ]);

ylabel('Component', 'parent', axesH);
xlabel(ydesc, 'parent', axesH);

set(axesH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);

if (icatb_get_matlab_version >= 2014)
    set(CB, 'color', FONT_COLOR);
end

function w = weighted_mean(x, nx)
w = sum(x.*nx)/sum(nx);

function CM = set_bar_color(BH, y2, ctype, CMfull)

minx = 0;
maxx = 1;
colorlength = ceil(size(CMfull, 1)/2);
if strmatch(ctype, 'r')
    CM = CMfull(colorlength+1:end,:);
else
    CM = CMfull(1:colorlength,:);
    CM = flipud(CM);
end

X = linspace(minx, maxx, colorlength);
for ii = 1:length(BH)
    if isnan(y2(ii))
        set(BH(ii),'facecolor','w')
    else
        for jj = 1:3
            YI(jj) = interp1(X,CM(:,jj),y2(ii));
        end
        set(BH(ii),'facecolor',YI)
    end
end


function outfiles = write_composite_sig_effects(V, terms, result_files, prefix, outputDir, thresh, threshdesc, timeNo)
%% Write composite files for covariates of interest
%

outfiles = repmat({''}, 1, length(terms));

for nT = 1:length(terms)
    
    term = terms{nT};
    
    B = zeros(V(1).dim(1:3));
    
    for ii = 1:length(result_files)
        if (timeNo > 0)
            load(fullfile(outputDir, result_files{ii}), 'time', 'mask_ind');
            UNI = time.UNI{timeNo};
        else
            load(fullfile(outputDir, result_files{ii}), 'UNI', 'mask_ind');
        end
        %mIND = strmatch(term, UNI.tests, 'exact');
        [mIND, con_no] = getTermAndConNum(term, UNI);
        tmp_mask_ind = mask_ind;
        
        if ~isempty(mIND)
            ps = UNI.p{mIND}(con_no, :);
            if (strcmpi(threshdesc, 'fdr'))
                p_masked = icatb_fdr(ps, thresh);
            else
                p_masked = thresh;
            end
            bad_inds = ps > p_masked;
            tmp_mask_ind(bad_inds) = [];
            ps = ps(bad_inds == 0);
            ps = -log10(ps + eps);
            ts = UNI.t{mIND}(con_no, :);
            d =  ps.*sign(ts(bad_inds == 0));
        else
            d = zeros(1, length(mask_ind));
        end
        
        TEMP = zeros(size(B));
        TEMP(tmp_mask_ind) = d;
        gmask = find(abs(TEMP) > abs(B));
        B(gmask) = TEMP(gmask);
        clear I mIND TEMP
    end
    
    tmpPath = fileparts(result_files{ii});
    
    if (timeNo > 0)
        V.fname = fullfile(outputDir, tmpPath, [prefix, '_sm_', term, '_sig_effects_time', num2str(timeNo), '.img']);
    else
        V.fname = fullfile(outputDir, tmpPath, [prefix, '_sm_', term, '_sig_effects.img']);
    end
    
    V.n(1) = 1;
    icatb_write_vol(V, B);
    outfiles{nT} = V.fname;
    
end


function subH = plotImage(subH, data, CLIM)
%% Function to plot the image at the specified position
%

image(rot90(data), 'parent', subH, 'CDataMapping', 'scaled');
set(subH, 'units', 'normalized');
imageAxisPos = get(subH, 'position');
oldWidth = imageAxisPos(3);
oldHeight = imageAxisPos(4);
yAxisRatio = imageAxisPos(4)/imageAxisPos(3);
xAxisRatio = imageAxisPos(3)/imageAxisPos(4);
if(yAxisRatio>1)
    yAxisRatio = 1;
else
    xAxisRatio = 1;
end
newWidth = imageAxisPos(3)*yAxisRatio;
newHeight = imageAxisPos(4)*xAxisRatio;

if (newHeight < oldHeight)
    imageAxisPos(2) = imageAxisPos(2) + 0.5*(oldHeight - newHeight);
end

imageAxisPos = [imageAxisPos(1), imageAxisPos(2), newWidth, newHeight];
set(subH, 'position', imageAxisPos);
%setImagePos(subH, imagePos);
set(subH, 'clim', CLIM); % set the axis positions to the specified
%axis(subH, 'square');
set(subH, 'XTick', []);
set(subH, 'XTickLabel', []);
set(subH, 'YTick', []);
set(subH, 'YTickLabel', []);



function CLim = newclim(BeginSlot,EndSlot,CDmin,CDmax,CmLength)
% Convert slot number and range to percent of colormap

PBeginSlot    = (BeginSlot - 1) / (CmLength - 1);
PEndSlot      = (EndSlot - 1) / (CmLength - 1);
PCmRange      = PEndSlot - PBeginSlot;

% Determine range and min and max of new CLim values
DataRange     = CDmax - CDmin;
ClimRange     = DataRange / PCmRange;
NewCmin       = CDmin - (PBeginSlot * ClimRange);
NewCmax       = CDmax + (1 - PEndSlot) * ClimRange;
CLim          = [NewCmin, NewCmax];


function C = plotUnivStats(T, f, comps, thresh, threshdesc, cmap, yLabelStr, xLabelStr, axesH, fnc)

icatb_defaults;
global FONT_COLOR;

A = T.logp;
ncomp = size(A,1);
A(isnan(A)) = 0;

pv = T.p(find(isnan(T.p) == 0));
if (strcmpi(threshdesc, 'fdr'))
    p_masked  = icatb_fdr( pv, thresh);
else
    p_masked = thresh;
end

if (strcmpi(threshdesc, 'fdr'))
    disp(p_masked)
    fdrlim = -log10(p_masked);
end

if (exist('fnc', 'var') && fnc)
    M = imagesc(1:length(comps),1:length(comps),A);
    set(axesH,'XTick',1:length(comps),'XTickLabel',T.y, 'TickDir', 'out')
    
else
    M = imagesc(f,1:length(comps),A);
end

CLIM = [-max(abs(A(:))) max(abs(A(:)))];

if (~isempty(find(CLIM ~= 0)))
    set(axesH, 'CLIM', CLIM);
end

set(axesH,'YTick',1:length(comps),'YTickLabel',T.y, 'TickDir', 'out')
%set(gca, 'CLim', [log10(eps) -log10(eps)])

colormap(cmap);
C = colorbar;
yt = get(C, 'YTick');
if (strcmpi(threshdesc, 'fdr'))
    set(C, 'YTick', sort([yt, -fdrlim fdrlim]));
end
ylabel(C, '-sign(t) log_1_0(p-value)', 'Interpreter', 'tex');
if (length(find(A==0)) == numel(A))
    T.c = ['No ', T.c];
end
title(T.c);
ylabel(yLabelStr, 'parent', axesH);
xlabel(xLabelStr, 'parent', axesH);

set(axesH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);

if (icatb_get_matlab_version >= 2014)
    set(C, 'color', FONT_COLOR);
end

function E = combine_clusters(D)

if isempty(D)
    E = D;
    return
end

nC = length(D); %number of clusters

% D =
%     cluster: [53x63x46 double]
%        mask: [29x1 double]
%           n: 29
%          fn: 0.0062
%        sign: 1
%        beta: 0.4424


for ii = 1:nC
    
    %        eval(['E.' F{jj} '{ii} = D{ii}.' F{jj} ';']);
    E.mask{ii} = D{ii}.mask;
    E.n(ii) = D{ii}.n;
    if isfield(D{1}, 'fn')
        E.fn(ii) = D{ii}.fn;
    end
    if isfield(D{1}, 'sign')
        E.sign(ii) = D{ii}.sign;
    end
    if isfield(D{1}, 'beta')
        E.beta(ii) = D{ii}.beta;
    end
    
    if isfield(D{1}, 'ICAmask')
        E.ICAmask{ii} = D{ii}.ICAmask;
    end
    
    if isfield(D{1}, 'allbeta')
        E.allbeta{ii} = D{ii}.allbeta;
    end
    
    if isfield(D{1}, 'allt')
        E.allt{ii} = D{ii}.allt;
    end
    
    if isfield(D{1}, 'allp')
        E.allp{ii} = D{ii}.allp;
    end
    
    
end


function openOrthoViews(hObject, event_data, hD)

figH = get(hObject, 'parent');

if (strcmpi(get(figH, 'selectionType'), 'open'))
    set(figH, 'pointer', 'watch');
    figName = get(figH, 'name');
    tag = 'Orth views';
    gH = icatb_getGraphics(figName, 'normal',tag, 'off');
    hD.currentFigure = gH;
    colormap(hD.cmap);
    axesH(1) = subplot(1, 3, 1);
    axesH(2) = subplot(1, 3, 2);
    axesH(3) = subplot(1, 3, 3);
    drawSlices(axesH, hD);
    set(gH, 'visible', 'on');
    set(gH, 'resize', 'on');
    set(figH, 'pointer', 'arrow');
end

function openUnivStats(hObject, event_data, vars)

figH = get(hObject, 'parent');

if (strcmpi(get(figH, 'selectionType'), 'open'))
    for nV = 1:length(vars)
        if (ishandle(vars{nV}))
            break;
        end
    end
    
    set(figH, 'pointer', 'watch');
    gH = icatb_getGraphics(get(figH, 'name'), 'normal', 'univ stats', 'off');
    colormap(get(figH, 'colormap'));
    sh = subplot(1, 1, 1);
    vars{nV} = sh;
    plotUnivStats(vars{:});
    set(gH, 'visible', 'on');
    set(gH, 'resize', 'on');
    set(figH, 'pointer', 'arrow');
end

function openUnivEffects(hObject, event_data, vars, ch)

figH = get(hObject, 'parent');

if (strcmpi(get(figH, 'selectionType'), 'open'))
    for nV = 1:length(vars)
        if (ishandle(vars{nV}))
            break;
        end
    end
    set(figH, 'pointer', 'watch');
    gH = icatb_getGraphics(get(figH, 'name'), 'normal', 'univ effects', 'off');
    colormap(get(figH, 'colormap'));
    sh = subplot(1, 1, 1);
    vars{nV} = sh;
    ch2 = plot_univariate_effects(vars{:});
    if (exist('ch', 'var'))
        set(ch2, 'YLim', get(ch, 'YLim'));
        set(ch2, 'YTick', get(ch, 'YTick'));
        set(ch2, 'YTickLabel', get(ch, 'YTickLabel'));
    end
    set(gH, 'visible', 'on');
    set(gH, 'resize', 'on');
    set(figH, 'pointer', 'arrow');
end


function [term_no, con_no] = getTermAndConNum(term_name, UNI)
%% Get term and contrast no
%

term_no = [];
con_no = [];
for nT = 1:length(UNI.tests)
    for nCon = 1:length(UNI.stats{nT}.Contrast)
        %tmp = [UNI.tests{nT}, ' ', UNI.stats{nT}.Contrast{nCon}];
        chk = (strcmpi(term_name, UNI.stats{nT}.Contrast{nCon}) || strcmpi(term_name, UNI.tests{nT}));
        if (chk)
            term_no = nT;
            con_no = nCon;
            return;
        end
    end
end

function B = getBetaWeights(UNIStats, con_no)
%% Get betas associated with the term
%

term = UNIStats.Term;
good_inds = [];
for nU = 1:length(UNIStats.Terms)
    if (isequal(UNIStats.Terms{nU}, term))
        good_inds = [good_inds, nU];
    end
end

B = UNIStats.B(good_inds, :);

if (con_no > size(B, 1))
    contrasts = UNIStats.Levels;
    B = B(contrasts(con_no, 1), :) - B(contrasts(con_no, 2), :);
else
    B = B(con_no, :);
end



function [levelNames, levelMeans, titleStr] = plotContextMenu(mancovanInfo, term_name, nF, timeNo)
%% Plot context menu (level means of categorical covariate)
%

levelMeans = {};
levelNames = {};
titleStr = '';

result_files = mancovanInfo.outputFiles(nF).filesInfo.result_files;
if (timeNo > 0)
    load(fullfile(mancovanInfo.outputDir, result_files{1}), 'time', 'fnc_corrs');
    UNI = time.UNI{timeNo};
else
    load(fullfile(mancovanInfo.outputDir, result_files{1}), 'UNI', 'fnc_corrs');
end

[term_no, con_no] = getTermAndConNum(term_name, UNI);

if (~isempty(term_no))
    
    covariateName = UNI.tests{term_no};
    matchIndex = strmatch(lower(covariateName), lower (cellstr(char(mancovanInfo.cov.name))), 'exact');
    
    if (~isempty(matchIndex))
        if (strcmpi(mancovanInfo.cov(matchIndex).type, 'categorical'))
            titleStr = mancovanInfo.cov(matchIndex).name;
            values = mancovanInfo.cov(matchIndex).value;
            levelNames = unique(values);
            if (isfield(mancovanInfo, 'time'))
                if (timeNo > 0)
                    % Load only those subjects and covariates
                    dat = fnc_corrs(mancovanInfo.time.subjects{timeNo}, :);
                else
                    dat = fnc_corrs(mancovanInfo.time.subjects{1}, :) - fnc_corrs(mancovanInfo.time.subjects{2}, :);
                end
                levelMeans = computeLevelMeans(dat, values(mancovanInfo.time.subjects{1}), levelNames);
            else
                levelMeans = computeLevelMeans(fnc_corrs, values, levelNames);
            end
        end
    end
    
end

function levelMeans = computeLevelMeans(fnc_corrs, values, levelNames)
%% Compute level means
%

levelMeans = cell(1, length(levelNames));
for nL = 1:length(levelNames)
    chkInd = strcmpi(levelNames{nL}, values);
    tmp = squeeze(mean(fnc_corrs(chkInd, :)));
    levelMeans{nL} = tmp;
end


function plotLevelMeans(hObject, event_data, mancovanInfo)
%% Plot level means
%

load icatb_colors coldhot;
ud = get(hObject, 'userdata');

levelNames = ud.levelNames;
levelMeans = ud.levelMeans;
featureName = ud.featureName;
titleStr = ud.title;

if (strcmpi(featureName, 'fnc correlations'))
    
    network_values = zeros(1, length(mancovanInfo.userInput.comp));
    for nV = 1:length(network_values)
        network_values(nV) = length(mancovanInfo.userInput.comp(nV).value);
    end
    network_names =  cellstr(char(mancovanInfo.userInput.comp.name));
    
    
    dd = [levelMeans{:}];
    CLIM = icatb_z_to_r(max(abs(dd(:))));
    
    gH = icatb_getGraphics(['Level Means (', titleStr, ')'], 'graphics',  ['level_means_', titleStr], 'on');
    
    for nL = 1:length(levelNames)
        sh = subplot(length(levelNames), 1, nL);
        if (~isempty(levelMeans{nL}))
            M = icatb_vec2mat(icatb_z_to_r(levelMeans{nL}), 1);
            %CLIM = max(abs(M(:)));
            fig_title = 'Features (FNC Correlations)';
            
            if (length(network_names) == 1)
                gH = icatb_plot_matrix(M, cellstr(num2str(mancovanInfo.comps(:))), cellstr(num2str(mancovanInfo.comps(:))), 'title', ['Corr (', levelNames{nL}, ')'], 'cmap', ...
                    coldhot, 'clim', [-CLIM, CLIM], 'ylabel', ['Corr (', levelNames{nL}, ')'], 'xlabel', 'Components', 'axesh', sh);
            else
                % gH = icatb_getGraphics(fig_title, 'graphics',  'FNC Correlations', 'on');
                set(gH, 'resize', 'on');
                %axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
                icatb_plot_FNC(M, [-CLIM, CLIM], cellstr(num2str(mancovanInfo.comps(:))), (1:length(mancovanInfo.comps)), gH, ['Corr (', levelNames{nL}, ')'], sh, network_values, network_names);
                colormap(coldhot);
            end
        end
        clear M;
    end
end

