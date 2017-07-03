function permutationTest(Params,nrSession,permSetIdx,winOn)

% This function generates null-distribution using random shuffling of the
% time-series before calculating the test statistic (mean correlation coefficient).
% Random (single-point) shuffling is used to preserve autocorrelation structure of  
% the time-series. Distribution is generated based on randomly picked voxels across 
% the whole brain from random frequency band (e.g., 10e8 voxels are randomly picked).
% For time-windowed data, also time-frames are selected randomly.
%
% Note: It is practical to run simulation in several parallel batches (e.g. in 10 sets of 
% 10e7 realizations) to make computation faster. To change parameters for the 
% number of batches and the number of realizations, change the following fields in the 
% parameter-struct before running this function:
% 1) PrivateParams.nrPermutationSets (number of batches)
% 2) PrivateParams.nrPermutations (number of realizations)
%
% Result of each batch is saved in separate mat-file (PrivateParams.statsDestination). From
% these files, thresholds can be calculated using calculateThreholds.m -function.
%
% Function inputs:
% Params - struct containing all necessary parameters
% nrSession - session index
% permSetIdx -  index of the bootstrap distribution batch
% winOn - 1 performs calculations for windowed data, 0 for across-session
% data
%
% See also:
% CALCULATETHRESHOLDS
% COMBINENULLVALS
% INITPARAMS
% RUNANALYSIS


% Modified 30.10.2013
% Jukka-Pekka Kauppi
% University of Helsinki
% Department of Computer Science
% jukka-pekka.kauppi@helsinki.fi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

showTime(1);


if nargin < 3
  error('At least three inputs must be given!!')

end
if nargin == 3
    winOn = 0;
end

Pub = Params.PublicParams;
Priv = Params.PrivateParams;

if ~Pub.calcStandard
    disp('Basic ISC analysis not selected, no resampling distribution generated...')
    return
end

if ~Pub.winOn && winOn
%  disp('Time-window ISC not selected...')
  return
end

if winOn
    if exist([Priv.statsDestination 'vals' num2str(permSetIdx) Priv.prefixSession num2str(nrSession) 'win.mat'],'file') == 2
        disp(['Resampling distribution for time-window ISC already exists...'])
        return
    end
else
    fn = [Priv.statsDestination 'vals' num2str(permSetIdx) Priv.prefixSession num2str(nrSession) '.mat'];
    if exist(fn,'file') == 2
        ww = whos('-file',fn);
        flag = [0 0];
        if Pub.calcPhase
            flag(2) = 0;
        else
            flag(2) = 1;
        end
        for m = 1:length(ww)
           if strcmp(ww(m).name,'val')
               flag(1) = 1;
           end
           if strcmp(ww(m).name,'ph')
               flag(2) = 1;
           end
        end
        if sum(flag) == 2
            disp(['Resampling distribution for basic ISC already exists...'])        
            if Pub.calcPhase
                disp(['Resampling distribution for phase synchronization maps already exists...'])                        
            end
            return
        end
    end
end

load([Pub.dataDestination 'memMaps'])
cDat = zeros([Priv.dataSize(nrSession,[4 2 3]), Priv.nrSubjects]);
%INDS = [];

if strcmp(Pub.fileFormat,'nii')
    bmask = load_nii(Priv.brainMask);
    bmask = single(bmask.img);
elseif strcmp(Pub.fileFormat,'mat')
    bmask = load(Priv.brainMask);
    fiel = fields(bmask);
    bmask = bmask.(fiel{1});
    bmask = single(bmask);
else
    error('Mask must be mat- or nii-file!')
end
bmask = logical(bmask);

INDS = find(triu(ones(Priv.nrSubjects,Priv.nrSubjects),1));

% add extra iterations to cope with NaNs:
nanCount = 0;
nrPermOrig = Priv.nrPermutations;
Priv.nrPermutations = round(2*Priv.nrPermutations);

iter = 0;
val = zeros(1,Priv.nrPermutations);
valFish = val;
if Pub.calcPhase && ~Pub.winOn
    ph = val*NaN;
end
brainVoxels = find(bmask);
nrVoxels = length(brainVoxels);
nrBands = Pub.nrFreqBands;
nrTimePoints = Priv.dataSize(nrSession,4);
nrTimePointsWin = Pub.windowSize;
nrTimeWindows = Priv.nrTimeIntervals(nrSession);

rand('state',sum(100*clock));

randInds = ceil(nrVoxels*rand(Priv.nrPermutations,1));
randVoxels = brainVoxels(randInds);
[x,y,z] = ind2sub(Priv.dataSize(nrSession,1:3),randVoxels);
[x inds] = sort(x);
y = y(inds);
z = z(inds);

randBands = ceil((nrBands+1)*rand(Priv.nrPermutations,1))-1;
randTimeWindows = ceil(nrTimeWindows*rand(Priv.nrPermutations,1));
randTimePoints = ceil(nrTimePoints*rand(Priv.nrPermutations,Priv.nrSubjects));
randTimePointsWin = ceil(nrTimePointsWin*rand(Priv.nrPermutations,Priv.nrSubjects));

for iter = 1:Priv.nrPermutations
    
    %     if mod(iter,ceil(nrPermOrig/100) ) == 0
    %         disp(['Permutation: ' num2str(iter)])
    %     end
    % get mapped source data of the subjects:
    if ( iter == 1 ) || ( x(iter) ~= x(iter-1) )
        % tic
        for k = 1:Priv.nrSubjects
            if randBands(iter) == 0 % load full-band data
                cDat(:,:,:,k) = memMaps.(Priv.origMapName).([Priv.prefixSession ...
                    num2str(nrSession)]).([Priv.prefixSubject ...
                    num2str(k)]).Data(x(iter)).tyz(:,:,:);
            else % load sub-band data
                cDat(:,:,:,k) = memMaps.(Priv.filtMapName).([Priv.prefixSession ...
                    num2str(nrSession)]).([Priv.prefixSubjectFilt ...
                    num2str(k)]).([Priv.prefixFreqBand num2str(randBands(iter))]).Data(x(iter)).tyz(:,:,:);
            end
        end
        %          toc
    end
    
    % obtain each subject's time series:
    ts = squeeze(cDat(:,y(iter),z(iter),:));
    tss = zeros(size(ts));
    
    if ~winOn
        % shuffle time-series:
        for k = 1:size(ts,2)
            tt = ts(:,k);
            tss(:,k) = tt([randTimePoints(iter,k)+1:end 1:randTimePoints(iter,k)]);
        end
        %   tic
        % calculate across-whole-session similarity values:
        %N = Priv.dataSize(nrSession,4);
        
        [n1,m1] = size(tss);
        xc = tss - repmat(sum(tss)/n1,n1,1);  % Remove mean
        c1 = (xc' * xc) / (n1-1); % calculate inner products
        d1 = sqrt(diag(c1)); % sqrt first to avoid under/overflow
        dd = d1*d1';
        dd(1:m1+1:end) = diag(c1); % remove roundoff on diag
        r1 = c1 ./ dd;
        r1 = r1(INDS);
        val(iter) = single(mean(r1));
        r1 = 0.5*(log(1+r1)./(1-r1));
        valFish(iter) = single(mean(r1));
        if isnan(val(iter))
            nanCount = nanCount + 1;
        end
        
        if Pub.calcPhase
            if randBands(iter) ~= 0
                subjPairs = Priv.nrSubjects*(Priv.nrSubjects-1)/2;
                Ts = tss';
                Ts = Ts - mean(Ts,2)*ones(size(Ts,2),1)'; % remove mean
                ite = 1;
                for kk = 1:Priv.nrSubjects
                    y1 = hilbert(Ts(kk,:));
                    for hh = 1:Priv.nrSubjects
                        if hh > kk
                            y2 = hilbert(Ts(hh,:));
                            Rm(ite,:) = angle(y1.*conj(y2));
                            ite = ite + 1;
                        end
                    end
                end
                if Priv.nrSubjects == 2
                    phh = abs(Rm);
                else
                    phh = sum(abs(Rm))/subjPairs;
                end
                phv = single(1 - phh/pi);
                rn = randperm(length(phv));
                ph(iter) = phv(rn(1));
                Fi = find(isnan(ph(1:iter)));
                if ~isempty(Fi)
                    if length(Fi) >= length(rn(2:end))
                        ph(Fi(1:(length(rn(2:end))))) = phv(rn(2:end));
                    else
                        ph(Fi) = phv(rn(2:(1+length(Fi))));
                    end
                end                
            end
        end
        
        %toc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        % calculate similarity values in randomly picked time-frame:
        %N = Pub.windowSize;
        
        ts_win = ts(Priv.startInds{nrSession}...
            (randTimeWindows(iter)):Priv.startInds{nrSession}...
            (randTimeWindows(iter))+Pub.windowSize-1,:);
        
        tss_win = zeros(size(ts_win));
        % shuffle time-series:
        for k = 1:size(ts_win,2)
            tt = ts_win(:,k);
            tss_win(:,k) = tt([randTimePointsWin(iter,k)+1:end 1:randTimePointsWin(iter,k)]);
        end
        
        %  tic
        % Mean of pairwise correlation:
        [n1,m1] = size(tss_win);
        xc = tss_win - repmat(sum(tss_win)/n1,n1,1);  % Remove mean
        c1 = (xc' * xc) / (n1-1); % calculate inner products
        d1 = sqrt(diag(c1)); % sqrt first to avoid under/overflow
        dd = d1*d1';
        dd(1:m1+1:end) = diag(c1); % remove roundoff on diag
        r1 = c1 ./ dd;
        r1 = r1(INDS);
        val(iter) = single(mean(r1));
        r1 = 0.5*(log(1+r1)./(1-r1));
        valFish(iter) = single(mean(r1));
        if isnan(val(iter))
            nanCount = nanCount + 1;
        end
        
        %      if mod(iter,ceil(Priv.nrPermutations/100)) == 0
        %        disp(['Permutation: ' num2str(iter) '/' num2str(nrPermOrig)])
        %        val(iter)
        %        valFish(iter)
        %      end
        
    end
    
    
    if mod(iter,ceil(nrPermOrig/100) ) == 0
        disp(['Realization number: ' num2str(iter) '/' num2str(nrPermOrig)])
        disp(['Resampled ISC value: ' num2str(val(iter))])
    end
    
    %   toc
    % break if the number of non-nan values equals original number of permutations:
    if ( iter - nanCount ) == nrPermOrig
        break
    end
end
% save only non-nan values:
val = val(1:iter);
nn = ~isnan(val);
val = val(nn);
valFish = valFish(1:iter);
nnFish = ~isnan(valFish);
valFish = valFish(nnFish);
if Pub.calcPhase && ~Pub.winOn
    ph = ph(1:iter);
    nnPh = ~isnan(ph);
    ph = ph(nnPh);
end


% x = x(1:iter);
% x = x(nn);
% y = y(1:iter);
% y = y(nn);
% z = z(1:iter);
% z = z(nn);
% randBands = randBands(1:iter);
% randBands = randBands(nn);
% randTimeWindows = randTimeWindows(1:iter);
% randTimeWindows = randTimeWindows(nn);
% randTimePoints = randTimePoints(1:iter);
% randTimePoints = randTimePoints(nn);
% randTimePointsWin = randTimePointsWin(1:iter);
% randTimePointsWin = randTimePointsWin(nn);

%save([Pub.dataDestination 'randInds' num2str(permSetIdx) Priv.prefixSession num2str(nrSession) 'win' num2str(winOn)],'x','y','z','randBands','randTimeWindows','randTimePoints','randTimePointsWin');
if winOn
    save([Priv.statsDestination 'vals' num2str(permSetIdx) Priv.prefixSession num2str(nrSession) 'win'],'val','valFish')
else
    if Pub.calcPhase
        save([Priv.statsDestination 'vals' num2str(permSetIdx) Priv.prefixSession num2str(nrSession)],'val','valFish','ph')
    else        
        save([Priv.statsDestination 'vals' num2str(permSetIdx) Priv.prefixSession num2str(nrSession)],'val','valFish')
    end    
end

showTime(0);
