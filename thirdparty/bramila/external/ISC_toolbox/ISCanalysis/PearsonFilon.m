function PearsonFilon(Params,nrSession,freqBlock)
% This function calculates sum ZPF statistics.

showTime(1);

Pub = Params.PublicParams;
Priv = Params.PrivateParams;

if ~Pub.corOn
    disp('Correlation criterion not selected, sum Z-Pearson-Filon map calculation skipped...')
    return
end

if ~Pub.calcStandard
    disp('Standard inter-subject synchronization analysis not specified, sum Z-Pearson-Filon map calculation skipped...')
    return
end

if Pub.nrFreqBands == 0
    disp('No frequency subbands selected, sum Z-Pearson Filon map calculation skipped...')
    return
end

% total number of frequency comparisons and subject pairs:
nrFreqComp = ((Priv.maxScale+2)^2-(Priv.maxScale+2))/2;
nrSubjPairs = ((Priv.nrSubjects)^2-(Priv.nrSubjects))/2;

% load pointers:
load([Pub.dataDestination 'memMaps'])
% fMRI-data (source data):
mMapFilt = memMaps.(Priv.filtMapName).([Priv.prefixSession num2str(nrSession)]);
mMapOrig = memMaps.(Priv.origMapName).([Priv.prefixSession num2str(nrSession)]);
% pearson-filon-data (destination data):
%mMapResultWhole = memMaps.(Priv.PFMapName).whole.([Priv.prefixSession num2str(nrSession)]).cor;
mMapmatResultWhole = memMaps.(Priv.PFmatMapName).whole.([Priv.prefixSession num2str(nrSession)]).cor;
clear memMaps

for fr = 1:length(freqBlock)
%    mMapmatResultWhole.([Priv.prefixFreqComp num2str(freqBlock(fr))]).Writable = true
    if ( mMapmatResultWhole.([Priv.prefixFreqComp num2str(freqBlock(fr))]).Writable == false ) %...
       % || ( mMapResultWhole.([Priv.prefixFreqComp num2str(freqBlock(fr))]).Writable == false ) 
        disp('Sum ZPF statistic computed already...')
        return
    end
end

% create look-up table for frequency band comparisons:
freqBandTable = zeros(3,length(freqBlock));
iter = 1;
it = 1;
for rr = 1:Priv.maxScale+2
    for cc = 1:Priv.maxScale+2
        if cc > rr
            IND = find(freqBlock == iter);
            if ~isempty(IND)
                freqBandTable(1,it) = rr;
                freqBandTable(2,it) = cc;
                freqBandTable(3,it) = iter;
                it = it + 1;
            end
            iter = iter + 1;
        end
    end
end
freqBandTable


VAL = zeros(1,nrSubjPairs);
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
% INDS = [];
% for hh = 1:Priv.nrSubjects-1
%   INDS = [INDS (1+hh*Priv.nrSubjects):(...
%   1+hh*Priv.nrSubjects+hh-1)];
% end

% significance thresholds for pairwise ZPF-statistic:
Th(1) = 1.96;
Th(2) = 2.807;
Th(3) = 3.2905;

%Z1 = zeros([Priv.dataSize(nrSession,1),Priv.dataSize(nrSession,2),Priv.dataSize(nrSession,3),length(freqBlock),7]);
Z2 = zeros([Priv.dataSize(nrSession,1),Priv.dataSize(nrSession,2),Priv.dataSize(nrSession,3),length(freqBlock),nrSubjPairs]);

iter = 0;
for xx = 1:Priv.dataSize(nrSession,1)
    disp(['x: ' num2str(xx) '/' num2str(Priv.dataSize(nrSession,1))])
    % init result matrix for across session data
    % process only non-zero x-slices:
    if sum(sum(squeeze(bmask(xx,:,:)))) > 0
        % init source time-series matrix:
        cDat = zeros([Priv.dataSize(nrSession,4),Priv.maxScale+2, Priv.nrSubjects,Priv.dataSize(nrSession,2),Priv.dataSize(nrSession,3)]);
        % read multi-band time-series of the subjects:
        for f = 0:Priv.maxScale+1
            if ~isempty(find(freqBandTable(1:2,:)==(f+1)))
                %   disp(['load band: ' num2str(f)])
                for k = 1:Priv.nrSubjects
                    if f == 0 % load full-band data
                        cDat(:,k,f+1,:,:) = mMapOrig.([Priv.prefixSubject num2str(k)]).Data(xx).tyz;
                    else % load sub-band data
                        cDat(:,k,f+1,:,:) = mMapFilt.([Priv.prefixSubjectFilt ...
                            num2str(k)]).([Priv.prefixFreqBand num2str(f)]).Data(xx).tyz;
                    end
                end
            end
        end
        for yy = 1:Priv.dataSize(nrSession,2)
            ziter = 0;
            for zz = 1:Priv.dataSize(nrSession,3)
                if bmask(xx,yy,zz)
                    ziter = ziter + 1;
                    ts = cDat(:,:,:,yy,zz);
                    
                    wfr = 1;
                    %  for wfr = 1:Priv.nrTimeIntervals(nrSession)
                    % calculate across-whole-session similarity values as a special case:
                    if wfr == 1
                        N = Priv.dataSize(nrSession,4);
                        % select time-series from different frequency bands:
                        for freqbInd = 1:length(freqBlock)
                            g = freqBandTable(1,freqbInd);
                            h = freqBandTable(2,freqbInd);
                            ts_v1 = ts(:,:,g);
                            ts_v2 = ts(:,:,h);
                            % select subject pairs:
                            subjpairInd = 0;
                            for m = 1:Priv.nrSubjects
                                for n = 1:Priv.nrSubjects
                                    if n > m
                                        % Compute Pearson-Filon statistic
                                        % for each subject pair:                                  
                                        subjpairInd = subjpairInd + 1;
                                        T = [ts_v1(:,m) ts_v1(:,n) ts_v2(:,m) ts_v2(:,n)];
                                        % calculate 4 by 4 correlation matrix (ref: Matlab corrcoef.m):
                                        [n1,m1] = size(T);
                                        suT = sum(T)/n1;
                                        xc = T - ones(n1,1)*suT;  % remove mean
                                %        xc2 = T - repmat(suT,n1,1);  % Remove mean
                                        c1 = (xc' * xc) / (n1-1); % calculate inner products
                                        d1 = sqrt(diag(c1)); % sqrt first to avoid under/overflow
                                        dd = d1*d1';
                                        dd(1:m1+1:end) = diag(c1); % remove roundoff on diag
                                        r1 = c1 ./ dd;
                                        % calculate Fisher's z-transform:
                                        z12 = 0.5*(log((1+r1(1,2))./(1-r1(1,2))));
                                        z34 = 0.5*(log((1+r1(3,4))./(1-r1(3,4))));
                                        % calculate Z-statistic (Pearson-Filon based on Fisher's z-transformation):
                                        Q = 0.5*( ...
                                            ( r1(1,3)-r1(3,2)*r1(1,2) ) * ( r1(4,2)-r1(3,2)*r1(3,4) ) + ...
                                            ( r1(1,4)-r1(1,3)*r1(3,4) ) * ( r1(2,3)-r1(1,3)*r1(1,2) ) + ...
                                            ( r1(1,3)-r1(1,4)*r1(3,4) ) * ( r1(2,4)-r1(1,4)*r1(1,2) ) + ...
                                            ( r1(1,4)-r1(1,2)*r1(2,4) ) * ( r1(2,3)-r1(2,4)*r1(3,4) ) );
                                        VAL(subjpairInd) = ( (z12-z34)*sqrt((N-3)/2) ) / ...
                                            sqrt( 1 - Q / ( (1-r1(1,2)^2)*(1-r1(3,4)^2) ) );
                                    end
                                end
                            end
                            Z2(xx,yy,zz,freqbInd,:) = VAL;
%                            Z1(xx,yy,zz,freqbInd,1) = sum(VAL >= Th(1));
%                            Z1(xx,yy,zz,freqbInd,2) = sum(VAL <= -Th(1));
%                            Z1(xx,yy,zz,freqbInd,3) = sum(VAL >= Th(2));
%                            Z1(xx,yy,zz,freqbInd,4) = sum(VAL <= -Th(2));
%                            Z1(xx,yy,zz,freqbInd,5) = sum(VAL >= Th(3));
%                            Z1(xx,yy,zz,freqbInd,6) = sum(VAL <= -Th(3));
%                            Z1(xx,yy,zz,freqbInd,7) = sum(VAL);
                        end
                    end % end across-session.
                    
                    % NOTE! PEARSON-FILON TEST FOR TIME WINDOW ISC 
                    % NOT CURRENTLY IN USE DUE TO HIGH COMPUTATION TIME
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % calculate similarity values for time-frames:
                    %            freqbInd = 0;
                    %            N = Pub.windowSize;
                    %            ts_win = ts(Priv.startInds{nrSession}...
                    %            (wfr):Priv.startInds{nrSession}...
                    %            (wfr)+Pub.windowSize-1,:,:);
                    %
                    %            for g = 1:Priv.maxScale + 2
                    %              for h = 1:Priv.maxScale + 2
                    %                if h > g
                    %                  freqbInd = freqbInd + 1;
                    %                  ts_v1 = squeeze(ts_win(:,g,:));
                    %                  ts_v2 = squeeze(ts_win(:,h,:));
                    %                  % select subject pairs:
                    %                  subjpairInd = 0;
                    %                  for m = 1:Priv.nrSubjects
                    %                    for n = 1:Priv.nrSubjects
                    %                      if n > m
                    %                        subjpairInd = subjpairInd + 1;
                    %                        T = [ts_v1(:,m) ts_v1(:,n) ts_v2(:,m) ts_v2(:,n)];
                    %                        % calculate 4 by 4 correlation matrix (ref: Matlab corrcoef.m):
                    %                        [n1,m1] = size(T);
                    %                        xc = T - repmat(sum(T)/n1,n1,1);  % Remove mean
                    %                        c1 = (xc' * xc) / (n1-1); % calculate inner products
                    %                        d1 = sqrt(diag(c1)); % sqrt first to avoid under/overflow
                    %                        dd = d1*d1';
                    %                        dd(1:m1+1:end) = diag(c1); % remove roundoff on diag
                    %                        r1 = c1 ./ dd;
                    %
                    %                        % calculate Fisher's z-transformation:
                    %                        z12 = 0.5*(log((1+r1(1,2))./(1-r1(1,2))));
                    %                        z34 = 0.5*(log((1+r1(3,4))./(1-r1(3,4))));
                    %                        % calculate Z-statistic (Pearson-Filon based on Fisher's z transformation):
                    %                        Q = 0.5*( ...
                    %                        ( r1(1,3)-r1(3,2)*r1(1,2) ) * ( r1(4,2)-r1(3,2)*r1(3,4) ) + ...
                    %                        ( r1(1,4)-r1(1,3)*r1(3,4) ) * ( r1(2,3)-r1(1,3)*r1(1,2) ) + ...
                    %                        ( r1(1,3)-r1(1,4)*r1(3,4) ) * ( r1(2,4)-r1(1,4)*r1(1,2) ) + ...
                    %                        ( r1(1,4)-r1(1,2)*r1(2,4) ) * ( r1(2,3)-r1(2,4)*r1(3,4) ) );
                    %                        %tic
                    %                        VAL(subjpairInd) = ( (z12-z34)*sqrt((N-3)/2) ) / ( 1-Q/(1-r1(1,2)^2)*(1-r1(3,4)^2) );
                    %
                    %                        %Q=toc;
                    %                        if ( wfr == 1 && n == 2 && m == 1 && h == 2 && g == 1 && mod(zz,10)==0 )
                    %                          disp(Q)
                    %                          VAL
                    %                          size(Zw)
                    %                          size(Z)
                    %                          freqbInd
                    %                          subjpairInd
                    %                          wfr
                    %                          yy
                    %                          zz
                    %                        end
                    %                       disp(['Frame' num2str(wfr) 'fb' num2str(g) num2str(h) 'subj' num2str(m) num2str(n) ' ' num2str(QQ) 's'])
                    %                       disp(['y: ' num2str(yy) ' z: ' num2str(zz) ' Freq: ' num2str(freqbInd) ' Subj: ' num2str(subjpairInd)'])
                    %                       Zw(yy,zz,freqbInd,subjpairInd,wfr)
                    %                       assignin('base',['TSfreq' num2str(g) num2str(h) 'subj' num2str(m) num2str(n)],T)
                    %                       %        assignin('base',['ZPFfreq' num2str(g) num2str(h) 'subj' num2str(m) num2str(n)],Z(yy,zz,freqbInd,subjpairInd))
                    %                     end
                    %                   end
                    %                 end
                    %                 Zw(zz,freqbInd,1,wfr) = sum(VAL >= 1.96);
                    %                 Zw(zz,freqbInd,2,wfr) = sum(VAL <= -1.96); % 1.645
                    %               end
                    %             end
                    %           end                              
                end
            end % end z
            % ...
            % ...
            % ...
            % ...
        end % end y 
    end
end % end x

for fr = 1:length(freqBlock)
    %           mMapResultWhole.([Priv.prefixFreqComp num2str(freqBlock(fr))...
    %               ]).Data.xyzc(:,:,:,:) = single(squeeze(Z1(:,:,:,fr,:)));
    mMapmatResultWhole.([Priv.prefixFreqComp num2str(freqBlock(fr))...
        ]).Data.xyzc = single(squeeze(Z2(:,:,:,fr,:)));
    %         for wfr = 1:Priv.nrTimeIntervals(nrSession)
    %           memMaps.(Priv.PFMapName).win.([Priv.prefixSession...
    %           num2str(nrSession)]).cor.([Priv.prefixFreqComp ...
    %           num2str(fr)]).Data(wfr).xyzc(xBlock(xx),yy,:,:) = ...
    %           single(squeeze(Zw(:,fr,:,wfr)));
    %         end
end

% Set data non-writable:
load([Pub.dataDestination 'memMaps'])
for fr = 1:length(freqBlock)
%    memMaps.(Priv.PFMapName).whole.([Priv.prefixSession num2str(nrSession)]).cor.([Priv.prefixFreqComp num2str(freqBlock(fr))]).Writable = false;
    memMaps.(Priv.PFmatMapName).whole.([Priv.prefixSession num2str(nrSession)]).cor.([Priv.prefixFreqComp num2str(freqBlock(fr))]).Writable = false;
end
save([Pub.dataDestination 'memMaps'],'memMaps')
clear memMaps

showTime(0);
