function PearsonFilonAcrossSessions(Params,bandNr,sessBlock)
% This function calculates sum ZPF statistics.

% See also: RUNANALYSIS, INITPARAMS

% Last updated: 21.6.2011 by Jukka-Pekka Kauppi
% Tampere University of Technology
% Department of Signal Processing
% e-mail: jukka-pekka.kauppi@tut.fi

showTime(1);


Pub = Params.PublicParams;
Priv = Params.PrivateParams;

if Pub.sessionCompOn == 0
    disp('Sum ZPF map computation across sessions not selected...')
    return
end


if ~Pub.corOn
    disp('Correlation coefficient needs to be selected in order to compute sum ZPF maps...')
    return
end

if ~Pub.calcStandard
    disp('Correlation coefficient needs to be selected in order to compute sum ZPF maps...')
    return
end

for xx = 1:Priv.nrSessions
    za(xx) = Priv.dataSize(xx,4);
end
if length(unique(za)) ~= 1
    disp('Each session must have same number of data points in order to compute sum ZPF maps across sessions...')
    return
end

if Priv.nrSessions < 2
    disp('Only one session specified, sumZPF across sessions cannot be computed...')
    return
end

% total number of subject pairs:
nrSubjPairs = ((Priv.nrSubjects)^2-(Priv.nrSubjects))/2;
% total number of session pairs:
nrSesComp = ((Priv.nrSessions)^2-(Priv.nrSessions))/2;

if nargin < 3
    sessBlock = 1:nrSesComp;
end

% load pointers:
load([Pub.dataDestination 'memMaps'])

for h = 1:Priv.nrSessions
    % fMRI-data (source data):
    if bandNr > 0
        mMap.(['sess' num2str(h)]) = memMaps.(Priv.filtMapName).([Priv.prefixSession num2str(h)]);
    else
        mMap.(['sess' num2str(h)]) = memMaps.(Priv.origMapName).([Priv.prefixSession num2str(h)]);
    end
end

% pearson-filon-data (destination data):
%mMapResultWhole = memMaps.(Priv.PFMapSessionName).whole.([Priv.prefixFreqBand num2str(bandNr)]).cor;
mMapmatResultWhole = memMaps.(Priv.PFmatMapSessionName).whole.([Priv.prefixFreqBand num2str(bandNr)]).cor;
clear memMaps

for fr = 1:length(sessBlock)
%    if ( mMapResultWhole.([Priv.prefixSessComp num2str(sessBlock(fr))]).Writable == false ) || ...
%            ( mMapmatResultWhole.([Priv.prefixSessComp num2str(sessBlock(fr))]).Writable == false )
  if ( mMapmatResultWhole.([Priv.prefixSessComp num2str(sessBlock(fr))]).Writable == false )        
      disp('Session comparison map already computed...')
        return
    end
end


% create look-up table for session comparisons:
sessBandTable = zeros(3,length(sessBlock));
iter = 1;
it = 1;
for rr = 1:Priv.nrSessions
    for cc = 1:Priv.nrSessions
        if cc > rr
            IND = find(sessBlock == iter);
            if ~isempty(IND)
                sessBandTable(1,it) = rr;
                sessBandTable(2,it) = cc;
                sessBandTable(3,it) = iter;
                it = it + 1;
            end
            iter = iter + 1;
        end
    end
end
sessBandTable


VAL = zeros(1,nrSubjPairs);

% if Pub.useTemplate
%     bmask = load_nii(Priv.brainMask);
%     bmask = single(bmask.img);
% else
%     bmask = load(Priv.brainMask);
%     fiel = fields(bmask);
%     bmask = bmask.(fiel{1});
%     bmask = single(bmask);
% end
% bmask = logical(bmask);

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

% significance thresholds for pairwise ZPF-statistic:
Th(1) = 1.96;
Th(2) = 2.807;
Th(3) = 3.2905;


% In each spatial x-location, we save time-series to the variable
% cDat(a,b,c,d,e), where:
% a = length of the time-series
% b = spatial y-dimension
% c = spatial z-dimension
% d = number of sessions
% e = number of subjects.
% For instance, if there are 12 subjects, 3 sessions, and 91x109x91x244 
% fMRI data, size(cData) = [244 109 91 3 12].
iter = 0;
for xx = 1:Priv.dataSize(1,1)
    disp(['x: ' num2str(xx) '/' num2str(Priv.dataSize(1,1))])
    % init result matrix for across session data
    % process only non-zero x-slices:
    if sum(sum(squeeze(bmask(xx,:,:)))) > 0
        % init source time-series matrix:
        cDat = zeros([Priv.dataSize(1,[4 2 3]),Priv.nrSessions, Priv.nrSubjects]);
        % read multi-band time-series of the subjects:
        for f = 1:Priv.nrSessions
            if ~isempty(find(sessBandTable(1:2,:)==f))
                %   disp(['load band: ' num2str(f)])
                for k = 1:Priv.nrSubjects
                    if bandNr == 0
                        cDat(:,:,:,f,k) = mMap.(['sess' num2str(f)]).([Priv.prefixSubject num2str(k)]).Data(xx).tyz; % time*y*z*sessions*subjects
                    else
                        cDat(:,:,:,f,k) = mMap.(['sess' num2str(f)]).([Priv.prefixSubjectFilt num2str(k)]).([Priv.prefixFreqBand num2str(bandNr)]).Data(xx).tyz;
                    end
                end
            end
        end

        %%%%%%%

        for yy = 1:Priv.dataSize(1,2)
  %          Z1 = zeros([Priv.dataSize(1,3),length(sessBlock),7]);
            Z2 = zeros([Priv.dataSize(1,3),length(sessBlock),nrSubjPairs]);
            ziter = 0;
            for zz = 1:Priv.dataSize(1,3)
                if bmask(xx,yy,zz)
                    ziter = ziter + 1;
                    ts = squeeze(cDat(:,yy,zz,:,:)); % time*sessions*subjects
                    N = Priv.dataSize(1,4); % number of time points
                    % select time-series from different sessions:
                    for sessInd = 1:length(sessBlock)
                        % take time-series from session g and h:
                        g = sessBandTable(1,sessInd);
                        h = sessBandTable(2,sessInd);
                        ts_v1 = squeeze(ts(:,g,:)); % time*subjects, session g
                        ts_v2 = squeeze(ts(:,h,:)); % time*subjects, session h

                        % select subject pairs:
                        subjpairInd = 0;
                        for m = 1:Priv.nrSubjects
                            for n = 1:Priv.nrSubjects
                                if n > m
                                    subjpairInd = subjpairInd + 1;
                                    T = [ts_v1(:,m) ts_v1(:,n) ts_v2(:,m) ts_v2(:,n)];
 
                                    % calculate 4 by 4 correlation matrix (ref: Matlab corrcoef.m):
                                    [n1,m1] = size(T);
                                    xc = T - repmat(sum(T)/n1,n1,1);  % Remove mean
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
                                    %                         disp(['y: ' num2str(yy) ' z: ' num2str(zz) ' Sess: ' num2str(sessInd) ' Subj: ' num2str(subjpairInd)'])
                                    VAL(subjpairInd) = ( (z12-z34)*sqrt((N-3)/2) ) / ...
                                        sqrt( 1 - Q / ( (1-r1(1,2)^2)*(1-r1(3,4)^2) ) );
                                    %                       toc
                                end
                            end
                        end
                        Z2(zz,sessInd,:) = VAL;
 %                       Z1(zz,sessInd,1) = sum(VAL >= Th(1));
 %                       Z1(zz,sessInd,2) = sum(VAL <= -Th(1));
 %                       Z1(zz,sessInd,3) = sum(VAL >= Th(2));
 %                       Z1(zz,sessInd,4) = sum(VAL <= -Th(2));
 %                       Z1(zz,sessInd,5) = sum(VAL >= Th(3));
 %                       Z1(zz,sessInd,6) = sum(VAL <= -Th(3));
 %                       Z1(zz,sessInd,7) = sum(VAL);
                    end

                end
            end % end z
            for fr = 1:length(sessBlock)
%                mMapResultWhole.([Priv.prefixSessComp num2str(sessBlock(fr))...
%                    ]).Data.xyzc(xx,yy,:,:) = single(squeeze(Z1(:,fr,:)));

                mMapmatResultWhole.([Priv.prefixSessComp num2str(sessBlock(fr))...
                    ]).Data.xyzc(xx,yy,:,:) = single(squeeze(Z2(:,fr,:)));

            end
        end

    end % end y
end % end x

% Set data non-writable:
load([Pub.dataDestination 'memMaps'])
for fr = 1:length(sessBlock)
    memMaps.(Priv.PFmatMapSessionName).whole.([Priv.prefixFreqBand num2str(bandNr)]).cor.([Priv.prefixSessComp num2str(sessBlock(fr))]).Writable = false;
end
save([Pub.dataDestination 'memMaps'],'memMaps')
clear memMaps

showTime(0);
