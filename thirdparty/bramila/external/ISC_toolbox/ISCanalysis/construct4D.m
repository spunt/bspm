function B = construct4D(Params,nrSession,nrBand,corMat,varargin)
%
%
% Return 4D data matrix of either single subject's time-series or
% inter-subject correlations.
%
% inputs:
% Params - struct of analysis parameters
% nrSession - session index
% nrBand - frequency band index
% corMat - type of data to load: 0 = time-series data, 1 = correlation
% matrices, 2 = phase maps
%
% Additional instructions depending on which data set you want to load.
%
% Loading time-series data:
%
% if you specify corMat = 0, then you must specify the following additional parameter:
% nrSubject - subject index
%
% Loading correlation matrices:
%
% if you specify corMat = 1, then you must specify the following additional parameter:
% winOn - load across session/time-window data: 0 = across session, 1 = windowed
%
% if you specify winOn = 1, then you must specify one additional parameter:
% timeInt - time-interval index
%
% Loading phase synchrony data:
%
% if you specify corMat = 2, no additional parameters are required
%
% Example usages:
% B = construct4D(Params,1,0,0,11); % session 1, original band, time-series data, subject 11
% B = construct4D(Params,1,4,1,0); % session 1, subband 4, correlation matrices, across-session
% B = construct4D(Params,1,0,1,1,23); % session 1, original band, correlation matrices, windowed, time-interval 23
% B = construct4D(Params,1,5,2); % session 1, subband 5, phase synchrony data
%
%
% Jukka-Pekka Kauppi, jukka-pekka.kauppi@tut.fi
% Tampere University of Technology

if corMat == 0
    if nargin ~= 5
        error('Number of input arguments must be 5!!')
        return
    end
    nrSubject = varargin{1};
end
if corMat == 1
    if nargin ~= 6
        error('Number of input arguments must be 6!!')
        return
    end
    winOn = varargin{1};
    if length(varargin) == 2
        timeInt = varargin{2};
    end
end
if corMat == 2
    if nargin ~= 4
        error('Number of input arguments must be 4!!')
        return
    end
end



Priv = Params.PrivateParams;
Pub = Params.PublicParams;

load([Pub.dataDestination 'memMaps'])
[pl,ms,en] = computer;

switch corMat
    case 0
        m = nrSubject;
        B = zeros([Priv.dataSize(nrSession,:)]);
        size(B)
        for xx = 1:Priv.dataSize(nrSession,1)
            if nrBand == 0
                if ~strcmp(en,Priv.computerInfo.endian)
                    A = swapbytes(memMaps.origMap.([Priv.prefixSession num2str(nrSession)]...
                        ).([Priv.prefixSubject num2str(m)]).Data(xx).tyz);
                else
                    A = memMaps.origMap.([Priv.prefixSession num2str(nrSession)]...
                        ).([Priv.prefixSubject num2str(m)]).Data(xx).tyz;
                end
            else
                if ~strcmp(en,Priv.computerInfo.endian)
                    A = swapbytes(memMaps.filtMap.([Priv.prefixSession num2str(nrSession)]...
                        ).([Priv.prefixSubjectFilt num2str(m)]).([Priv.prefixFreqBand ...
                        num2str(nrBand)]).Data(xx).tyz);
                else
                    A = memMaps.filtMap.([Priv.prefixSession num2str(nrSession)]...
                        ).([Priv.prefixSubjectFilt num2str(m)]).([Priv.prefixFreqBand ...
                        num2str(nrBand)]).Data(xx).tyz;
                end
            end
            %  disp(num2str(xx))
            B(xx,:,:,:) = shiftdim(A,4);
        end

    case 1
        if ~strcmp(en,Priv.computerInfo.endian)
            if winOn
                B = swapbytes(memMaps.cormatMap.win.([Priv.prefixFreqBand num2str(nrBand)]...
                    ).([Priv.prefixSession num2str(nrSession)]).cor.([Priv.prefixTimeVal num2str(timeInt)]).Data.xyzc);
            else
                B = swapbytes(memMaps.cormatMap.whole.([Priv.prefixFreqBand num2str(nrBand)]...
                    ).([Priv.prefixSession num2str(nrSession)]).cor.Data.xyzc);
            end
        else
            if winOn
                B = memMaps.cormatMap.win.([Priv.frefixFreqBand num2str(nrBand)]...
                    ).([Priv.prefixSession num2str(nrSession)]).cor.([Priv.prefixTimeVal num2str(timeInt)]).Data.xyzc;
            else
                B = memMaps.cormatMap.whole.([Priv.prefixFreqBand num2str(nrBand)]...
                    ).([Priv.prefixSession num2str(nrSession)]).cor.Data.xyzc;
            end
        end

    case 2

        B = zeros([Priv.dataSize(nrSession,:)]);
        size(B)
        for xx = 1:Priv.dataSize(nrSession,1)
            if ~strcmp(en,Priv.computerInfo.endian)
                A = swapbytes(memMaps.(Priv.phaseMapName).([...
                    Priv.prefixSession num2str(nrSession)]).([...
                    Priv.prefixFreqBand num2str(nrBand)]).Data(xx).tyz);
            else
                A = memMaps.(Priv.phaseMapName).([...
                    Priv.prefixSession num2str(nrSession)]).([...
                    Priv.prefixFreqBand num2str(nrBand)]).Data(xx).tyz;
            end
            B(xx,:,:,:) = shiftdim(A,4);
        end

    otherwise
        error('Invalid matrix type!')
        return
end

