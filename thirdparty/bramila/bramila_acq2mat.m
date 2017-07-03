function filename=bramila_acq2mat(biopacfile)
% acq2mat - reads acq file from biopac and store it in drifter-friendly
%       format data are saved as filename_partN.mat into same folder as
%       filename.acq
%
% Input:
%       biopacfile - structure
%       biopacfile.name - path to biopac data folder with 'filename.acq' files
%       channel inputs are 0 or 1:
%       biopacfile.CH1 - logical[1,1] usually breath data from belt
%       biopacfile.CH2 - logical[1,1] usually ECG
%       biopacfile.CH3 - logical[1,1] usually GSR
%       biopacfile.CH4 - logical[1,1] usually heart rate from pulse oxymeter_left
%       biopacfile.CH5 - logical[1,1] usually heart rate from pulse oxymeter_right
%       biopacfile.CH35 - logical[1,1] usually MRI scan off/on information
%       biopacfile.dt - double[1,1] sampling interval in seconds/sample
%       biopacfile.dtMRI - double[1,1] sampling niterval of fMRI in seconds/sample
%       biopacfile.controlplot - logical[1,1] 0/1 plots all data in control
%       plot
%       biopacfile.filter - logical[1,1] 0/1 filters biopac data with FIR
%       bandpass filter
% Output:
%       filename - cell of strings [1,N] filenames of .mat files containig
%       data which coresponds to epi files, N - number of parts, same as number
%       of epi files (eg. stimuli settings)
% Done by Lenka: lenka.vondrackova@aalto.fi
channels=zeros(1,6);
channels(1)=biopacfile.CH1;
channels(4)=biopacfile.CH4;
channels(5)=biopacfile.CH5;
channels(6)=biopacfile.CH35;
channels=logical(channels);

%% load acq file
acq = load_acq(biopacfile.name);

%% channel setting

if acq.hdr.graph.num_channels~=6
    dat=zeros(size(acq.data,1),6);
    headerx={'N/A','N/A','N/A','N/A','N/A','N/A'};
    for i=1:acq.hdr.graph.num_channels;
        ind=acq.hdr.per_chan_data(i).num;
        if ind==35
            ind=6;
        elseif ind>5;
            continue
        end
        dat(:, ind)=acq.data(:,i);
        headerx{ind}=acq.hdr.per_chan_data(i).comment_text;
    end
    acq.data=dat;
    
else
    for i=1:acq.hdr.graph.num_channels
        headerx{i}=acq.hdr.per_chan_data(i).comment_text
    end
end

%visualize
if biopacfile.controlplot==1
    disp('plotting control plot...')
    h=figure('Visible','Off');
    plot(acq.data)
    legend(headerx)
    title(biopacfile.name)
    saveas(h,[biopacfile.name(1:end-4) '_controlplot.png'],'png');
end

%% reduce data to be used with drifter
data=acq.data(:,channels);
header=headerx(channels);

clear headerx
clear acq.data
% % Commented it out, because on cluster there is no plotting
% figure
% plot(data)
% legend(header)
% title(biopacfile.name)


if biopacfile.CH35==0;
    for i=1:sum(channels)
        refdata{i}.data=data(:,i);
        refdata{i}.dt=biopacfile.dt;
        refdata{i}.header=header(i);
    end
    warning(['CH35 was set to zero, cannot split measurements!!! Data from other' ...
        ' chanels are stored in: ' biopacfile.name(1:end-4) '_refdata.mat'])
    filename=[biopacfile.name(1:end-4) '_refdata.mat'];
    save(filename, 'data')
    
else
    %% run only if MRI scanner active == CH35 active
    disp('counting measurements...')
    %extract number of measurements, starting time and sapmlig freq of MRI
    scanneractive=find(data(:,3)~=0);
    if(isempty(scanneractive))
        disp('It looks like we do not have markers for the measurements, we cannot process these data')
        filename='';
        return
    end
    times=diff(scanneractive(diff(scanneractive)~=1));
    mn=mean(times);
    dtMRI=round(mean(times(times<=mn)));
    startOfMeasurement=scanneractive([1; find(diff(scanneractive)>mn)+1]);
    dt=acq.hdr.graph.sample_time/1000;
    if dt~=biopacfile.dt;
        warning(['diference between calculated sampling time(%i)'...
            ,' and set sampling time(%i)'],dt,biopacfile.dt)        
    else
        dt=biopacfile.dt;
    end
    fprintf(['Period of fMRI scan is : %f samples, \n' ...
        'it is %f seconds\n'],dtMRI, dtMRI*dt)
    numOfMeasurements=size(startOfMeasurement,1);
    
    %% parse measurements
    fprintf('parsing into %d measurements...\n', numOfMeasurements)
    volume=zeros(1,numOfMeasurements);
    
    if numOfMeasurements==1
        volume=sum(diff(data(startOfMeasurement:end,3))>1);
        parsed=data(startOfMeasurement:startOfMeasurement+(volume*dtMRI),:);
        refdata{1}.data=parsed(:,1);
        refdata{1}.dt=dt;
        refdata{2}.data=parsed(:,2);
        refdata{2}.dt=dt;
        filename{1}=[biopacfile.name(1:end-4) '.mat'];
        save(filename{1},'refdata');
        disp(['saved as: ' filename{1}])
    else
        for i=1:numOfMeasurements
            if i~=numOfMeasurements
                volume(i)=sum(diff(data(startOfMeasurement(i):startOfMeasurement(i+1),3))>1);
            else
                volume(i)=sum(diff(data(startOfMeasurement(i)-1:end,3))>1);
            end
            % if vol==vol....
            parsed=data(startOfMeasurement(i)-1:startOfMeasurement(i)+(volume(i)*dtMRI),:);
            refdata{1}.data=parsed(:,1);
            refdata{1}.dt=dt;
            refdata{2}.data=parsed(:,2);
            refdata{2}.dt=dt;
            filename{i}=[biopacfile.name(1:end-4) '_part_' num2str(i) '.mat'];
            save(filename{i}, 'refdata');
            disp(['saved as: ' filename{i}])
        end
    end
end

%% lincense for load_acq
% Copyright (c) 2009, Jimmy Shen
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
