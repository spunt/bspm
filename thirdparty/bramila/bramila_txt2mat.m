function filename=txt2mat(biopacfile)
% txt2mat - reads txt file from biopac and store it in drifter-friendly 
%       format data are saved as filename_partN.mat into same folder as
%       filename.txt
%  
% Input:
%       biopacfile - structure
%       biopacfile.name - path and filename of biopac data eg.'measurement_1-4.txt'
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

%% settings of channels
% channel 2 and 3 can be switched on, so txt2 mat can work also as parser
% current script doesnt support this settings, script can also work more
% with header file to load number of chanels and then work with data trough
% channel names 'PPG', 'Respiratory' etc. but now the settings is rigid and
% for processing are used numbers of columns [1, 4 or 5, 6];

channels=zeros(1,6);
channels(1)=biopacfile.CH1;
channels(4)=biopacfile.CH4;
channels(5)=biopacfile.CH5;
channels(6)=biopacfile.CH35;
channels=logical(channels);

%open file and read channels
disp('reading from file ...')
FID=fopen(biopacfile.name);
sf=textscan(FID,'%f%s',1, 'headerlines',1);
xheader=textscan(FID,'%s',6, 'headerlines',14);
N=textscan(FID,'%f %f %f %f %f %f',1);
x=textscan(FID,'%f %f %f %f %f %f');
fclose(FID);
data=cell2mat(x);
header=cellstr(xheader{1});
if biopacfile.controlplot==1
    disp('plotting control plot...')
    h=figure
    plot(data)
    legend(header)
    savefig(h,[biopacfile.filename(1:end-4) '_controlplot']);
end
header=header(channels);
data=data(1:min(cell2mat(N))-1,channels);

%visualize

figure
plot(data) 
legend(header)
title(biopacfile.name(1:end-4))


%run only if MRI scanner active CH35 active
if biopacfile.CH35==0;
    for i=1:sum(channels)
        refdata{i}.data=data(:,i);
        refdata{i}.dt=biopacfile.dt;
        refdata{i}.header=header(i);
    end
    warning(['CH35 was set to zero, cannot split measurements!!! Data from other' ...
        ' chanels are stored in: ' biopacfile.name(1:end-4) '_refdata.mat'])
    filename=[biopacfile.name(1:end-4) '_refdata.mat'];
    save(filename, 'refdata')
    
else
    disp('counting measurements...')
    %extract number of measurements, starting time and sapmlig freq of MRI
    scanneractive=find(data(:,3)~=0);
    times=diff(scanneractive(diff(scanneractive)~=1));
    mn=mean(times);
    dtMRI=round(mean(times(times<mn)));
    startOfMeasurement=scanneractive([1; find(diff(scanneractive)>mn)+1]);
    if strcmp(sf{2},cellstr('msec/sample'))
        dt=sf{1}/1000;
        if dt~=biopacfile.dt;
            warning(['diference between calculated sampling time(%i)'...
                ,' and set sampling time(%i)'],dt,biopacfile.dt)
        end
    else
        dt=biopacfile.dt;
    end
    fprintf(['Period of fMRI scan is : %f samples, \n' ...
        'it is %f seconds\n'],dtMRI, dtMRI*dt)
    numOfMeasurements=size(startOfMeasurement,1);
    
    %parse measurements
    fprintf('parsing into %d measurements...\n', numOfMeasurements)
    volume=zeros(1,numOfMeasurements);
    
    if numOfMeasurements==1
        volume=sum(diff(data(startOfMeasurement:end,3))>1);
        parsed=data(startOfMeasurement:startOfMeasurement+(volume*dtMRI),:);
        refdata{1}.data=parsed(:,1);
        refdata{1}.dt=dt;
        refdata{2}.data=parsed(:,2);
        refdata{2}.dt=dt;
        bramila_biopacFilter(refdata);
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
            bramila_biopacFilter(refdata);
            fn=sprintf('_refdata_part_%d.mat',i);
            filename{i}=[biopacfile.name(1:end-4) fn];
            save(filename{i},'refdata');
            disp(['saved as: ' filename(i)])
        end
    end
end

end
