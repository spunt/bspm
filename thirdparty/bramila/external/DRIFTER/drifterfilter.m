function filtered = drifterfilter( refdata)
%DRIFTERFILTER provides filtration of reference signal
%   refdata - D-dim array of structures this fields are requied
%              .name      = String with reference signal name
%              .dt        = sampling interval
%              .data      = 1xN vector
%              .downdt    = dt to downsample to
%              .freqlist  = List of possible frequencies for IMM
%                            in beats/min.
%              .downsampled = 1xM vector
%              .filter 0 or 1 for filtration only or amplitude normaisation
%	       .outfolder (optional for writing files in a separate folder)  

if ~isfield(refdata,'downdt'), refdata.downdt = refdata.dt; end
if ~isfield(refdata,'downsampled'), refdata.downsampled = refdata.data; end
if refdata.filter==1 || refdata.filter==0;
   %% FIR filtration
    % calculate frequencies to keep
    bpmLow=refdata.freqlist(1); % lower bpm
    bpmUp=refdata.freqlist(end);% upper bpm
    bpmHalf=0.5*(bpmUp-bpmLow);% middle frequency
    
    fs=1/refdata.downdt; %sampling frequency
    % values for bandpass filter parameter estim
    fLlow=(bpmLow-bpmHalf)/(60); %lower freq lowercut start
    fUup=(bpmUp+bpmHalf)/(60); %upper freq uppercut end
    fLow=(bpmLow)/(60); %lower freq lowercut end
    fUp=(bpmUp)/(60); %upper freq uppercut start
    a=[0,1,0]; % bandpass setting
    dev=[0.05 0.05 0.05]; % alowed error
    [n, Wn, beta, ftype]=kaiserord([fLlow, fLow, fUp, fUup], a, dev, fs); % fir paramter estimation
    b = fir1(n,Wn,ftype,kaiser(n+1,beta),'scale'); % filter
    refdata.filtered=filter(b,1,refdata.downsampled); % filtered data
    filtered=refdata.filtered;
    %% normalisation
    if refdata.filter==1
        %% bad region detect
        [localMax, Maxlocation]=findpeaks(refdata.filtered); % local maxima
        [localmin, minlocation]=findpeaks(-refdata.filtered); % local minima
        localmin=-localmin;
        
        meandat=mean(refdata.filtered); % mean of signal
        stddat=2*std(refdata.filtered); % 2*std of signal
        % selection of wrong/good data by signal mena and std
        goodMax=(localMax<meandat+stddat); %points that satisfy
        meanMax=mean(localMax(goodMax)); % mean of max
        goodmin=(localmin>meandat-stddat); %points that satisfy
        meanmin=mean(localmin(goodmin)); % mean of min
        
        % plot part 1
        h=figure('Visible', 'off');
        plot(refdata.filtered)
        hold on
        plot(Maxlocation(goodMax),localMax(goodMax),'g.',minlocation(goodmin), localmin(goodmin), 'y.')
        plot(Maxlocation(~goodMax),localMax(~goodMax),'r.',minlocation(~goodmin), localmin(~goodmin), 'r.')
        
        
        %% padding wrong parts with artificial signal
        % artificial parts are at the same frequency as original signal only
        % amplitude normalised
        localMax(~goodMax)=meanMax;
        localmin(~goodmin)=meanmin;
        refdata.filtered(Maxlocation)=localMax; % set wrong peaks to local Max mean
        refdata.filtered(minlocation)=localmin; % set wrong peaks to local min mean
        % mask for interpoation
        refdata.interp=zeros(size(refdata.filtered)); % not a peak = 0
        location=sort([Maxlocation, minlocation]); % peak positions
        refdata.interp(location)=1; % mask - peaks = 1
        refdata.interp(Maxlocation(~goodMax))=-1; % wrong peaks = -1
        refdata.interp(minlocation(~goodmin))=-1; % wrong peaks = -1
        
        %  interpolation
        for j=2:numel(location)-1; %go trough all peaks
            if refdata.interp(location(j))==-1 % if bad interpolate between start, mid,  end
                Start=location(j-1);
                Mid= location(j);
                End=location(j+1);
                refdata.filtered(Start:End)=interp1([Start, Mid, End], ...
                    [refdata.filtered(Start),refdata.filtered(Mid),...
                    refdata.filtered(End)],Start:End);
            end
        end
        %% final filtration
        refdata.filtered=filter(b,1,refdata.filtered); % filtered data
        
        % plot part 2
        plot(refdata.filtered,'r')
        plot([1,Maxlocation(end)],[meanMax, meanMax],'g',[1,Maxlocation(end)],[meanmin, meanmin],'y')
	
	% new code to store plots in different folder
	[temp_path,temp_name,temp_ext]=fileparts(refdata.name);

	if(isfield(refdata,'outfolder'))
		outpath=refdata.outfolder;
	else
		outpath=temp_path;
	end
	% check if we can write in that folder
	[temp_status, temp_permissions]=fileattrib(outpath);
	if(temp_permissions.UserWrite~=1)
		error(['DRIFTER cannot write in folder: ' outpath]);
	else
		disp(['Storing DRIFTER plots in ' outpath])
	end

        %saveas(h,[refdata.name,'_filter'],'jpg')
        saveas(h,fullfile(outpath,[temp_name  '_filter']),'jpg')
        filtered=refdata.filtered;
    end
end
end



