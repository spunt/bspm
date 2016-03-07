function SPM8_FirstLevel(fn,pname,ord)
%%% Look at the embedded DefaultParams function (at the bottom of this 
%%% script) for paramter and sequence specification.  The DefaultParams 
%%% function will be called when there is no Params.m file in the working 
%%% directory.  This script was designed to be used with SPM8.
%%%
%%% All configuration options should be done in the Params m-file.
%%%
%%% FN: is a Session Name / Folder Name where the data for a particular
%%% Session can be found.
%%%
%%% PNAME: is the name of the Params file that you want to use.  If this is
%%% not specified it will default to Params.m
%%%
%%% ORD:  This is an optional input that you can use to respcify the order
%%% of processing without changing the order in the Params.m file. This
%%% will specify the order of execution of P.PO (see the Parameter m-file).
%%%
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%% Copyright (C) 2014,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.
P = [];
if nargin == 0 || isempty(fn);
    fn = setdir;
end

if nargin<2
    pname = 'Params';
    if exist(pname)
        eval(['P = ' pname ';']);
        pname = 'Params';
    else
        disp('Warning! No Params.m file was found.  A new Params.m file has been written in the current directory.  Please check the settings in this file, and then rerun the script.');
        DefaultParams;
        edit Params.m
        return
    end
else
    disp(['P = ' pname ';'])
    eval(['P = ' pname ';']);
end

if ~exist([P.root filesep fn filesep P.DestFold]);
    mkdir([P.root filesep fn filesep P.DestFold]);
end
cd([P.root filesep fn filesep P.DestFold]);

if exist([pname '.mat'])==0;
    save([pname '.mat'],'P');
end

if nargin<3
   ord = 1:length(P.PO); 
end

tmp = P.root;
cd(P.root);


for qq = 1:length(ord)
    eval(P.PO{ord(qq)});
end

try
    P1 = load([P.root filesep fn filesep P.DestFold filesep pname '.mat']);
    P1.P.UserTime = UserTime;
    P1.P.root = P.root;
    P1.P.DestFold = P.DestFold;
    P1.P.TR = P.TR;
    P1.P.List = P.List;
    P1.P.BadRuns = P.BadRuns;
    P1.P.Runs = P.Runs;
    P1.P.nVols = P.nVols;
    P = P1.P;
catch
    cd(P.root);
    disp('Something Went Wrong With This Analysis!');
    fprintf('\n\n\n\n');
    return
end
 
save([P.root filesep fn filesep P.DestFold filesep pname '.mat'],'P');
    
    function Bad_Vols(fn)
        fprintf('\n%s\n', 'Examing Image Files and Creating Bad Volume Regressors:');
        eval(['P = ' pname ';']);
        if isempty(P.nVols)
            ch = ReadInFile([P.root filesep fn filesep P.DestFold '/' P.List],'\t');
            tmp = dir_wfp([P.root filesep fn filesep P.bv.SourceFold filesep P.bv.SourcePrefix '*' P.bv.SourceSuffix]);
            for ii = 1:length(ch); if strcmpi(ch{ii},'NA'); P.nVols(ii) = NaN;  else; P.nVols(ii) = numel(spm_vol(tmp{contains(ch{ii},tmp)})); end; end;
        end
        %if ~isempty(contains('aschultz',{UserTime})); keyboard; end
        
        warning off
        delete([P.root filesep fn filesep P.DestFold filesep 'bad_runs.txt']);
        delete([P.root filesep fn filesep P.DestFold filesep 'NoGo*.mat']);
        warning on
        
        nnn = [P.root filesep fn filesep P.DestFold filesep P.bv.LogFileName];
        S = [];
        WriteDataToText(S, nnn, 'w', '\t');
        
        if P.bv.do == 1;
            P1 = load([P.root filesep fn filesep P.DestFold filesep pname '.mat']);
            cd([P.root filesep fn filesep P.DestFold]);
            
            %%% Get List of Run Identifiers and create the RunIndex;
            list = ReadInFile([P.root filesep fn filesep P.DestFold filesep P.List],'\t');
            list = list(:,1);
            origList = list;
            which = contains('NA',list);
            if exist('RunInfo.csv'); 
                dd = ReadInFile('RunInfo.csv',',',1);
                missDat = unique([dd{2:end,1}]);
                br = setdiff(1:length(list),missDat);
                
                for ii = 1:length(br)
                    disp(['Junking Run ' num2str(br(ii)) ' due to missing task data.']);
                    S =[];
                    save(['NoGo_MissingTaskData' sprintf('%0.2d',br(ii)) '.mat'], 'S');
                    S{1} = ['Junking Runs ' num2str(br(ii)) ' due to missing task data.'];
                    WriteDataToText(S, nnn, 'a', '\t');
                end
                
                dlmwrite([P.root filesep fn filesep P.DestFold filesep P.BadRuns],br,'delimiter','\t');
                br = [];
            end
            
            try
                br = load([P.root filesep fn filesep P.DestFold filesep P.BadRuns]);
            catch
                br = [];
            end
            
            if isfield(P.cm,'IgnoreBadRuns')
                if P.cm.IgnoreBadRuns == 1
                    br = [];
                end
            end
            
            which = unique([br which]);
            ind = setdiff(1:length(list),which);
            
            %%% Experimental: added on 08-24ish-2011
            %ind = intersect(ind,P.Runs);
            
            list = list(ind);
            RunIndex = ind;
            
            %%% Collect Global SNR Values, for writing to analysis log and
            %%% screening for bad runs.
            if P.bv.CheckGlobalSNR
                %%% Get list of file names that have the SNR data
                [ll l2] = dir_wfp([P.root filesep fn filesep P.bv.SNR_SourceFold filesep P.bv.SNR_prefix '*' P.bv.SNR_suffix]);
                ind = [];

%                 try
                for ii = 1:length(list); ind(ii) = contains(['.*' list{ii} '.*'],l2); end
%                 catch
%                    keyboard;
%                 end
                ll = ll(ind); l2 = l2(ind);
                %%% Collect the SNR infro and write it to the Analysis Log.
                
                snr = [];
                for ii = 1:length(ll);
                    snr = [snr; dlmread(ll{ii},'\t',1,1)];
                end
                disp(snr)
                S = [];
                S(1).a1 = 'From Run:';
                S(1).a2 = 'Mean:';
                S(1).a3 = 'SD:';
                S(1).a4 = 'SNR:';
                for ii = 1:size(snr,1);
                    S(ii+1).a1 = l2{ii};
                    S(ii+1).a2 = snr(ii,1);
                    S(ii+1).a3 = snr(ii,2);
                    S(ii+1).a4 = snr(ii,3);
                end
                %keyboard
                WriteDataToText(S, nnn, 'a', '\t');
                %%% Screen SNR values for abnormally low values.
                %keyboard;
%                 try
                ind = find(snr(:,3)<P.bv.SNRthresh | isnan(snr(:,3)));
%                 catch
%                     keyboard;
%                 end
                if numel(ind)>0;
                    disp(['Junking Runs ' num2str(RunIndex(ind)) ' (' l2{ind} ')' ' due to poor SNR values.']);
                    S =[];
                    save(['NoGo_BadSNR' sprintf('%0.2d',RunIndex(ind)) '.mat'], 'S');
                    S.a1 = ['Junking Runs ' num2str(ind') ' (' l2{ind} ')' ' due to poor SNR values.'];
                    WriteDataToText(S, nnn, 'a', '\t');
                    try
                        tmp = load([P.root filesep fn filesep P.DestFold filesep P.BadRuns]);
                    catch
                        tmp = [];
                    end
                    
                    tmp = sort(unique([tmp(:); ind(:)]))';
                    dlmwrite([P.root filesep fn filesep P.DestFold filesep P.BadRuns],tmp,'delimiter','\t');
                end

                S = [];
                S.a1 = '  ';
                WriteDataToText(S, nnn, 'a', '\t');
            end
            
            %try br = load([P.root filesep fn filesep P.DestFold filesep P.BadRuns]); catch; br = []; end %#ok<CTCH>
            %RunIndex = setdiff(RunIndex,br);
            %list = origList(setdiff(RunIndex,br));
            
            %%% Get the List of Input files for the 2nd level analysis as
            %%% well as the realignmnet parameters for each run.
            
            %%% Make sure that everything matches up.
            [nn nn_1] = dir_wfp([P.root filesep fn filesep P.bv.SourceFold filesep P.bv.SourcePrefix '*.nii']);
            %keyboard
            for ii = 1:length(list); i1(ii) = contains(['.*' list{ii} '.*'],nn_1); end
            [nn2 nn2_1] = dir_wfp([P.root filesep fn filesep P.bv.RealignDir filesep P.bv.RealignPrefix '*.txt']);

            for ii = 1:length(list); i2(ii) = contains(['.*' list{ii} '.*'],nn2_1); end
            %%% Look for bad volumes using both the global signal and the
            %%% movement parameters
            count = 0;
            All = [];
            All2 = [];
            for kk = 1:length(list);
                flag = 1;
                count = count+1;
                x = i1(kk);
                x2 = i2(kk);
                %keyboard;
                disp(['Examining Run ' num2str(RunIndex(kk)) ': ' nn_1{x}]);
                
                MM = openIMG(nn{x});
                GM = [];
                for jj = 1:size(MM,4);
                    M = MM(:,:,:,jj);
                    indd = find(M>((mean(M(1:end)))/8));
                    GM(jj) = mean(M(indd));
                end
                GM1 = diff(GM);
                ind1 = find(abs(GM1)>(std(GM1)*P.bv.GlobalSigThresh))+1;
                
                dat = load(nn2{x2});
                
                if isfield(P.bv, 'MeanMovementThresh')
                    mv = mean(sqrt(sum(diff(dat(:,1:3)).^2,2)));
                    WriteDataToText({['Mean Movement = ' sprintf('%0.4f',mv)]}, nnn, 'a', '\t');
                    disp(['Mean Movement for Run ' num2str(RunIndex(kk)) ' = ' sprintf('%0.4f',mv)]);
                    if mv > P.bv.MeanMovementThresh
                        disp(['Junking Run ' num2str(RunIndex(kk)) ' due to excessive mean movement.']);
                        S = {''; ['Junking Run ' num2str(RunIndex(kk)) ' due to excessive mean movement. MV = ' sprintf('%0.3f',mv)]};
                        WriteDataToText(S, nnn, 'a', '\t');
                        save(['NoGo_ExcessiveMeanMovement' sprintf('%0.2d',RunIndex(kk)) '.mat'], 'S');
                        try tmp = load([P.root filesep fn filesep P.DestFold filesep P.BadRuns]); catch; tmp = []; end
                        if size(tmp,1)>size(tmp,2); tmp = tmp'; end
                        tmp = sort(unique([tmp RunIndex(kk)]));
                        dlmwrite([P.root filesep fn filesep P.DestFold filesep P.BadRuns],tmp,'delimiter','\t');
                        WriteDataToText({' '}, nnn, 'a', '\t');
                        flag = 0;
                    end
                    WriteDataToText({' '}, nnn, 'a', '\t');
                end
                
                tm = DistMat(dat(:,1:3),0);
                tm = max(tm(:));
                tr = DistMat(dat(:,4:6),0);
                tr = max(tr(:));

                S = {''; ['Total Movement = ' num2str(round(1000*(tm))/1000) ' milimeters']; ['Total Rotation = ' num2str(round(1000*(tr))/1000) ' degrees']};
                if tm > P.bv.TotalMovementThresh
                    disp(['Junk run # ' num2str(RunIndex(kk)) ' due to too much overall movement.']);
                    S{end+1} = ['Junk run # ' num2str(RunIndex(kk)) ' due to too much overall movement.'];
                    try tmp = load([P.root filesep fn filesep P.DestFold filesep P.BadRuns]); catch; tmp = []; end
                    if size(tmp,1)>size(tmp,2); tmp = tmp'; end
                    tmp = sort(unique([tmp RunIndex(kk)]));
                    dlmwrite([P.root filesep fn filesep P.DestFold filesep P.BadRuns],tmp,'delimiter','\t');
                    S{end+1} = '';
                    disp(S{2}); disp(S{3});
                    WriteDataToText(S, nnn, 'a', '\t');
                    flag = 0;
                end
                if tr > P.bv.TotalRotationThresh
                    disp(['Junk run # ' num2str(RunIndex(kk)) ' due to too much overall rotation.']);
                    S{end+1} = ['Junk run # ' num2str(RunIndex(kk)) ' due to too much overall rotation.'];
                    try tmp = load([P.root filesep fn filesep P.DestFold filesep P.BadRuns]); catch; tmp = []; end
                    if size(tmp,1)>size(tmp,2); tmp = tmp'; end
                    tmp = sort(unique([tmp RunIndex(kk)]));
                    dlmwrite([P.root filesep fn filesep P.DestFold filesep P.BadRuns],tmp,'delimiter','\t');
                    S{end+1} = '';
                    disp(S{2}); disp(S{3});
                    WriteDataToText(S, nnn, 'a', '\t');
                    flag = 0;
                end
                S{end+1} = '';
                disp(S{2}); disp(S{3});
                WriteDataToText(S, nnn, 'a', '\t');
                
                tmp = diff(dat).^2;
                pos = sqrt(sum(tmp(:,1:3),2));
                ori = sqrt(sum(tmp(:,4:6),2)) * (360/(2*pi));
                %if kk ==6; keyboard; end;
                ind2 = find(pos>P.bv.MoveThresh)+1;
                ind3 = find(ori>P.bv.RotThresh)+1;

                ind = unique([ind1'; ind2; ind3]);
                
                if numel(ind)>P.bv.BadVolThresh
                    disp(['Junking Runs ' num2str(RunIndex(kk)) '(' nn_1{x} ')' ' due to too many bad volumes.']);
                    try
                        tmp = load([P.root filesep fn filesep P.DestFold filesep P.BadRuns]);
                    catch
                        tmp = [];
                    end
                    
                    tmp = sort(unique([tmp RunIndex(kk)]));
                    dlmwrite([P.root filesep fn filesep P.DestFold filesep P.BadRuns],tmp,'delimiter','\t');
                    S = [];
                    save(['NoGo_TooManyBadVolumes' sprintf('%0.2d',RunIndex(kk)) '.mat'], 'S');
                    S.a1 = ['Junking Runs ' num2str(RunIndex(kk)) '(' nn_1{x} ')' ' due to too many bad volumes.'];
                    WriteDataToText(S, nnn, 'a', '\t');
                    flag = 0;
                end
                
                tim2 = (ind)+sum(P.nVols(1:RunIndex(count)-1));
                %All = [All; tim2];
                %tim3 = (ind)+sum(P.nVols(RunIndex(1:count-1)));
                %All2 = [All2; tim3];
                S = [];
                S.a1 = ['Negating ' num2str(length(tim2)) ' volumes from run ' num2str(RunIndex(kk)) '; With local indices of: '  num2str(ind') ' and global indices of: ' num2str(tim2')];
                WriteDataToText(S, nnn, 'a', '\t');
                
                %if flag
                R{RunIndex(kk)} =  zeros(P.nVols(RunIndex(kk)),numel(ind));
                for qq = 1:length(ind);
                    R{RunIndex(kk)}(ind(qq),qq) = 1;
                end
                %end
            end
            S = [];
            S.a1 = '  ';
            WriteDataToText(S, nnn, 'a', '\t');
            type(nnn);
            
            %keyboard; 
            save([P.root filesep fn filesep P.DestFold filesep P.bv.BadVolRegsName], 'R')                        

            %%% Finish up.
            P1.P.bv = P.bv;
            P1.P.bv.UserTime = UserTime;
            P = P1.P;
            save([P.root filesep fn filesep P.DestFold filesep pname '.mat'], 'P');
            
            cd(P.root);
        else
            cd(P.root);
            fprintf('\n%s\n', 'Skipping run/volume assessment.');
        end
    end

    function Create_Model(fn)
        fprintf('\n%s\n', 'Create and Estimate the First Level Model:');
        eval(['P = ' pname ';']);
        
        %if ~isempty(contains('aschultz',{UserTime})); keyboard; end
        
        if isempty(P.nVols)
            ch = ReadInFile([P.root filesep fn filesep P.DestFold '/' P.List],'\t');
            tmp = dir_wfp([P.root filesep fn filesep P.bv.SourceFold filesep P.bv.SourcePrefix '*' P.bv.SourceSuffix]);
            for ii = 1:length(ch); if strcmpi(ch{ii},'NA'); P.nVols(ii) = NaN;  else; P.nVols(ii) = numel(spm_vol(tmp{contains(ch{ii},tmp)})); end; end;
        end
        
        if isempty(P.Runs)
            P.Runs = 1:numel(P.nVols);
        end
        
        if P.cm.do == 1;
            P1 = load([P.root filesep fn filesep P.DestFold filesep pname '.mat']);

            cd([P.root filesep fn filesep P.DestFold]);

            list = ReadInFile([P.root filesep fn filesep P.DestFold filesep P.List],'\t');
            list = list(:,1);

            which = contains('NA',list);
            try
                br = load([P.root filesep fn filesep P.DestFold filesep P.BadRuns]);
            catch
                br = [];
            end
            
            if isfield(P.cm,'IgnoreBadRuns')
                if P.cm.IgnoreBadRuns == 1
                    br = [];
                end
            end
            
            which = unique([which br]);
            RunIndex = setdiff(P.Runs,which);
            
            if isfield(P.cm, 'ngr')
                if numel(RunIndex)<P.cm.ngr;
                    warning('Not enough good runs to analyze data!');
                    return
                end
            else
                if numel(RunIndex)<3
                    warning('Not enough good runs to analyze data!');
                    %return;
                    %cd(P.root);
                end
            end
            
            [ll l2] = dir_wfp([P.root filesep fn filesep P.cm.SourceFold filesep P.cm.SourcePrefix '*' P.cm.SourceSuffix]);
            
            ind = [];
            for ii = RunIndex;
                ind(ii) = contains(['.*' list{ii} '.*'],l2);
            end
            ll = ll(ind(ind~=0));
            P.cm.SourceList = l2(ind(ind~=0));
            
            SPM.xY.RT=P.TR;
            SPM.xY.P = char(ll);
            
            %%% Put together SPM.xBF
            spm('defaults','FMRI')
            % The following lines have the side effect of modifying the global
            % defaults variable. This is necessary to pass job.timing.fmri_t to
            % spm_hrf.m. The original values are saved here and restored at the end
            % of this function, after the design has been specified. The original
            % values may not be restored if this function crashes.
            global defaults
            olddefs.stats.fmri.fmri_t=spm_get_defaults('stats.fmri.fmri_t');
            olddefs.stats.fmri.fmri_t0=spm_get_defaults('stats.fmri.fmri_t0');
            defaults.stats.fmri.t =  P.cm.MicrotimeRes;
            defaults.stats.fmri.t0 = P.cm.MicrotimeOnset;
            SPM.xBF.UNITS =          P.cm.Units;
            SPM.xBF.dt =             P.TR/P.cm.MicrotimeRes;
            SPM.xBF.T     =          P.cm.MicrotimeRes;
            SPM.xBF.T0    =          P.cm.MicrotimeOnset;
            
            if strcmpi('hrf',P.cm.basis_func)
                SPM.xBF.name = 'hrf';
            elseif strcmpi('hrf (with time derivative)',P.cm.basis_func)
                SPM.xBF.name = 'hrf (with time derivative)';
            elseif strcmpi('fir',P.cm.basis_func)
                SPM.xBF.name = 'Finite Impulse Response';
                SPM.xBF.length  = P.cm.length;
                SPM.xBF.order   = P.cm.order;
                %keyboard
            end
            
            
            
            SPM.xBF = spm_get_bf(SPM.xBF);
            
            SPM.xBF.Volterra = P.cm.Volterra;
            
            %keyboard;
            SPM.nscan = sum(P.nVols(RunIndex));

            SPM.xGX.iGXcalc = P.cm.iGXcalc;
            SPM.xGX.sGXcalc = P.cm.sGXcalc;
            SPM.xGX.sGMsca =  P.cm.sGMsca;

            SPM.xVi.form = P.cm.AR;
            
            for ii = 1:length(RunIndex);
                if numel(P.cm.HP_filt)==1
                    SPM.xX.K(ii).HParam = P.cm.HP_filt;
                else
                    SPM.xX.K(ii).HParam = P.cm.HP_filt(RunIndex(ii));
                end
            end
            %%% Will need something in here to choose between hard coded
            %%% onset and duration information and a config file.

            if strcmpi(P.cm.OnsetsFile(end-3:end),'.csv')    
                dat = ReadInFile([P.root filesep fn filesep P.DestFold filesep P.cm.OnsetsFile],',',1);
            else
                dat = ReadInFile([P.root filesep fn filesep P.DestFold filesep P.cm.OnsetsFile],'\t',1);
            end
            %keyboard
            t1 = [dat{2:end,1}];  %%% First Column is the Run#
            Cond = dat(2:end,2);  %%% Second Column is Condition Name
            t2 = [dat{2:end,3}];  %%% Third Column is Onset.
            t3 = [dat{2:end,4}];  %%% Fourth Column is Duration.
            
            c = 0;
            for ii = RunIndex;
                c = c+1;
                SPM.nscan(c) = P.nVols(ii);
                i1 = find(t1 == ii);
                %try
                uni = unique(Cond,'stable');
                %catch
                    %keyboard
                %end
                
                if strcmpi(P.cm.Units ,'scans')
                    mx = P.nVols(ii);
                else
                    mx = P.nVols(ii)*P.TR;
                end
                
                cc = 0;
                for jj = 1:length(uni);
                    i2 = contains(uni{jj}, Cond(i1));
                    if isempty(i2);
                        continue;
                    end
                    %if ~isempty(contains('aschultz',{UserTime})); keyboard; end;
                    if all(t2(i1(i2))<mx)
                        cc = cc+1;
                        SPM.Sess(c).U(cc).name{1} = uni{jj};
                        SPM.Sess(c).U(cc).ons =     t2(i1(i2));
                        SPM.Sess(c).U(cc).dur =     t3(i1(i2));
                        SPM.Sess(c).U(cc).P.name =  P.cm.TempMod;
                        SPM.Sess(c).U(cc).P.h = 0;      %%?
                        SPM.Sess(c).C.C = [];           %%?
                        SPM.Sess(c).C.name = cell(0);   %%?
                    else
                        %keyboard;
                        error('Onsets go beyond the end of the run');
                    end
                end
            end
            
            c = 0;
            if P.cm.addBadVolRegs == 1;
                if exist([P.root filesep fn filesep P.DestFold filesep P.cm.BadVolRegs])>0
                    X = load([P.root filesep fn filesep P.DestFold filesep P.cm.BadVolRegs]);
                    %keyboard;
                    for ii = RunIndex
                        c = c+1;
                        tmp = X.R{ii};
                        for jj = 1:size(tmp,2);
                            SPM.Sess(c).C.name{end+1} = ['bv_' num2str(jj)];
                        end
                        SPM.Sess(c).C.C = tmp;
                    end
                end
            end
            if P.cm.addMotionRegressors == 1;
                c = 0;
                for ii = RunIndex
                    c = c+1;
                    [l1 l2] = dir_wfp([P.root filesep fn filesep P.cm.MotionFold filesep P.cm.MotionPrefix '*.txt']);
                    mot = load(l1{contains(['.*' list{ii} '.*'],l2)});
                    for kk = 1:size(mot,2);
                        SPM.Sess(c).C.name{end+1} = ['mot_' num2str(kk)];
                    end
                    SPM.Sess(c).C.C = [SPM.Sess(c).C.C mot];
                end
            end

            mkdir([P.root filesep fn filesep  P.cm.DestFold]);
            cd([P.root filesep fn filesep  P.cm.DestFold]);
            
                       
            delete('beta*');
            delete('mask*');
            delete('ResMS*');
            delete('RPV*');
            delete('SPM.mat');
            delete('ess*'); 
            delete('con*');
            delete('spmF*');
            delete('spmT*');

            SPM = spm_fmri_spm_ui(SPM); %% Will cause compiled script to fail.
            
            if P.cm.ScreenTaskCorrMot
                S = [];
                ccc = 0;
                mot = [];
                for ii = RunIndex;
                    ccc = ccc+1;
                    [l1 l2] = dir_wfp([P.root filesep fn filesep P.cm.MotionFold filesep P.cm.MotionPrefix '*.txt']);
                    mot = load(l1{contains(['.*' list{ii} '.*'],l2)});
                    
                    tmp = unique(Cond);
                    xx = [];
                    c = 0;
                    for jj = 1:length(tmp);
                        th1 = contains(['.*' num2str(ccc) '.*' tmp{jj} '.*bf\(1\).*'],SPM.xX.name);
                        if isempty(th1);
                            continue;
                        end
                        rows =  1+sum(P.nVols(1:ccc-1)):sum(P.nVols(1:ccc));
                        c = c+1;
                        xx(:,c) = SPM.xX.X(rows,th1);
                    end
                    
                    sig = mean(xx,2);
                    ind = find(sig~=0);
                    coLin{ii} = corr(sig(ind),mot(ind,:)).^2;
                end
                
                %%%  Write This data to the Analysis Log
                S = [];
                ccc = 1;
                S(ccc).Label = [' '];
                S(ccc).Corr1 = 'X';
                S(ccc).Corr2 = 'Y';
                S(ccc).Corr3 = 'Z';
                S(ccc).Corr4 = 'Yaw';
                S(ccc).Corr5 = 'Pitch';
                S(ccc).Corr6 = 'Roll';
                for ii = RunIndex;
                    ccc = ccc+1;
                    S(ccc).Label = ['Run #' num2str(ii) ' (' P.cm.SourceList{ccc-1} ')'];
                    S(ccc).Corr1 = coLin{ii}(1);
                    S(ccc).Corr2 = coLin{ii}(2);
                    S(ccc).Corr3 = coLin{ii}(3);
                    S(ccc).Corr4 = coLin{ii}(4);
                    S(ccc).Corr5 = coLin{ii}(5);
                    S(ccc).Corr6 = coLin{ii}(6);
                end
                
                nn = [P.root filesep fn filesep P.DestFold filesep P.bv.LogFileName];
                WriteDataToText(S, nn, 'a', '\t');

                badMotion = [];
                for ii = 1:length(coLin);
                    badMotion(ii) = any(coLin{ii}>P.cm.TaskCorrThresh);
                end
                badMotion = find(badMotion>0);

                if ~isempty(badMotion)
                    for ii = 1:length(badMotion)
                        %try
                        whRu(ii) = find(RunIndex==badMotion(ii));
                        disp(['Junking Runs ' num2str(RunIndex(whRu(ii))) ' (' P.cm.SourceList{whRu(ii)} ')' ' due to high task correlated motion.']);
                        %catch; keyboard; end
                        PPP(ii).message = ['Junking Runs ' num2str(RunIndex(whRu(ii))) ' (' P.cm.SourceList{whRu(ii)} ')' ' due to high task correlated motion.'];
                        S = [];
                        save(['NoGo_TaskCorrelatedMotion' sprintf('%0.2d',br(ii)) '.mat'], 'S');
                    end
                    try
                        tmp = load([P.root filesep fn filesep P.DestFold filesep P.BadRuns]);
                    catch
                        tmp = [];
                    end
                    if size(tmp,1)>size(tmp,2); tmp = tmp'; end
                    tmp = sort(unique([tmp RunIndex(whRu(ii))]));
                    dlmwrite([P.root filesep fn filesep P.DestFold filesep P.BadRuns],tmp,'delimiter','\t');
                    WriteDataToText(PPP, nn, 'a', '\t');
                    cd(P.root);
                    Create_Model(fn)
                    return;
                end
                S = []; S.a = ' ';
                WriteDataToText(S, nn, 'a', '\t');
                 
                %edit(nn)
            end
            
            
            if P.cm.DisableThresholdMasking == 1;
                %SPM.xM.T(:) = -Inf;  %% disable threshold masking
                SPM.xM.TH(:) = -Inf;  %% disable threshold masking
            else
                SPM.xM.I = 1;
            end
            
            if isfield(P.cm, 'ExplicitMask');
                SPM.xM.VM = P.cm.ExplicitMask;
            end
            %spm_print
            SPM = spm_spm(SPM);

            spm_get_defaults('stats.fmri.fmri_t',olddefs.stats.fmri.fmri_t); % Restore old timing
            spm_get_defaults('stats.fmri.fmri_t0',olddefs.stats.fmri.fmri_t0); % parameters

            
            P1.P.cm = P.cm;
            P1.P.cm.UserTime = UserTime;
            P = P1.P;
            save([P.root filesep fn filesep P.DestFold filesep pname '.mat'], 'P');
            
            cd(P.root);
        else
            cd(P.root);
            disp('Skipping Model Formation and Estimation');
        end
    end

    function Make_Cons(fn)
        fprintf('\n%s\n', 'Creating Contrast Images:');
        eval(['P = ' pname ';']);
        if isempty(P.nVols)
            ch = ReadInFile([P.root filesep fn filesep P.DestFold '/' P.List],'\t');
            tmp = dir_wfp([P.root filesep fn filesep P.bv.SourceFold filesep P.bv.SourcePrefix '*' P.bv.SourceSuffix]);
            for ii = 1:length(ch); if strcmpi(ch{ii},'NA'); P.nVols(ii) = NaN;  else; P.nVols(ii) = numel(spm_vol(tmp{contains(ch{ii},tmp)})); end; end;
        end
        
        if P.ct.do == 1;
            try
                cd([P.root filesep fn filesep P.cm.DestFold]);
            catch
                disp('No Analysis Folder Found!');
                return;
            end
            delete('ess*');
            delete('con*');
            delete('spmF*');
            delete('spmT*');
            P1 = load([P.root filesep fn filesep P.DestFold filesep pname '.mat']);
            
            x = load('SPM.mat');
            SPM = x.SPM;
            Con = [];
            cc = 0;

            for ii = 1:length(P.ct.con);
                %if ~isempty(contains('aschultz',{UserTime})); keyboard; end
                
                [vec L R] = CreateConVec(P.ct.con(ii).left,P.ct.con(ii).right,SPM,P.ct.con(ii).WeightWithin,P.ct.con(ii).BlockThresh);
                %disp([sum(L.counts) sum(R.counts)]);
                if isempty(P.ct.con(ii).right);
                    if sum(L.counts)<P.ct.MinEvents
                        nnn = [P.root filesep fn filesep P.DestFold filesep P.bv.LogFileName];
                        S(1).message = ' ';
                        S(2).message = ['Not Enough Events to Form the ' P.ct.con(ii).name ' Contrast.'];
                        WriteDataToText(S, nnn, 'a', '\t');
                        continue;
                    end
                elseif isempty(P.ct.con(ii).left);
                    if sum(R.counts)<P.ct.MinEvents
                        nnn = [P.root filesep fn filesep P.DestFold filesep P.bv.LogFileName];
                        S(1).message = ' ';
                        S(2).message = ['Not Enough Events to Form the ' P.ct.con(ii).name ' Contrast.'];
                        WriteDataToText(S, nnn, 'a', '\t');
                        continue;
                    end
                else
                    if sum(L.counts)<P.ct.MinEvents || sum(R.counts)<P.ct.MinEvents
                        nnn = [P.root filesep fn filesep P.DestFold filesep P.bv.LogFileName];
                        S(1).message = ' ';
                        S(2).message = ['Not Enough Events to Form the ' P.ct.con(ii).name ' Contrast.'];
                        WriteDataToText(S, nnn, 'a', '\t');
                        continue;
                    end
                end
                
                if isempty(vec)
                    continue
                end
                
                cc = cc+1;
                Con(cc).name = P.ct.con(ii).name;
                if size(vec,1)==1
                    Con(cc).STAT = 'T';
                else
                    Con(cc).STAT = 'F';
                end
                Con(cc).c = vec';
            end

            if length(Con)==0
               error('No Contrasts met the Minimum Event# criteria'); 
            end
            
            for ii = 1:length(Con)
                xCon(ii) = spm_FcUtil('Set',Con(ii).name,Con(ii).STAT,'c',Con(ii).c,SPM.xX.xKXs);
            end
            
            SPM.xCon = xCon;
            SPM = spm_contrasts(SPM,1:length(SPM.xCon));
            
            load SPM.mat
            for ii = 1:length(SPM.xCon)
                curname = SPM.xCon(ii).Vcon.fname(1:end-4);
                newname = SPM.xCon(ii).name;
                SPM.xCon(ii).Vcon.fname = [curname(1:4) newname '.img'];
                
                movefile([curname '.hdr'], [curname(1:4) newname '.hdr'])
                movefile([curname '.img'], [curname(1:4) newname '.img'])
                
                f = '0000';
                f2 = num2str(ii);
                f(end-length(f2)+1:end) = f2;
                
                if SPM.xCon(ii).STAT == 'T';
                    copyfile(['spmT_' f '.hdr'], ['spmT_' newname '.hdr']);
                    copyfile(['spmT_' f '.img'], ['spmT_' newname '.img']);
                end
                if SPM.xCon(ii).STAT == 'F';
                    copyfile(['spmF_' f '.hdr'], ['spmF_' newname '.hdr']);
                    copyfile(['spmF_' f '.img'], ['spmF_' newname '.img']);
                end
            end
            save SPM.mat SPM

            if P.ct.MoveContrasts == 1
                disp('Copying Contrast Images');
                
                
                %%% Remove any existing contrasts for this session;
                p = genpath([P.root filesep P.ct.GroupConFold filesep]);
                p = [':' p];
                ind = find(p==':');
                
                for ii = 1:length(ind)-1;
                    tmp = p(ind(ii)+1:ind(ii+1)-1);
                    tmp = [tmp filesep fn '*'];
                    try; ls(tmp); end
                    delete(tmp);
                end
                
                nPref = [];
                if (isfield(P.ct,'Normalize') && ~isempty(P.ct.Normalize))
                    D = [];
                    eval(['D = ' P.ct.Normalize ';']);
                    if ~(isfield(D.nn.rflags,'skipThis') && D.nn.rflags.skipThis==1)
                        error('Something is off with the parameters, and the normalization of contrast images.');
                    end
                    delete nn_*
                    df = [P.root filesep fn filesep D.nn.DestFold '/deformations.mat'];
                    nn = dir_wfp('*.img');
                    spm_write_sn(char(nn), df, D.nn.rflags);
                    nPref = D.nn.rflags.prefix;
                end
                
               
                for ii = 1:length(SPM.xCon)
                    if exist([P.root filesep P.ct.GroupConFold filesep SPM.xCon(ii).name])==0
                        mkdir([P.root filesep P.ct.GroupConFold filesep SPM.xCon(ii).name]);
                    end
                    
                    copyfile([nPref SPM.xCon(ii).Vcon.fname(1:end-3) 'img'], [P.root filesep P.ct.GroupConFold filesep Con(ii).name filesep fn SPM.xCon(ii).Vcon.fname(4:end-3) 'img']);
                    copyfile([nPref SPM.xCon(ii).Vcon.fname(1:end-3) 'hdr'], [P.root filesep P.ct.GroupConFold filesep Con(ii).name filesep fn SPM.xCon(ii).Vcon.fname(4:end-3) 'hdr']);
                end
            end
            
            
            P1.P.ct = P.ct;
            P1.P.ct.UserTime = UserTime;
            P = P1.P;
            
            save([P.root filesep fn filesep P.DestFold filesep pname '.mat'], 'P');
            
            cd(P.root);
        else
            cd(P.root);
            disp('Skipping Contrast Estimation');
        end
    end

    function persisText(txt,ii)
        % CUMDISP persistent disp
        % cumdisp; initializes persistent display
        % cumdisp(text); displays persistent text
        %
        if nargin==2
            if ii == 1
                fprintf('\n');
            end
        end
        persistent oldtxt;
        if nargin<1
            oldtxt='';
            fprintf(1,'\n');
        else
            fprintf(1,[repmat('\b',[1,length(oldtxt)]),txt]);
            oldtxt=sprintf(txt);
        end
    end

    function out = UserTime
        tmp = pwd;
        cd ~
        user = pwd;
        cd(tmp);
        
        ind = find(user == filesep);
        if ind(end)==numel(user);
            user = user(ind(end-1)+1:ind(end)-1);
        else
            user = user(ind(end)+1:end);
        end
        out = ['last run by ' user ' on ' datestr(clock)];
    end

end







