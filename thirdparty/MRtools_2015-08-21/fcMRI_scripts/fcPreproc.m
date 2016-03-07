function nfn = fcPreproc(fn,pname,ord)
%%% All configuration options should be specified in the Param file.
%%%
%%% FN: is a 4D nifti file containing the resting state scans.
%%%
%%% PNAME: is the name of the Params file that you want to use.  (note do
%%% not include .m at the end of the param file name).
%%%
%%% ORD:  This is an optional input that you can use to respcify the order
%%% of processing without changing the order in the Params.m file. THis
%%% will specify the order of execution of P.PO.
%%%
%%% Written by Aaron Schultz, May 5th, 2010 (aschultz@martinos.org)
%%%
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%%
%%% Last Updated Dec. 11 2012;
%%%
%%% Copyright (C) 2012,  Aaron P. Schultz
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

if ischar(fn)
    fn = {fn};
end

nfn = fn;
if nargin<2
    pname = 'Params';
    if exist(pname)
        eval(['P = ' pname ';']);
        pname = 'Params';
    else
        disp('Warning! No Params.m file was found.  A new Params.m file has been written in the current directory.  Please check the settings in this file, and then rerun the script.');
        return
    end
else
    eval(['P = ' pname ';']);
end

if exist([pname '.mat'])==0;
    save([pname '.mat'],'P');
end

if nargin<3
   ord = 1:length(P.PO); 
end

for qq = 1:length(ord)
   eval(P.PO{ord(qq)});
end

try
    P1 = load([pname '.mat']);
    P1.P.UserTime = UserTime;
    P = P1.P;
catch
    cd(P.root);
    disp('Something Went Wrong With This Analysis!');
    fprintf('\n\n\n\n');
    return
end

save([pname '.mat'],'P');

    function nfn = drop_vols(fn)
        for zz = 1:numel(fn);
            eval(['P = ' pname ';']);
            P1 = load([pname '.mat']);
            if P.dv.do == 1
                fprintf('\n%s\n', ['DROPPING ' num2str(length(P.dv.dropVols)) ' VOLUMES']);
                [M V] = openIMG(fn{zz});
                vec = setdiff(1:length(V), P.dv.dropVols);
                c = 0;
                for ii = vec
                    c = c+1;
                    persisText(['Writing Volume #' num2str(c)],ii);
                    V = spm_vol( [fn{zz} ',' num2str(ii)]);
                    V.fname = [P.dv.prefix fn{zz}];
                    V.n = [c 1];
                    spm_write_vol(V,M(:,:,:,ii));
                end
                nfn{zz} = [P.dv.prefix fn{zz}];
                
                P1.P.dv = P.dv;
                P1.P.dv.UserTime = UserTime;
                P = P1.P;
                save([pname '.mat'], 'P');
            else
                nfn{zz} = fn;
            end
        end
    end

    function nfn = slice_time(fn)      
        fprintf('\n%s\n', 'PERFORMING SLICE TIME CORRECTION:');
        eval(['P = ' pname ';']);
        P1 = load([pname '.mat']);
        nfn = fn;
        
        if P.st.do == 1;
            for zz = 1:numel(fn)
                spm_slice_timing(fn{zz}, P.st.sliceorder, P.st.refslice, P.st.timing, P.st.prefix);
                nfn{zz} = [P.st.prefix fn{zz}];
            end
            P1.P.st = P.st;
            P1.P.st.UserTime = UserTime;
            P = P1.P;
            save([pname '.mat'], 'P');
        end
    end

    function nfn = realign(fn)
        fprintf('\n%s\n', 'REALIGNING IMAGES:');
        eval(['P = ' pname ';']);
        P1 = load([pname '.mat']);
        
        if P.rr.do == 1;
            %keyboard;
            if P.rr.RealignmentType == 1;
                inria_realign(fn,P.rr.RealignParsI);
            else
                spm_realign(fn,P.rr.RealignParsS);
            end
            spm_reslice(fn,P.rr.ReslicePars);
            
            if P.rr.ReslicePars.which==2 || P.rr.ReslicePars.which==1
                nfn = [P.rr.ReslicePars.prefix fn];
            else
                nfn = fn;
            end
            
            
            P1.P.rr = P.rr;
            P = P1.P;
            save([pname '.mat'], 'P');             
        else
            nfn = fn;
        end
    end

    function nfn = normalize(fn) 
        fprintf('\n%s\n', 'NORMALIZING IMAGES:');
        eval(['P = ' pname ';']);
        P1 = load([pname '.mat']);
        nfn = fn;
        if P.nn.do == 1
            nn = dir([P.nn.source]);
            if numel(nn)>1; disp('There is more than one mean image in this folder'); return; end;
            source = [nn.name];
            template = P.nn.template;
            
            spm_normalise(template,source,'deformations.mat','','',P.nn.NormPars);
            
            for zz = 1:numel(fn)
                spm_write_sn(fn{zz}, 'deformations.mat', P.nn.rflags);
                nfn{zz} = [P.nn.rflags.prefix fn{zz}];
                Data_QC(nfn{zz});
            end
            
            P1.P.nn = P.nn;
            P1.P.nn.UserTime = UserTime;
            P = P1.P;
            save([pname '.mat'], 'P');
        end
    end

    function nfn = smooth(fn)
        fprintf('\n%s\n', 'SMOOTHING IMAGES:');
        eval(['P = ' pname ';']);
        P1 = load([pname '.mat']);
        nfn = fn;
        
        if P.ss.do == 1
            for zz = 1:numel(fn)
                spm_smooth(fn{zz},[P.ss.prefix fn{zz}],P.ss.kernel,16);
                nfn{zz} = [P.ss.prefix fn{zz}];
                
                if isfield(P, 'bv');
                    if P.bv.do == 1;
                        Bad_Vols(nfn{zz});
                    end
                end
            end
            P1.P.ss = P.ss;
            P1.P.ss.UserTime = UserTime;
            P = P1.P;
            save([pname '.mat'], 'P');
        end
    end

    function Bad_Vols(fn)
        fprintf('\n%s\n', 'Examing Image Files and Creating Bad Volume Regressors:');
        eval(['P = ' pname ';']);
        %keyboard;
        if P.bv.do == 1;
            P1 = load([pname '.mat']);

            nnn = P.bv.LogFileName;
            %%% Collect Global SNR Values, for writing to analysis log and
            %%% screening for bad runs.
            if P.bv.CheckGlobalSNR
                %%% Get list of file names that have the SNR data
                [ll l2] = dir_wfp(['SNR_Images' filesep P.bv.SNR_prefix '*' fn(4:end-4) '*' P.bv.SNR_suffix]);

                if(numel(ll)>1)
                    error('Too Many SNR data files!');
                end
                %%% Collect the SNR infro and write it to the Analysis Log.
                
                snr = dlmread(ll{1},'\t',1,1);
                
                disp(snr(1,:))
                S = [];
                S = {'From Run:' 'Mean:' 'SD:' 'SNR:'; l2{1} snr(1,2) snr(1,2) snr(1,3); ' ' ' ' ' ' ' '};
                
                WriteDataToText(S, nnn, 'w', '\t');

                %%% Screen SNR values for abnormally low values.
                if snr(3) < P.bv.SNRthresh || isnan(snr(3));
                    disp(['Junking Resting Run due to poor SNR values.']);
                    S =[];
                    S.a1 = ['Junking Resting Run due to poor SNR values.'];
                    WriteDataToText(S, nnn, 'a', '\t');
                    save NoGo_Bad_SNR.mat nnn;
                end

                S = [];
                WriteDataToText({' '}, nnn, 'a', '\t');
            end
            
            [nn2 nn2_1] = dir_wfp([P.bv.RealignPrefix '*' fn(7:end-4) '*.txt']);
            if numel(nn2)>1
                error('Too Many Realignment Files!');
            end
            %%% Look for bad volumes using both the global signal and the
            %%% movement parameters
            
            disp(['Examining Resting Run']);
            
            tmp = fn(4:end);
            
            MM = openIMG(tmp);
            GM = [];
            for jj = 1:size(MM,4);
                M = MM(:,:,:,jj);
                indd = find(M>((mean(M(1:end)))/8));
                GM(jj) = mean(M(indd));
            end
            GM1 = diff(GM);
            ind1 = find(abs(GM1)>(std(GM1)*P.bv.GlobalSigThresh))+1;
            %keyboard;
            dat = load(nn2{1});
            
            mv = mean(sqrt(sum(diff(dat(:,1:3)).^2,2)));
            WriteDataToText({['Mean Movement = ' sprintf('%0.4f',mv)]}, nnn, 'a', '\t');
            if mv > P.bv.MeanMovementThresh
                disp('Junking Rest Run due to excessive mean movement.');
                S = {''; ['Junking Rest Run due to excessive mean movement. MV = ' sprintf('%0.3f',mv)]};
                WriteDataToText(S, nnn, 'a', '\t');
                save NoGo_Too_Much_Movement.mat mv;
            end
            WriteDataToText({' '}, nnn, 'a', '\t');
            
            tm = DistMat(dat(:,1:3),0);
            tm = max(tm(:));
            tr = DistMat(dat(:,4:6),0);
            tr = max(tr(:));
            
            
            S = {''; ['Total Movement = ' num2str(round(1000*(tm))/1000) ' milimeters']; ['Total Rotation = ' num2str(round(1000*(tr))/1000) ' degrees']};
            if tm > P.bv.TotalMovementThresh
                disp('Junking Rest Run due to too much overall movement.');
                S{end+1} = 'Junking Rest Run due to too much overall movement.';
                save NoGo_Bad_Movement.mat nnn;
            end
            if tr > P.bv.TotalRotationThresh
                disp(['Junking Rest Run due to too much overall rotation.']);
                S{end+1} = ['Junking Rest Run due to too much overall rotation.'];
                save NoGo_Bad_Rotation.mat nnn;
            end
            S{end+1} = '';
            disp(S{2});
            disp(S{3});
            WriteDataToText(S, nnn, 'a', '\t');
            
            tmp = diff(dat).^2;
            pos = sqrt(sum(tmp(:,1:3),2));
            ori = sqrt(sum(tmp(:,4:6),2)) * (360/(2*pi));
            ind2 = find(pos>P.bv.MoveThresh)+1;
            ind3 = find(ori>P.bv.RotThresh)+1;
            
            ind = unique([ind1'; ind2; ind3]);
            R =  zeros(size(MM,4),numel(ind));
            for qq = 1:length(ind);
                R(ind(qq),qq) = 1;
            end
            
            
            S = [];
            S{1} = ['Negating ' num2str(length(ind)) ' volumes; With local indices of: '  num2str(ind')];
            WriteDataToText(S, nnn, 'a', '\t');
            
            if numel(ind)>P.bv.BadVolThresh
                disp(['Junking Rest Run due to too many bad volumes.']);
                S = [];
                S.a1 = ['Junking Rest Run due to too many bad volumes.'];
                WriteDataToText(S, nnn, 'a', '\t');
                save NoGo_Too_Many_Bad_Vols.mat nnn;
            end
                
            S = [];
            WriteDataToText({' '}, nnn, 'a', '\t');
            type(nnn);
            
            save([P.bv.BadVolRegsName(1:end-4) fn(end-4) '.mat'], 'R')                        

            %%% Finish up.
            P1.P.bv = P.bv;
            P1.P.bv.UserTime = UserTime;
            P = P1.P;
            save([pname '.mat'], 'P');
            
        else
            fprintf('\n%s\n', 'Skipping run/volume assessment.');
        end
    end

    function nfn = filter_data(fn)        
        fprintf('\n%s\n', 'FILTERING DATA:');
        eval(['P = ' pname ';']);
        P1 = load([pname '.mat']);
        
        nfn = fn;
        
        if P.ff.do == 1
            for zz = 1:numel(fn)
                fprintf('\n%s\n', 'Reading in image data.');
                [M V] = openIMG(fn{zz});
                
                fprintf('\n%s\n', 'Reshaping Data.');
                ss = size(M);
                M2 = double(reshape(M, prod(ss(1:3)),ss(4))');
                
                ss2 = size(M2);
                
                if P.ff.Detrend == 1
                    fprintf('\n%s\n', 'Detrending Data.');
                    M2 = detrend(M2);
                end
                
                fprintf('\n%s\n', 'Filtering Data.');
                if P.ff.LowCut == 0
                    M3 = ft_preproc_highpassfilter(M2',1/P.TR,P.ff.HighCut,P.ff.FilterOrder,'but','twopass' ,'no')';
                else
                    M3 = ft_preproc_bandpassfilter(M2',1/P.TR,[P.ff.LowCut P.ff.HighCut],P.ff.FilterOrder,'but','twopass' ,'no')';
                end
                
                M5 = reshape(M3', ss(1), ss(2), ss(3), ss(4));
                
                fprintf('\n%s\n', 'Writing Out Filtered Data.');
                for ii = 1:length(V);
                    persisText(['Writing Volume #' num2str(ii)],ii);
                    V = spm_vol([fn{zz} ',' num2str(ii)]);
                    V.fname = [P.ff.prefix fn{zz}];
                    V.n = [ii 1];
                    V.dt = [64 0];
                    spm_write_vol(V,M5(:,:,:,ii));
                end
                nfn{zz} = [P.ff.prefix fn{zz}];
            end
            fprintf('\n');
            P1.P.ff = P.ff;
            P1.P.ff.UserTime = UserTime;
            P = P1.P;
            save([pname '.mat'], 'P');
        end
    end

    function nfn = motion_regress(fn)
        fprintf('\n%s\n', 'REMOVING MOTION VARIANCE:');
        eval(['P = ' pname ';']);
        P1 = load([pname '.mat']);

        nfn = fn;
        
        if P.mot.do == 1
            fprintf('\n%s\n', 'Loading Image Data.');
            for zz = 1:numel(fn);
                %%% Load in the actual Resting State Data                
                ti = find(fn{zz}=='_');
                nn = dir([P.bv.RealignPrefix fn{zz}(ti(2)+1:end-4) '.txt']);
                [dat V]  = openIMG(fn{zz});
                
                try
                    fprintf('\n%s\n', 'Loading Motion Regressors and Computing First Derivitive.');
                    %%% Load in motion regressors and compute the derivitive
                    mot1 = load(nn(1).name);
                    if P.mot.MotionDeriv == 1;
                        mot = [mot1 [0 0 0 0 0 0; diff(mot1)]];
                    end
                    
                    %%% These are the Regressors
                    Regs = [ones(size(dat,4),1) mot];
                    fprintf('\n%s\n', [num2str(size(Regs,2)-1) ' Regressors are being factored out.']);
                    
                    fprintf('\n%s\n', 'Reshaping Data');
                    ss = size(dat);
                    D = double(reshape(dat, prod(ss(1:3)),ss(4))');
                    ss2 = size(D);
                    
                    fprintf('\n%s\n', 'Regressing Out Nuisance Variables.');
                    %%% Run the Regression analysis to remove variability attributable to the
                    %%% regressors.
                    
                    b = pinv(Regs)*D;
                    pred = Regs*b;
                    M3= D-pred;
                end
                M5 = reshape(M3', ss(1), ss(2), ss(3), ss(4));
                
                fprintf('\n%s\n', 'Writing Out Data.');
                V = spm_vol(fn{zz});
                
                for ii = 1:length(V);
                    persisText(['Writing Volume #' num2str(ii)],ii);
                    V = spm_vol([fn{zz} ',' num2str(ii)]);
                    V.fname = [P.mot.prefix fn{zz}];
                    V.n = [ii 1];
                    V.dt = [64 0];
                    spm_write_vol(V,M5(:,:,:,ii));
                end
                nfn{zz} = [P.mot.prefix fn{zz}];
            end
            fprintf('\n');
            P1.P.mot = P.mot;
            P1.P.mot.UserTime = UserTime;
            P = P1.P;
            save([pname '.mat'], 'P');
        end
    end

    function nfn = nuisance_regress(fn)
        fprintf('\n%s\n', 'REMOVING NUISANCE VARIANCE:');
        
        eval(['P = ' pname ';']);
        P1 = load([pname '.mat']);
        
        nfn = fn;
        
        if P.res.do == 1
            for zz = 1:numel(fn);
                fprintf('\n%s\n', 'Loading Image Data.');
                %%% Load in the actual Resting State Data
                [dat V]  = openIMG(fn{zz});
                
                if P.res.Motion
                    fprintf('\n%s\n', 'Loading Motion Regressors and Computing First Derivitive.');
                    %%% Load in motion regressors and compute the derivitive
                    try
                        if isfield(P.mot,'file')
                            nn = dir(P.mot.file);
                        else
                            nn = dir('realignment_params_*.txt');
                        end
                        mot = load(nn(1).name);
                    catch
                        [t,sts] = spm_select(1,'any','could not find motion regressors. please selct file with all regressors');
                        mot = load(t);
                    end
                    if P.res.MotionDeriv == 1;
                        mot = [mot [0 0 0 0 0 0; diff(mot)]];
                    end
                    
                else
                    mot = [];
                end
                
                if ~isempty(P.res.Masks)
                    fprintf('\n%s\n', 'Gather and Compute Physiological Regressors');
                    
                    %%% Get the white matter, ventricle and whole brain masks.
                    clear ind
                    for ii = 1:length(P.res.Masks);
                        RegLabs{ii} = P.res.Masks{ii}(1:end-4);
                        v = spm_vol(P.res.Masks{ii});
                        Nv = resizeVol(v,V(1));
                        %Nv = SliceAndDice(v,V(1),V(1),V(1),0,[]);
                        ind{ii} = find(isnan(Nv)==0);
                    end
                    
                    %%% compute average signal series for each mask
                    for ii = 1:size(dat,4)
                        tmp = dat(:,:,:,ii);
                        for jj = 1:length(ind)
                            BS(ii,jj) = mean(tmp(ind{jj}));
                        end
                    end
                    WriteDataToText(RegLabs, ['PhysioRegressors_' fn{zz}(1:end-4) '.txt'], 'w', '\t', 0, 4);
                    WriteDataToText(BS, ['PhysioRegressors_' fn{zz}(1:end-4) '.txt'], 'a', '\t', 0, 4);
                    
                    if P.res.MaskDeriv == 1
                        BS = [BS [zeros(1,size(BS,2)); diff(BS)]];
                    end
                else
                    BS = [];
                end
                
                %%% These are the Regressors
                %keyboard
                Regs = [ones(size(dat,4),1) mot BS];
                fprintf('\n%s\n', [num2str(size(Regs,2)-1) ' Regressors are being factored out.']);
                
                fprintf('\n%s\n', 'Reshaping Data');
                ss = size(dat);
                D = double(reshape(dat, prod(ss(1:3)),ss(4))');
                ss2 = size(D);
                
                
                fprintf('\n%s\n', 'Regressing Out Nuisance Variables.');
                %%% Run the Regression analysis to remove variability attributable to the
                %%% regressors.
                
                b = pinv(Regs)*D;
                pred = Regs*b;
                M3= D-pred;
                
                M5 = reshape(M3', ss(1), ss(2), ss(3), ss(4));
                
                fprintf('\n%s\n', 'Writing Out Data.');
                V = spm_vol(fn{zz});
                
                for ii = 1:length(V);
                    persisText(['Writing Volume #' num2str(ii)],ii);
                    V = spm_vol([fn{zz} ',' num2str(ii)]);
                    V.fname = [P.res.prefix fn{zz}];
                    V.n = [ii 1];
                    V.dt = [64 0];
                    spm_write_vol(V,M5(:,:,:,ii));
                end
                nfn{zz} = [P.res.prefix fn{zz}];
                
                
            end
            fprintf('\n');
            P1.P.res = P.res;
            P1.P.res.UserTime = UserTime;
            P = P1.P;
            save([pname '.mat'], 'P');
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







