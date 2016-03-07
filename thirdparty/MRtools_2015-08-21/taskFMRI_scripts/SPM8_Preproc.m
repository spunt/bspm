function SPM8_Preproc(fn,pname,ord);
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

nfn = fn;
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
    eval(['P = ' pname ';']);
end
cd([P.root filesep fn filesep P.DestFold]);
% P.MatDir

if exist([pname '.mat'])==0;
    save([pname '.mat'],'P');
end

if nargin<3
   ord = 1:length(P.PO); 
end

cd(P.root);
for qq = 1:length(ord)
    eval(P.PO{ord(qq)});
end

P1 = load([P.root filesep fn filesep P.DestFold filesep pname '.mat']);
P1.P.UserTime = UserTime;
P1.P.TR = P.TR;
P1.P.root = P.root;
P1.P.DestFold = P.DestFold;
P1.P.List = P.List;
P1.P.Runs = P.Runs;
P = P1.P;

save([P.root filesep fn filesep P.DestFold filesep pname '.mat'],'P');

    function nfn = drop_vols(fn)
        eval(['P = ' pname ';']);
        if P.dv.do == 1
            
            fprintf('\n%s\n', ['DROPPING ' num2str(length(P.dv.dropVols)) ' VOLUMES']);
            P1 = load([P.root filesep fn filesep P.DestFold filesep pname '.mat']);
            warning off
            mkdir([P.root filesep fn filesep P.dv.DestFold]);
            warning on;
            cd([P.root filesep fn filesep P.dv.SourceFold]);
            
            list = ReadInFile([P.root filesep fn filesep P.DestFold filesep P.List],'\t');
            list = list(:,1);
            
            if ~isempty(P.Runs);
                if numel(list)~=numel(P.Runs); error('number of RunIDs do not match number of specified runs!  Fix the RunIDs and then try rerunning.'); end
            end
            
            ind = setdiff(1:length(list), contains('NA',list));
            list = list(ind);
            
            [trash nn] = dir_wfp([P.dv.SourcePrefix '*' P.dv.SourceSuffix]);
            ind = [];
         
            for ii = 1:length(list); ind(ii) = contains(list{ii}, nn); end
            nn = nn(ind);
            
            for ii = 1:numel(nn);
                [M V] = openIMG(nn{ii});
                vec = setdiff(1:length(V), P.dv.dropVols);
                c = 0;
                for jj = vec
                    c = c+1;
                    persisText(['Writing Volume #' num2str(c)],jj);
                    V = spm_vol( [nn{ii} ',' num2str(jj)]);
                    V.fname = [P.dv.prefix nn{ii}];
                    V.n = [c 1];
                    spm_write_vol(V,M(:,:,:,jj));
                end
            end
            
            if ~strcmp(pwd, [P.root filesep fn filesep P.dv.DestFold])
                movefile([P.dv.prefix '*'], [P.root filesep fn filesep P.dv.DestFold]);
            end
            
            cd(P.root);
            nfn = fn;
            P1.P.dv = P.dv;
            P = P1.P;
            save([P.root filesep fn filesep P.DestFold filesep pname '.mat'], 'P');
        else
            fprintf('\n%s\n', 'SKIPPING DROP VOLUMES');
            nfn = fn;
        end
    end

    function nfn = slice_time(fn)
        
        fprintf('\n%s\n', 'PERFORMING SLICE TIME CORRECTION:');
        eval(['P = ' pname ';']);
        
        if P.st.do == 1;
            P1 = load([P.root filesep fn filesep P.DestFold filesep pname '.mat']);
            warning off
            mkdir([P.root filesep fn filesep P.st.DestFold]);
            warning on
            cd([P.root filesep fn filesep P.st.SourceFold]);
            
            list = ReadInFile([P.root filesep fn filesep P.DestFold filesep P.List],'\t');
            list = list(:,1);

            if ~isempty(P.Runs)
                if numel(list)~=numel(P.Runs); error('number of RunIDs do not match number of specified runs!  Fix the RunIDs and then try rerunning.'); end
            end
            
            ind = setdiff(1:length(list), contains('NA',list));
            list = list(ind);
            
            nn = dir_wfp([P.st.SourcePrefix '*' P.st.SourceSuffix]);
            ind = [];
            %keyboard;
            for ii = 1:length(list); ind(ii) = contains(list{ii}, nn); end
            nn = nn(ind);
            
            for ii = 1:length(nn)
                spm_slice_timing(nn{ii}, P.st.sliceorder, P.st.refslice, P.st.timing, P.st.prefix);
            end
            
            if ~strcmp(pwd, [P.root filesep fn filesep P.st.DestFold])
                movefile([P.st.prefix '*'], [P.root filesep fn filesep P.st.DestFold]);
            end
            
            cd(P.root);
            nfn = fn;
            P1.P.st = P.st;
            P = P1.P;
            save([P.root filesep fn filesep P.DestFold filesep pname '.mat'], 'P');
            
            
        else
            fprintf('\n%s\n', 'SKIPPING SLICE TIMING');
            nfn = fn;
        end
    end

    function nfn = realign(fn)
        fprintf('\n%s\n', 'REALIGNING IMAGES:');
        eval(['P = ' pname ';']);
        
        if P.rr.do == 1;
            %keyboard;
            P1 = load([P.root filesep fn filesep P.DestFold filesep pname '.mat']);
            cd([P.root filesep fn filesep P.rr.SourceFold]);
            
            list = ReadInFile([P.root filesep fn filesep P.DestFold filesep P.List],'\t');
            list = list(:,1);
            
            if ~isempty(P.Runs)
                if numel(list)~=numel(P.Runs); error('number of RunIDs do not match number of specified runs!  Fix the RunIDs and then try rerunning.'); end
            end
            
            ind = setdiff(1:length(list), contains('NA',list));
            list = list(ind);
                        
            nn = dir_wfp([P.rr.SourcePrefix '*' P.rr.SourceSuffix]);
            ind = [];

            for ii = 1:length(list); ind(ii) = contains(list{ii}, nn); end
            nn = nn(ind);
            
            if P.rr.RealignSeparate == 1;
                for ii = 1:length(nn);
                    if P.rr.RealignmentType == 1;
                        inria_realign(nn{ii},P.rr.RealignParsI);
                    else
                        spm_realign(nn{ii},P.rr.RealignParsS);
                    end
                    spm_reslice(nn{ii},P.rr.ReslicePars);
                    %movefile('mean*', [P.root filesep fn filesep P.rr.DestFold]);
                end
            end
            if P.rr.RealignSeparate == 2;
                if P.rr.RealignmentType == 1;
                    inria_realign(nn,P.rr.RealignParsI);
                else
                    %disp('here?');
                    spm_realign(nn,P.rr.RealignParsS);
                end                
                spm_reslice(nn,P.rr.ReslicePars);
            end
            
            if P.rr.ReslicePars.which==2;
                mkdir([P.root filesep fn filesep P.rr.DestFold]);
                movefile([P.rr.ReslicePars.prefix '*'], [P.root filesep fn filesep P.rr.DestFold]);
                
                if P.rr.RealignSeparate == 1
                tmp = dir_wfp('mean*.nii');
                for ii = 1:numel(tmp);
                    movefile(tmp{ii}, [P.root filesep fn filesep P.rr.DestFold '/']);
                end
                elseif P.rr.RealignSeparate == 2
                    movefile(['mean*'],[P.root filesep fn filesep P.rr.DestFold]);
                end
                
                if P.rr.RealignmentType == 1
                    movefile(['realignment_*'],[P.root filesep fn filesep P.rr.DestFold]);
                elseif P.rr.RealignmentType == 2
                    movefile(['rp*'],[P.root filesep fn filesep P.rr.DestFold]);
                end
            end
            
            delete SPMtmp.mat
            nfn = fn;
            P1.P.rr = P.rr;
            P = P1.P;
            save([P.root filesep fn filesep P.DestFold filesep pname '.mat'], 'P');
             
            cd(P.root);
        else
            nfn = fn;
        end
    end

    function nfn = normalize(fn) 
        fprintf('\n%s\n', 'NORMALIZING IMAGES:');
        eval(['P = ' pname ';']);
        if P.nn.do == 1
            P1 = load([P.root filesep fn filesep P.DestFold filesep pname '.mat']);
            cd([P.root filesep fn filesep P.nn.SourceFold]);     

            list = ReadInFile([P.root filesep fn filesep P.DestFold filesep P.List],'\t');
            list = list(:,1);

            if ~isempty(P.Runs)
                if numel(list)~=numel(P.Runs); error('number of RunIDs do not match number of specified runs!  Fix the RunIDs and then try rerunning.'); end
            end
            
            ind = setdiff(1:length(list), contains('NA',list));
            list = list(ind); 

            nn1 = dir_wfp([P.nn.SourcePrefix '*' P.nn.SourceSuffix]);

            ind = [];

            for ii = 1:length(list); ind(ii) = contains(list{ii}, nn1); end
            nn1 = nn1(ind);

            template = P.nn.template;
            nn2 = dir_wfp(P.nn.source);
            ind = [];
            
            if P.rr.RealignmentType == 1
                if P.rr.RealignSeparate ==1
                    for ii = 1:length(list); ind(ii) = contains(list{ii}, nn2); end
                    nn2 = nn2(ind);
                    for ii = 1:length(nn2)
                        source = nn2{ii};
                        spm_normalise(template,source,['deformations_' list{ii} '.mat'],'','',P.nn.NormPars);
                        spm_write_sn(nn1{ii}, ['deformations_' list{ii} '.mat'], P.nn.rflags);
                    end
                else
                    source = nn2{1};
                    spm_normalise(template,source,'deformations.mat','','',P.nn.NormPars);
                    for ii = 1:length(nn1)
                        spm_write_sn(nn1{ii}, 'deformations.mat', P.nn.rflags);
                    end
                end
            end
            if P.rr.RealignmentType == 2
                source = nn2{1};
                spm_normalise(template,source,'deformations.mat','','',P.nn.NormPars);
                if ~isfield(P.nn.rflags,'skipThis') || P.nn.rflags.skipThis==0
                    spm_write_sn(char(nn1), 'deformations.mat', P.nn.rflags);
                end
            end
            
            if ~strcmp(pwd, [P.root filesep fn filesep P.nn.DestFold])
                mkdir([P.root filesep fn filesep P.nn.DestFold]);
                try
                    movefile([P.nn.rflags.prefix '*'],[P.root filesep fn filesep P.nn.DestFold]);
                end
                movefile('deformations*',[P.root filesep fn filesep P.nn.DestFold]);
                
            end
            %%% Should add in move file commands here and the cd before the
            %%% SNR measurements
            
            if ~isfield(P.nn.rflags,'skipThis') || P.nn.rflags.skipThis==0
                cd([P.root filesep fn filesep P.nn.DestFold]);
                nn = dir_wfp([P.nn.rflags.prefix '*.nii']);
                for ii = 1:length(nn)
                    Data_QC(nn{ii});
                end
            else
                cd([P.root filesep fn filesep P.rr.DestFold]);
                for ii = 1:length(nn1)
                    Data_QC(nn1{ii});
                end
                cd('SNR_Images');
                nn = dir_wfp('*.nii');
                def = dir_wfp([P.root filesep fn filesep P.nn.DestFold '/deformations.mat']);
                if numel(def)==1
                    for ii = 1:numel(nn);
                        spm_write_sn(nn{ii}, def{1}, P.nn.rflags);
                    end
                    mkdir('NativeSpace');
                    movefile([P.rr.ReslicePars.prefix '*.nii'],'NativeSpace/')
                    system(['mv ' pwd ' ' P.root filesep fn filesep P.nn.DestFold '/']);
                    cd([P.root filesep fn filesep P.nn.DestFold]);
                end
            end
           
            P1.P.nn = P.nn;
            P1.root = P.root;
            P = P1.P;
            save([P.root filesep fn filesep P.DestFold filesep pname '.mat'], 'P');
             
            cd(P.root);
            
            nfn = fn;
        else
            nfn = fn;
        end
    end

    function nfn = smooth(fn)
        fprintf('\n%s\n', 'SMOOTHING IMAGES:');
        eval(['P = ' pname ';']);

        if P.ss.do == 1
            P1 = load([P.root filesep fn filesep P.DestFold filesep pname '.mat']);
            cd([P.root filesep fn filesep P.ss.SourceFold]);
            
            list = ReadInFile([P.root filesep fn filesep P.DestFold filesep P.List],'\t');
            list = list(:,1);
            
            if ~isempty(P.Runs)
                if numel(list)~=numel(P.Runs); error('number of RunIDs do not match number of specified runs!  Fix the RunIDs and then try rerunning.'); end
            end
            
            ind = setdiff(1:length(list), contains('NA',list));
            list = list(ind);
            
            [trash nn] = dir_wfp([P.ss.SourcePrefix '*' P.ss.SourceSuffix]);
            ind = [];
            for ii = 1:length(list); ind(ii) = contains(list{ii}, nn); end
            nn = nn(ind);
            
            
            
            for ii = 1:length(nn)
                spm_smooth(nn{ii},[P.ss.prefix nn{ii}],P.ss.kernel,0);
            end
            
            if ~strcmp(pwd, [P.root filesep fn filesep P.ss.DestFold])
                mkdir([P.root filesep fn filesep P.ss.DestFold]);
                movefile('ss*',[P.root filesep fn filesep P.ss.DestFold]);
            end
            
            P1.P.ss = P.ss;
            P = P1.P;
            save([P.root filesep fn filesep P.DestFold filesep pname '.mat'], 'P');
             
            cd(P.root);
            
            nfn = fn;
        else
            nfn = fn;
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

end







