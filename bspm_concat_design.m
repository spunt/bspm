
%     'Sn(1) preWhy*bf(1)'
%     'Sn(1) preWhyxt-1^1*bf(1)'
%     'Sn(1) preWhyxITI_Duration^1*bf(1)'
%     'Sn(1) preWhyxRT^1*bf(1)'
%     'Sn(1) preHow*bf(1)'
%     'Sn(1) preHowxt-1^1*bf(1)'
%     'Sn(1) preHowxITI_Duration^1*bf(1)'
%     'Sn(1) preHowxRT^1*bf(1)'
%     'Sn(1) preMath*bf(1)'
%     'Sn(1) preMathxt-1^1*bf(1)'
%     'Sn(1) preMathxITI_Duration^1*bf(1)'
%     'Sn(1) preMathxRT^1*bf(1)'
clear all
studydir = '/Users/bobspunt/Documents/fmri/mbc';
analysis = 'pretrialpm_ve_ctrlyesno_itidur_tmin1_sdcut3';
subdir = files([studydir filesep 'MBC*']);
cond = {'preWhy' 'preHow' 'preMath'};
pm = {'t-1' 'ITI' 'RT'};
count = 1;
for p = 1:length(pm)
    for c = 1:length(cond)
        pmnam{count} = [cond{c} 'x' pm{p}];
        count=count+1;
        condnam{c} = [cond{c} '**bf'];
    end;
end;
allcond = [condnam pmnam];
nicecondnam = [cond pmnam];

for s = 1:length(subdir)
    
    spmmat = files([subdir{s} filesep 'analysis' filesep analysis filesep 'SPM*mat']);

    % get X pmmatrix and names
    load(spmmat{1})
    xmatrix = SPM.xX.xKXs.X; % filtered and whitened design matrix
    xname = SPM.xX.name; % regressor names
    for i = 1:length(allcond)

        condidx =find(~cellfun('isempty',regexp(xname,allcond{i})));
        tmpx = xmatrix(:,condidx);
        ccx(:,i) = sum(tmpx')';

    end
    
    outdir = [subdir{s} filesep 'behav'];
    currentfilename=[outdir filesep 'cc_X_ve.txt'];
    save(currentfilename,'ccx','-ascii');
    clear ccx SPM
    save([outdir filesep 'cc_X_ve_condnames.mat'], 'nicecondnam');
    
end


