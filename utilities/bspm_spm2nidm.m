function bspm_spm2nidm()

command = sprintf('bet %s %s -f %2.1f',input, output, f);

curl -sLo -  https://github.com/incf-nidash/nidmresults-spm/archive/12.575ac2c.tar.gz | tar xzvf -
export exporter_path=`pwd`/nidmresults-spm-12.575ac2c/exporter
cp -p $SPMDIR/config/spm_cfg_results.m $SPMDIR/config/spm_cfg_results_ORIGINAL.m
cp -p $SPMDIR/config/spm_run_results.m $SPMDIR/config/spm_run_results_ORIGINAL.m
cp -p $exporter_path/spm_cfg_results.m $SPMDIR/config/spm_cfg_results.m
cp -p $exporter_path/spm_run_results.m $SPMDIR/config/spm_run_results.m
addpath('<PATH_TO_EXPORTER>')


matlabbatch{1}.spm.stats.results.spmmat = {'/Users/bobspunt/Documents/fmri/dog/_groupstats_/SURF1_ls6w2bs_Pmodby_None_WLS_100s_May_19_2016/OSTT_N17_PCTIN85_NOMASK_20160523/Dog_-_Scramble/SPM.mat'};
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = 'pdf';
matlabbatch{1}.spm.stats.results.write.none = 1;
matlabbatch{1}.spm.stats.results.export.nidm.subjects.subject = 1;
matlabbatch{1}.spm.stats.results.export.nidm.modality = 2;
matlabbatch{1}.spm.stats.results.export.nidm.refspace = 2;


for i = 1:length(cmd)
 [status, result] = system(cmd);
 end