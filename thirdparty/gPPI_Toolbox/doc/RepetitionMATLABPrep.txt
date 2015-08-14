%% Traditional PPI #1
create_sphere_image('SPM.mat',[0 -28 37],{'Tutorial_PCC'},6)
P.subject='Repetition_Tutorial';
P.directory=pwd;
P.VOI=[pwd filesep 'Tutorial_PCC_mask.nii'];
P.Region='PCC_traditional1'
P.Estimate=1;
P.contrast=0;
P.extract='eig';
P.Tasks={'N1'  'N2'  'F1'  'F2'};
P.Weights=[1 -1 -1 1];
P.analysis='psy';
P.method='trad';
P.CompContrasts=0;
P.Weighted=0;

%% Traditional PPI #2
create_sphere_image('SPM.mat',[0 -28 37],{'Tutorial_PCC'},6)
P.subject='Repetition_Tutorial';
P.directory=pwd;
P.VOI=[pwd filesep 'Tutorial_PCC_mask.nii'];
P.Region='PCC_traditional1'
P.Estimate=1;
P.contrast=0;
P.extract='eig';
P.Tasks={'N1'  'N2'};
P.Weights=[1 -1];
P.analysis='psy';
P.method='trad';
P.CompContrasts=0;
P.Weighted=0;

%% gPPI Example
create_sphere_image('SPM.mat',[0 -28 37],{'Tutorial_PCC'},6)
P.subject='Repetition_Tutorial';
P.directory=pwd;
P.VOI=[pwd filesep 'Tutorial_PCC_mask.nii'];
P.Region='PCC_gPPI_Repetition'
P.Estimate=1;
P.contrast=0;
P.extract='eig';
P.Tasks={'1' 'N1'  'N2'  'F1'  'F2'};
P.Weights=[];
P.analysis='psy';
P.method='cond';
P.CompContrasts=1;
P.Weighted=0;
P.ConcatR=0;
P.preservevarcorr=0;
P.Contrasts(1).left={'N1'};
P.Contrasts(1).right={'N2'};
P.Contrasts(1).STAT='T';
P.Contrasts(1).Weighted=0;
P.Contrasts(1).MinEvents=5;
P.Contrasts(1).name='N1_minus_N2';
P.Contrasts(2).left={'F1'};
P.Contrasts(2).right={'F2'};
P.Contrasts(2).STAT='T';
P.Contrasts(2).Weighted=0;
P.Contrasts(2).MinEvents=5;
P.Contrasts(2).name='F1_minus_F2';
P.Contrasts(3).left={'N1'  'F2'};
P.Contrasts(3).right={'N2'  'F1'};
P.Contrasts(3).STAT='T';
P.Contrasts(3).Weighted=0;
P.Contrasts(3).MinEvents=5;
P.Contrasts(3).name='N1_F2_minus_N2_F1';
