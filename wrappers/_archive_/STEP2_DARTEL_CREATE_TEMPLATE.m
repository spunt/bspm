studydir = '/Users/bobspunt/Documents/fmri/conte';
subpat = 'CC*';
epipat = 'EP*';
t1pat = 'GR*T1_MEAN';
anatpat = [studydir filesep subpat filesep 'raw' filesep t1pat];
bob_spm_dartel_create_template(anatpat);