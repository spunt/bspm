%% Results for NT modules

% computes the Fearson's correlation for reference subnetworks as defined
% in Yeo2011 and Power 2011 and returns the highest correlating network label
% as reported in Table S1

cfg=[];
cfg.reference='power';
types={'yeo','power','cole','gordon'}; % we will only use the first 2

for n=1:12 % for each NT module
    fprintf(['Module ' num2str(n) '|']) %
    for atlas=1:2
        cfg.reference=types{atlas};
        cfg.infile=['/Users/enrico/Documents/phd/writings/ongoing/ASN/git/hfasdmodules/mesoscopic/mod_NT' num2str(n) '_nointerp.nii'];
        [scores labels]=bramila_compareModules(cfg);
        wincorr=find(max((scores(:,1)))==(scores(:,1)));
        fprintf([cfg.reference ' ' labels{scores(wincorr,2)} ' (' num2str(scores(wincorr,1),2) ') | ' ]);
    end
    disp(' ')
end