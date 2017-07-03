function [cout pval]=bramila_spatialCorr(cfg)
% compare two spatial maps using bootstrap resampling based on 3D fourier transform
%
% use as 
%   [cout, pval] = bramila_spatialCorr(cfg)
% cfg has to have
% cfg.infileA= full path to single volume nifti file to test (condition A)
% cfg.infileB= full path to single volume nifti file to test (condition B)
% cfg.mask= full path to nifti mask (volume with only 0 and 1)
% cfg.niter=5000 %number of permutations
%
% function returns
%   cout = within mask spatial correlations using Pearson and Spearman
%   pval = within mask estimated pvalue using Pearson and Spearman 
% if correlations are negative please check 1-pval as a result


NITER=cfg.niter;
niiA=load_nii(cfg.infileA);
niiB=load_nii(cfg.infileB);
niimask=load_nii(cfg.mask);
mask=sign(niimask.img);
if (any(find(mask<0)))
    error('mask should only contain ones or zeros');
end

inmask=find(mask>0);
dataA=double(niiA.img);
dataB=double(niiB.img);

% used for bootstrap resampling
dil_ft=fftn(abs(dataA)); % three dimensional fourier transform 
a_dil_ft=angle(dil_ft); % phase
m_dil_ft=abs(dil_ft);   % magnitude


cout(1)=corr(dataA(inmask),dataB(inmask));
cout(2)=corr(dataA(inmask),dataB(inmask),'type','Spearman');

% do permutations for both pearson and spearman
surro=zeros(NITER,2);
parfor iter=1:NITER
	tempA=a_dil_ft(:);
	tempA=tempA(randperm(length(tempA)));
	tempA=reshape(tempA,size(a_dil_ft));
	tempI=real(ifftn(m_dil_ft.*exp(1i*tempA)));
		cpear=corr(tempI(inmask),dataB(inmask));
		cspear=corr(tempI(inmask),dataB(inmask),'type','Spearman');            
	surro(iter,:)=[cpear cspear];
end
for pear_spear=1:2
	[fi xi]=ksdensity(surro(:,pear_spear),'function','cdf','npoints',200);
	pval(pear_spear)=1-interp1([-1 xi 1],[0 fi 1],cout(pear_spear));    % trick to avoid NaNs
end

