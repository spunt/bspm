function [out pval]=bramila_CMcorr(xSim,modelMat,iter,type)
	% [out pval]=bramila_mantel(matrix1,matrix2,iter,type)
	%
	% Mantel test for (dis)similarity matrices. The matrices must be squared, of same size, symmetrical and with ones in the main diagonal for the similarity case or zeros in the main diagonal 
	% for the distance (dissimilarity) matrix case. If distance matrix, then all values should be positive.
	% Mantel test is performed by correlating the top triangle between the two matrices. P values is obtained with permutations
	% Input parameter "type" can only be 'pearson' or 'spearman'
	%
	% e.g.:
	% 	a = corr(randn(100,10)); % a confusion matrix
	% 	b = corr(randn(100,10)); % another confusion matrix
	% 	[out pval] = bramila_CMcorr(a,b,5000,'spearman')
	%
	% (c) Enrico Glerean 2017 - Brain and Mind Laboratory Aalto University http://becs.aalto.fi/bml/

    if(strcmp(type,'pearson')==0 && strcmp(type,'spearman')==0)
		error('types allowed are only pearson and spearman');
	end
	if(iter<=0)
		warning('setting permutations to 1e5');
		iter=1e5;
	end

	% test that they are square
	kk1=size(xSim);
	kk2=size(modelMat);
	if(kk1(1) ~= kk1(2))	error('matrix 1 is not square'); end
	if(kk2(1) ~= kk2(2))	error('matrix 2 is not square'); end
	if(kk1(1) ~= kk2(1))	error('matrices are not of same size'); end

	% we do not test that they are symmetrical matrices
	% we do not test that they both are similarity or dissimilarity matrices


	ids=find(ones(kk1(1))); % all the IDs, in mantel test it would be the top triangle off diagonal
    out=corr(xSim(ids),modelMat(ids),'type',type);
	
    surro=zeros(iter,1);
    parfor i=1:iter
		pe=randperm(size(xSim,1));
		temp=xSim(pe,pe); % swap rows and columns
		% question: does it need to normalize so that sum along rows is = 1 (i.e. null confusion matrices).
		surro(i)=corr(temp(ids),modelMat(ids),'type',type);
   	end
	[fi xi]=ksdensity(surro,'function','cdf','npoints',200);
	pval_left=interp1([-1 xi 1],[0 fi 1],out);    % trick to avoid NaNs
	pval=1-pval_left;
end
    

        
