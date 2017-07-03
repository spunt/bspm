function [out pval]=bramila_mantel(xSim,modelMat,iter,type)
	% [out pval]=bramila_mantel(matrix1,matrix2,iter,type)
	%
	% Mantel test for (dis)similarity matrices. The matrices must be squared, of same size, symmetrical and with ones in the main diagonal for the similarity case or zeros in the main diagonal 
	% for the distance (dissimilarity) matrix case. If distance matrix, then all values should be positive.
	% Mantel test is performed by correlating the top triangle between the two matrices. P values is obtained with permutations
	% Input parameter "type" can only be 'pearson' or 'spearman'
	%
	% e.g.:
	% 	a = corr(randn(100,10));
	% 	b = corr(randn(100,10));
	% 	[out pval] = bramila_mantel(a,b,5000,'spearman')
	%
	% (c) Enrico Glerean 2013 - Brain and Mind Laboratory Aalto University http://becs.aalto.fi/bml/

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

	% test that they are symmetrical matrices
	temp=xSim-xSim';
	if(sum(temp(:))~=0) error('matrix 1 is not symmetrical'); end
	temp=modelMat -modelMat';
	if(sum(temp(:))~=0) error('matrix 2 is not symmetrical'); end
	% test that they both are similarity or dissimilarity matrices
	diag1=sum(diag(xSim));
	diag2=sum(diag(modelMat));
	if(diag1~=kk1(1) && diag1~= 0) error('matrix 1 is not a similarity or a distance matrix'); end 
	if(diag2~=kk2(1) && diag2~= 0) error('matrix 2 is not a similarity or a distance matrix'); end 
	if(diag1~=diag2)	error('matrices are not of the same type (one is similarity and the other distance)'); end


	ids=find(triu(ones(kk1(1)),1));
    out=corr(xSim(ids),modelMat(ids),'type',type);
	
    surro=zeros(iter,1);
    parfor i=1:iter
		pe=randperm(size(xSim,1));
		temp=xSim(pe,pe);
		surro(i)=corr(temp(ids),modelMat(ids),'type',type);
   	end
	[fi xi]=ksdensity(surro,'function','cdf','npoints',200);
	pval_left=interp1([-1 xi 1],[0 fi 1],out);    % trick to avoid NaNs
	pval=1-pval_left;
end
    

        
