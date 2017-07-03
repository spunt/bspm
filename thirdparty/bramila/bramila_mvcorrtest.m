function pval = bramila_mvcorrtest(reg_cfg,iterations,corrtype)
%BRAMILA_TEST_TISSUE_ISC Summary of this function goes here
% %   Detailed explanation goes here

L=length(reg_cfg);
[total_X,total_subjID,total_typeID] = get_total_design(reg_cfg);
total_X = zscore(total_X); % removes the need for a constant term!

T=size(total_X,1);
[count,val]=hist(total_typeID,unique(total_typeID));
[~,i]=max(count);
M=rank(total_X(:,total_typeID==val(i)));

if nargin<3
   if T/M<2
       corrtype = 2;
   else
       corrtype = 1;
   end
end

if corrtype == 2
    fprintf('\nTesting multivariate ISC (pair-wise)\n');
elseif corrtype == 1
    fprintf('\nTesting multivariate ISC (full design) \n');
else
    corrtype
    error('Correlation type must be 1 (full) or 2 (pair-wise)');
end

shifts = randi(T,iterations,L,'uint16');
shifts(end+1,1:L) = ones(1,L,'uint16');

printint=round(linspace(iterations/10,iterations*9/10,9));


pval=nan(1,3);

for type=1:3            
    
    XXX = total_X(:,total_typeID==type);
    subj_id = total_subjID(total_typeID==type);
    dist=nan(1,iterations+1);
    kk=1;    
    if size(XXX,2)>0
        
        if type==1
            fprintf('..global signal: ');
        elseif type==2
            fprintf('..WM signal: ');
        elseif type==3
            fprintf('..CSF signal: ');
        end
        
        for iter=1:(iterations+1)
            
            if kk<10 && printint(kk)==iter
               fprintf('%2.0f%% ',100*iter/iterations);
               kk=kk+1;
            end
            
            XX = subjectwise_shift(XXX,L,subj_id,shifts(iter,:));
            
            %var_perc = zeros(1,L-1);
            k=0;
            for i=1:L
                Y=XX(:,subj_id==i);
                
                if corrtype==2
                    
                     for j=1:L
                        if i~=j
                        X=XX(:,subj_id==j);
                        if size(X,2)>0 && size(Y,2)>0                           
                            var_exp = 100*compute_variance(X,Y);
                            k=k+1;
                            var_perc(k)=var_exp;
                        end
                        end
                    end
                    
                else
                    
                    X=XX(:,subj_id~=i);
                    if size(X,2)>0 && size(Y,2)>0
                        var_exp = 100*compute_variance(X,Y);
                        k=k+1;
                        var_perc(k)=var_exp;
                    end
                    
                end
                
            end
            
            dist(iter)=mean(var_perc);
                        
        end
        
        if nnz(isnan(dist))>0
            warning('Distribution has NaN values!');
        end
        pval(type)=nnz(dist(iterations+1)<dist(1:iterations))/iterations;
        fprintf('100%% (pval=%3.2f)\n',pval(type));        
        
    end
    
end

pval(pval==0)=1/(iterations+1);

end

function XX = subjectwise_shift(XX,L,subj_id,shifts)

for i=1:L
    ind = (subj_id==i);
    XX(:,ind) = XX([shifts(i):end,1:(shifts(i)-1)],ind);
end

end

function var = compute_variance(X,Y)

beta = X \ Y;
var = sum(sum(abs(X*beta).^2,1))./ sum(sum(abs(Y).^2,1));

end

function [total_X,subj_id,type_id] = get_total_design(reg_cfg)

total_X = [];
subj_id =[];
type_id =[];

for subj = 1:length(reg_cfg)
    total_X = cat(2,total_X,reg_cfg{subj}.reg);
    subj_id = cat(2,subj_id,subj*ones(1,length(reg_cfg{subj}.id)));
    type_id = cat(2,type_id,reg_cfg{subj}.id);
end

end
