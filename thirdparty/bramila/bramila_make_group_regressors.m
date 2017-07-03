function [reg_cfg,var_expl,cv_data] = bramila_make_group_regressors(reg_cfg,method,selected_tissues)
% removes inter-subject correlations from individual tissue regressors
% method = 'full' does normal full-model regression
% method = 'regularized' does partial least squares regression (PLS) with
% 10-fold cross validation to avoid possible overfitting

if nargin<3
   selected_tissues=[1,1,1];
end

var_expl = nan(length(reg_cfg),3);
cv_data=[];

if any(selected_tissues==1)
   fprintf('\nRemoving shared signal from nuisance regressors\n');
else
   fprintf('\nNo tissues selected, skippping\n');
   return;
end

[total_X,subj_id,type_id,labels] = get_total_design(reg_cfg);
total_X = zscore(total_X);

total_X_global = total_X(:,type_id==1);
subj_id_global = subj_id(type_id==1);
labels_global = labels(type_id==1);

total_X_wm = total_X(:,type_id==2);
subj_id_wm = subj_id(type_id==2);
labels_wm = labels(type_id==2);

total_X_csf = total_X(:,type_id==3);
subj_id_csf = subj_id(type_id==3);
labels_csf = labels(type_id==3);

s = RandStream('mt19937ar','Seed',100*sum(clock()));
RandStream.setGlobalStream(s);
saved_stream = s.State;

for subj = 1:length(reg_cfg)
       
    % all subjects have the same cv partitioning
    s.State = saved_stream;
        
    X_global=total_X_global(:,subj_id_global~=subj);
    Y_global=total_X_global(:,subj_id_global==subj);
    global_model_size = size(X_global,2);
    if selected_tissues(1)==1 && size(X_global,2)>0 && size(Y_global,2)>0     
        if strcmp(method,'pls')
            [order,vari,cv_mse1]=get_PLS_order(X_global,Y_global);
            perc_global = vari(order);
            [~,~,~,~,~,~,~,stats] = plsregress(X_global,Y_global,order);
            global_order=order;
        elseif strcmp(method,'pca')
            [X_global,order,~,cv_mse1]=get_PCA_order(X_global,Y_global);
            [~,~,~,~,~,vari,~,stats] = plsregress(X_global,Y_global,size(X_global,2));
            perc_global = 100*sum(vari(2,1:order));        
            global_order=order;
        elseif strcmp(method,'full')            
            [~,~,~,~,~,PCTVAR,~,stats] = plsregress(X_global,Y_global,size(X_global,2));
            var_exp = 100*cumsum(PCTVAR(2,:));
            perc_global=var_exp(end);
            global_order=size(X_global,2);
        else
             error('Unknown regression method: ''%s''',method);
        end
        Y_global=stats.Yresiduals;
    else
        global_order=0;
        perc_global=0;
        cv_mse1=[];
    end        
    
    X_wm=total_X_wm(:,subj_id_wm~=subj);
    Y_wm=total_X_wm(:,subj_id_wm==subj);
    wm_model_size = size(X_wm,2);
    if selected_tissues(2)==1 && size(X_wm,2)>0 && size(Y_wm,2)>0      
        if strcmp(method,'pls')            
            [order,vari,cv_mse2]=get_PLS_order(X_wm,Y_wm);
            perc_wm = vari(order);
            [~,~,~,~,~,~,~,stats] = plsregress(X_wm,Y_wm,order);
            wm_order=order;
        elseif strcmp(method,'pca')
            [X_wm,order,~,cv_mse2]=get_PCA_order(X_wm,Y_wm);            
            [~,~,~,~,~,vari,~,stats] = plsregress(X_wm,Y_wm,size(X_wm,2));
            perc_wm = 100*sum(vari(2,1:order));
            wm_order=order;
        elseif strcmp(method,'full')            
            [~,~,~,~,~,PCTVAR,~,stats] = plsregress(X_wm,Y_wm,size(X_wm,2));
            var_exp = 100*cumsum(PCTVAR(2,:));
            perc_wm=var_exp(end);
            wm_order=size(X_wm,2);
        else
             error('Unknown regression method: ''%s''',method);
        end
        Y_wm=stats.Yresiduals;
    else
        wm_order=0;
        perc_wm=0;
        cv_mse2=[];
    end
    
    X_csf=total_X_csf(:,subj_id_csf~=subj);
    Y_csf=total_X_csf(:,subj_id_csf==subj);
    csf_model_size = size(X_csf,2);
    if selected_tissues(3)==1 && size(X_csf,2)>0 && size(Y_csf,2)>0    
        if strcmp(method,'pls')
            [order,vari,cv_mse3]=get_PLS_order(X_csf,Y_csf);
            perc_csf = vari(order);
            [~,~,~,~,~,~,~,stats] = plsregress(X_csf,Y_csf,order);
            csf_order=order;
        elseif strcmp(method,'pca')
            [X_csf,order,~,cv_mse3]=get_PCA_order(X_csf,Y_csf);            
            [~,~,~,~,~,vari,~,stats] = plsregress(X_csf,Y_csf,size(X_csf,2));
            perc_csf = 100*sum(vari(2,1:order));  
            csf_order=order;
        elseif strcmp(method,'full')
            [~,~,~,~,~,PCTVAR,~,stats] = plsregress(X_csf,Y_csf,size(X_csf,2));
            var_exp = 100*cumsum(PCTVAR(2,:));
            perc_csf=var_exp(end);
            csf_order=size(X_csf,2);
        else
             error('Unknown regression method: ''%s''',method);
        end
        Y_csf=stats.Yresiduals;
    else
        csf_order=0;
        perc_csf=0;
        cv_mse3=[];
    end
    
    if strcmp(method,'pls') || strcmp(method,'pca')
        if perc_global>0
            fprintf('..subject %i: Global=%4.1f%% (%i/%i), WM=%4.1f%% (%i/%i), CSF=%4.1f%% (%i/%i)\n',subj,perc_global,global_order,global_model_size,perc_wm,wm_order,wm_model_size,perc_csf,csf_order,csf_model_size);
            cv_data{subj,1}=cv_mse1;
            cv_data{subj,2}=cv_mse2;
            cv_data{subj,3}=cv_mse3;
        else
            fprintf('..subject %i: WM=%4.1f%% (%i/%i), CSF=%4.1f%% (%i/%i)\n',subj,perc_wm,wm_order,wm_model_size,perc_csf,csf_order,csf_model_size);
            cv_data{subj,1}=[];
            cv_data{subj,2}=cv_mse2;
            cv_data{subj,3}=cv_mse3;            
        end
    elseif strcmp(method,'full')
        if perc_global>0
            fprintf('..subject %i: Global=%4.1f%% (%i), WM=%4.1f%% (%i), CSF=%4.1f%% (%i)\n',subj,perc_global,global_model_size,perc_wm,wm_model_size,perc_csf,csf_model_size);
        else
            fprintf('..subject %i: WM=%4.1f%% (%i), CSF=%4.1f%% (%i)\n',subj,perc_wm,wm_model_size,perc_csf,csf_model_size);
        end
    else
        error('Unknown regression method: ''%s''',method);
    end
    
    var_expl(subj,1:3)=[perc_global,perc_wm,perc_csf];    
    
    % replace regressors with clean ones
    reg_cfg{subj}.reg=[Y_global,Y_wm,Y_csf];
    % regressor id's (1=global,2=white,3=csf)
    reg_cfg{subj}.id=[ones(1,size(Y_global,2)),2*ones(1,size(Y_wm,2)),3*ones(1,size(Y_csf,2))];
    % regressor labels
    reg_cfg{subj}.labels = [labels_global(subj_id_global==subj),labels_wm(subj_id_wm==subj),labels_csf(subj_id_csf==subj)];

end
end

function [total_X,subj_id,type_id,labels] = get_total_design(reg_cfg)

total_X = [];
subj_id =[];
type_id =[];
labels=[];

for subj = 1:length(reg_cfg)
    total_X = cat(2,total_X,reg_cfg{subj}.reg);
    subj_id = cat(2,subj_id,subj*ones(1,length(reg_cfg{subj}.id)));
    type_id = cat(2,type_id,reg_cfg{subj}.id);
    labels = cat(2,labels,reg_cfg{subj}.labels);
end

end

function [order,var_exp,cv_mse] = get_PLS_order(X,Y)

[~,~,~,~,~,PLSPctVar,cv_data] = bramila_plsregress(X,Y,size(X,2),'cv',10,'mcreps',10);

% get response variance
var_exp = 100*cumsum(PLSPctVar(2,:));
% convert between [0,100]
var_exp_perc = 100*var_exp/var_exp(end);
if abs(100-var_exp_perc(end))>0.00001
    error('Cumulative re-scaled variance percentage not 100%% (BUG)')
end

% get MSE related to response
cv_mse = cv_data.MSE(2,:);
cv_se = cv_data.SE(2,:);

[limit_cv,limit_var]=compute_limits(var_exp_perc,cv_mse,cv_se);

order = min(limit_cv,limit_var);
order=max(1,order);

end

function [X_out,order,var_exp,cv_mse] = get_PCA_order(X,Y)

% global isdebug;
[~,PCAScores,PLSPctVar,cv_data] = bramila_pcaregress(X,Y,size(X,2),10,10);

% get response variance (already cumulative)
var_exp = 100*(PLSPctVar(2,:));
% convert between [0,100]
var_exp_perc = 100*var_exp/var_exp(end);
if abs(100-var_exp_perc(end))>0.00001
    error('Cumulative re-scaled variance percentage not 100%% (BUG)')
end
% get MSE related to response
cv_mse = cv_data.MSE;
cv_se = cv_data.SE;

%a=mean(cv_mse);
%plot(0:(length(cv_mse)-1),cv_mse-a,'x-','Color',0.2+0.8*rand(1,3));hold on;
%plot(0:(length(cv_mse)-1),cv_mse-a+cv_se,'g--');
[limit_cv,limit_var]=compute_limits(var_exp_perc,cv_mse,cv_se);

% get final model order
order = min(limit_cv,limit_var);
order = max(1,order);

X_out = PCAScores(:,1:order);

end

function [limit_cv,limit_var]=compute_limits(var_exp_perc,cv_mse,cv_se)
% Estimate optimal model order based on MSE, MSE-SE and total response variance
% explained

if length(var_exp_perc)+1 ~= length(cv_mse)
    error('Variance vector length incorrect (BUG)')
end
if length(cv_mse) ~= length(cv_se)
    error('length of MSE and SE vectors mismatch (BUG)')
end

% compute extrema distance
mi=min(cv_mse);
ma=max(cv_mse);
dev = ma-mi;

% locate all local minima
[val,loc] = findpeaks(-[inf,cv_mse,inf]);
val = -val;
loc(loc==1 | loc==length(cv_mse)+2)=[];
loc=loc-1;

% remove minima that lie beyond high MSE barrier
for i=1:(length(loc)-1)
    maxdev = -inf;
    for j = (loc(i)+1):loc(i+1)
        maxdev = max(maxdev,cv_mse(j)-cv_mse(loc(i)));
    end
    if maxdev>dev/4
        loc((i+1):end)=[];
        val((i+1):end)=[];        
        break;
    end
end

% get CV minima within one SE
mse = cv_mse(1:loc(end));
se = cv_se(1:loc(end));
minMSE = min(mse);
minIx = find(mse==minMSE,1,'first');
minplus = mse(minIx) + se(minIx);
limit_cv = find((mse(1:minIx) <= minplus),1,'first') - 1;

% get upper limit
limit_var = find(var_exp_perc>=95,1,'first');

end
