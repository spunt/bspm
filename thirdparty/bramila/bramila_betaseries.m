function betas=bramila_betaseries(cfg)
% bramila_betaseries
% see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4019671
% see https://www.ncbi.nlm.nih.gov/pubmed/21924359
%
% betas=bramila_betaseries(cfg)
%   cfg.y= time series of ROIs (time 1st dimension)
%   cfg.events = cell matrix, each cell contains a list of time points of
%   event, one cell per condition
%   cfg.TR= TR
%   cfg.hrf = hrf (see bramila_hrf)

T=size(cfg.y,1); % number of time points
NC=length(cfg.events); % number of conditions
R=size(cfg.y,2); % number of ROIs

allEvents=[];
allEventsClasses=[];
for cc=1:NC
    eventVec=cfg.events{cc};
    if(size(eventVec,2)~=1) error('Event is not a column vector'); end
    allEvents=[allEvents;eventVec];
    allEventsClasses=[allEventsClasses;cc*ones(size(eventVec))];
end

betavec=zeros(length(allEvents),R);
for ev=1:length(allEvents)
    % create a regressor for the event and a regressor with all other
    % events
    curr_ev=zeros(T,1);
    all_other_ev=zeros(T,1);
    for t=1:length(allEvents)
        if(t==ev)
            curr_ev(round(allEvents(t)/cfg.TR))=1;
        else
            all_other_ev(round(allEvents(t)/cfg.TR))=1;
        end
    end
    
    x1=conv(curr_ev,cfg.hrf);
    x1=x1(1:T);
    
    x2=conv(all_other_ev,cfg.hrf);
    x2=x2(1:T);
    
    % run regression
    for r=1:R 
        b=regress(cfg.y(:,r),[x1 x2 ones(T,1)]);
        betavec(ev,r)=b(1);
    end
end

% reshape output
for cc=1:NC
    temp=betavec(find(allEventsClasses==cc),:);
    betas{cc}=temp;
    clear temp;
end
    
