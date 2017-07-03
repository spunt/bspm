function setDataPaths(Params)

Pub = Params.PublicParams;                                                                                                                              
Priv = Params.PrivateParams;                                                                                                                            

% check/create paths:                                                                                                                                   
if exist(Pub.dataDestination) ~= 7                                                                                                                      
  mkdir(Pub.dataDestination)                                                                                                                          
end                                                                                                                                                     
if exist(Priv.subjectDestination) ~= 7                                                                                                                  
  mkdir(Priv.subjectDestination)                                                                                                                      
end                                                                                                                                                     
if exist(Priv.subjectFiltDestination) ~= 7                                                                                                              
  mkdir(Priv.subjectFiltDestination)                                                                                                                  
end                                                                                                                                                     
if exist(Priv.resultsDestination) ~= 7                                                                                                                  
  mkdir(Priv.resultsDestination)                                                                                                                      
end                                                                                                                                                     
if exist(Priv.statsDestination) ~= 7                                                                                                                    
  mkdir(Priv.statsDestination)                                                                                                                        
end                                                                                                                                                     
if exist(Priv.PFDestination) ~= 7                                                                                                                       
  mkdir(Priv.PFDestination)                                                                                                                           
end                                                                                                                                                     
if exist(Priv.withinDestination) ~= 7                                                                                                                   
  mkdir(Priv.withinDestination)                                                                                                                       
end                                                                                                                                                     
if exist(Priv.phaseDifDestination) ~= 7                                                                                                                 
  mkdir(Priv.phaseDifDestination)                                                                                                                     
end                                             
%if exist(Priv.logDestination) ~= 7                                                                                                                 
%  mkdir(Priv.logDestination)                                                                                                                     
%end  
Tag = Pub.dataDescription;
save([Pub.dataDestination Pub.dataDescription],'Params');
save([Pub.dataDestination 'Tag'],'Tag');
