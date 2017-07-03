function showTime(startTime)
  
if startTime
  str = 'Started';
else
  str = 'Finished';
end
  
wt = clock;
if wt(5) < 10
  mins = ['0' num2str(wt(5))];
else
  mins = num2str(wt(5));
end
disp([str ' at: ' date ', ' num2str(wt(4)) ':' mins ':' num2str(wt(6))])  
