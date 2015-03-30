epi             = grabfiles('RA*/raw', 'EP*/ua*nii');
[~, subname]    = files('RA*'); 
for i = 1:length(epi)
   printmsg(subname{i}); 
   bfsl_bet_epi(epi{i}, .3, 'b');  
end
