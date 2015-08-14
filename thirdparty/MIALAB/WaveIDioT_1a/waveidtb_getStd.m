% function [sighhh_hhl, siglhl_lhh, sighll_hlh , siglll_llh ] = getStd( hhh , hhl , llh )
%  
%          
%     sigma = mean([(median(abs(hhl{1}(:)))/0.6745) (median(abs(hhh{1}(:)))/0.6745)]);
%     
% %     sighhh_hhl = [1.00 0.22 0.09 0.04]*sigma;
% %     siglhl_lhh = [0.79 0.25 0.11 0.05]*sigma;
%     
%     sighhh_hhl = [1.00 0.22 0.09 0.04]*sigma;
%     siglhl_lhh = [0.79 0.25 0.11 0.05]*sigma;
%     
%     sighll_hlh = siglhl_lhh;
%     siglll_llh = sighll_hlh;


function [sighhh_hhl, siglhl_lhh, sighll_hlh , siglll_llh ] = getStd( lhl, lhh , hll , hlh , hhl , hhh , lll, llh  )
 
    sigInput_hh = mean( [ median(abs(hhh{1}(:)))/0.6745, median(abs(hhl{1}(:)))/0.6745  ] );
    sigInput_hl = mean( [ median(abs(hlh{1}(:)))/0.6745, median(abs(hll{1}(:)))/0.6745  ] );
    sigInput_lh = mean( [ median(abs(lhh{1}(:)))/0.6745, median(abs(lhl{1}(:)))/0.6745  ] );
    sigInput_ll = mean( [ median(abs(lll{1}(:)))/0.6745, median(abs(llh{1}(:)))/0.6745  ] );
         
   
    sighhh_hhl = [1 1 1 1]*sigInput_hh;
    siglhl_lhh = [1 1 1 1]*sigInput_lh;
    sighll_hlh = [1 1 1 1]*sigInput_hl;
    siglll_llh = [1 1 1 1]*sigInput_ll;
