function nii_jaccard_batch
%Applies nii_jaccard to many images

out= 'results.tab';
 subj = strvcat('IBSR_01.nii','IBSR_02.nii','IBSR_03.nii','IBSR_04.nii','IBSR_05.nii','IBSR_06.nii','IBSR_07.nii','IBSR_08.nii',...
     'IBSR_09.nii','IBSR_10.nii','IBSR_11.nii','IBSR_12.nii','IBSR_13.nii','IBSR_14.nii','IBSR_15.nii','IBSR_16.nii','IBSR_17.nii','IBSR_18.nii');
 dir = '/Users/rorden/segfinal/ibsr/manual/';
 b = zeros(2,length(subj));
 binary = true;

dirt = '/Users/rorden/segfinal/ibsr/usn/';
for i=1:size(subj,1)
   b(1,i) = nii_jaccard(fullfile(dir,['s' deblank(subj(i,:))]),fullfile(dirt,['c1' deblank(subj(i,:))]),binary,2 );
   b(2,i) = nii_jaccard(fullfile(dir,['s' deblank(subj(i,:))]),fullfile(dirt,['c2' deblank(subj(i,:))]),binary,3 );
end;
dlmwrite(out,'usn','-append', 'delimiter', '');
dlmwrite(out, b, '-append', 'delimiter', '\t');

dirt = '/Users/rorden/segfinal/ibsr/newseg/';
for i=1:size(subj,1)
   b(1,i) = nii_jaccard(fullfile(dir,['s' deblank(subj(i,:))]),fullfile(dirt,['c1' deblank(subj(i,:))]),binary,2 );
   b(2,i) = nii_jaccard(fullfile(dir,['s' deblank(subj(i,:))]),fullfile(dirt,['c2' deblank(subj(i,:))]),binary,3 );
end;
dlmwrite(out,'newseg_nocleanup','-append', 'delimiter', '');
dlmwrite(out, b, '-append', 'delimiter', '\t');


dirt = '/Users/rorden/segfinal/ibsr/newseg/';
for i=1:size(subj,1)
   b(1,i) = nii_jaccard(fullfile(dir,['s' deblank(subj(i,:))]),fullfile(dirt,['ec1' deblank(subj(i,:))]),binary,2 );
   b(2,i) = nii_jaccard(fullfile(dir,['s' deblank(subj(i,:))]),fullfile(dirt,['ec2' deblank(subj(i,:))]),binary,3 );
end;
dlmwrite(out,'newseg','-append', 'delimiter', '');
dlmwrite(out, b, '-append', 'delimiter', '\t');

dirt = '/Users/rorden/segfinal/ibsr/etpm/';
for i=1:size(subj,1)
   b(1,i) = nii_jaccard(fullfile(dir,['s' deblank(subj(i,:))]),fullfile(dirt,['ec1' deblank(subj(i,:))]),binary,2 );
   b(2,i) = nii_jaccard(fullfile(dir,['s' deblank(subj(i,:))]),fullfile(dirt,['ec2' deblank(subj(i,:))]),binary,3 );
end;
dlmwrite(out,'etpm','-append', 'delimiter', '');
dlmwrite(out, b, '-append', 'delimiter', '\t');

subj = strvcat('1_24.nii','2_4.nii','4_8.nii','5_8.nii','6_10.nii','7_8.nii','8_4.nii','11_3.nii','12_3.nii','13_3.nii','15_3.nii','17_3.nii','16_3.nii','100_23.nii','110_3.nii','111_2.nii','112_2.nii','191_3.nii','202_3.nii','205_3.nii');
dir = '/Users/rorden/segfinal/ibsrold/manual/';

dirt = '/Users/rorden/segfinal/ibsrold/usn/';
for i=1:size(subj,1)
   b(1,i) = nii_jaccard(fullfile(dir,['c' deblank(subj(i,:))]),fullfile(dirt,['c1' deblank(subj(i,:))]),binary,193 );
   b(2,i) = nii_jaccard(fullfile(dir,['c' deblank(subj(i,:))]),fullfile(dirt,['c2' deblank(subj(i,:))]),binary,255 );
end;
dlmwrite(out,'usn','-append', 'delimiter', '');
dlmwrite(out, b, '-append', 'delimiter', '\t');

dirt = '/Users/rorden/segfinal/ibsrold/newseg/';
for i=1:size(subj,1)
   b(1,i) = nii_jaccard(fullfile(dir,['c' deblank(subj(i,:))]),fullfile(dirt,['c1' deblank(subj(i,:))]),binary,193 );
   b(2,i) = nii_jaccard(fullfile(dir,['c' deblank(subj(i,:))]),fullfile(dirt,['c2' deblank(subj(i,:))]),binary,255 );
end;
dlmwrite(out,'newseg_nocleanup','-append', 'delimiter', '');
dlmwrite(out, b, '-append', 'delimiter', '\t');

dirt = '/Users/rorden/segfinal/ibsrold/newseg/';
for i=1:size(subj,1)
   b(1,i) = nii_jaccard(fullfile(dir,['c' deblank(subj(i,:))]),fullfile(dirt,['ec1' deblank(subj(i,:))]),binary,193 );
   b(2,i) = nii_jaccard(fullfile(dir,['c' deblank(subj(i,:))]),fullfile(dirt,['ec2' deblank(subj(i,:))]),binary,255 );
end;
dlmwrite(out,'newseg','-append', 'delimiter', '');
dlmwrite(out, b, '-append', 'delimiter', '\t');

dirt = '/Users/rorden/segfinal/ibsrold/etpm/';
for i=1:size(subj,1)
   b(1,i) = nii_jaccard(fullfile(dir,['c' deblank(subj(i,:))]),fullfile(dirt,['ec1' deblank(subj(i,:))]),binary,193 );
   b(2,i) = nii_jaccard(fullfile(dir,['c' deblank(subj(i,:))]),fullfile(dirt,['ec2' deblank(subj(i,:))]),binary,255 );
end;
dlmwrite(out,'etpm','-append', 'delimiter', '');
dlmwrite(out, b, '-append', 'delimiter', '\t');







