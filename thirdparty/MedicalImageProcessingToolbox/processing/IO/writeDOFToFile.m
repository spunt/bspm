function writeDOFToFile(DOF,DOFfn,n)

% DOF: matrix with the dof values
% DOFfn: filename
% n: number of DOF's
    
     fid=fopen(DOFfn,'w','native');
     fprintf(fid,['DOF: ' num2str(n) '\n'] );
     for i =1:n
        fprintf(fid,'%G \t %G \t %G \n', DOF(i,1),DOF(i,2),DOF(i,3));
     end
     
     fclose(fid);
end