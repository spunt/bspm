a=lui_spm_select;
%a=lui_spm_load;

numfiles = size(a,1);

for j = 1:numfiles,
   disp(['working on file#' num2str(j)]);
   [V] = lui_spm_vol(a(j,:));
   data = lui_spm_read_vols(V);
   Lmsk = zeros(size(data));
   Rmsk = zeros(size(data));
   
   xdim = size(data,1);
   
   L=data(1:floor(xdim/2),:,:);
   R=data(end:-1:ceil(xdim/2)+1,:,:);
   
%   data(1:floor(xdim/2),:,:) = (R-L)./(R+L+eps);
%   data(1:floor(xdim/2),:,:) = atan(L./(R+eps))*180/pi;
   data(1:floor(xdim/2),:,:) = (L-R);
   
   data(end:-1:ceil(xdim/2)+1,:,:)=(R-L);
   
   if (xdim/2)~=round(xdim/2),
      data(floor(xdim/2)+1,:,:)=data(floor(xdim/2),:,:);
   end;
   
   [dr,name,ext,versn]=fileparts(a(j,:));
   fn = [name ext versn];

   V.fname = fullfile(dr,['lat_' fn]);
  %% V.dim(4)=16;
   lui_spm_write_vol(V,data);
   
   if (j==1),
      Lavg = L/numfiles;
      Ravg = R/numfiles;
   else,
      Lavg = Lavg+L/numfiles;
      Ravg = Ravg+R/numfiles;
   end;
   
   Lmsk(1:floor(xdim/2),:,:) = L;
   Rmsk(1:floor(xdim/2),:,:) = R;
   Lmsk(end:-1:ceil(xdim/2)+1,:,:)=L;
   Rmsk(end:-1:ceil(xdim/2)+1,:,:)=R;

   V.fname = fullfile(dr,['L_' fn]);
  %% V.dim(4)=16;
   lui_spm_write_vol(V,Lmsk);
   
   V.fname = fullfile(dr,['R_' fn]);
  %% V.dim(4)=16;
   lui_spm_write_vol(V,Rmsk);
   
end;

Lmsk(1:floor(xdim/2),:,:) = Lavg;
Rmsk(1:floor(xdim/2),:,:) = Ravg;
Lmsk(end:-1:ceil(xdim/2)+1,:,:)=Lavg;
Rmsk(end:-1:ceil(xdim/2)+1,:,:)=Ravg;

%masks
V.fname = fullfile(dr,['mask_Rp_Lp_' fn]);
%V.dim(4)=16;
lui_spm_write_vol(V,1*((Rmsk>0)&(Lmsk>0)));

V.fname = fullfile(dr,['mask_Rp_Ln_' fn]);
%V.dim(4)=16;
lui_spm_write_vol(V,1*((Rmsk>0)&(Lmsk<0)));

V.fname = fullfile(dr,['mask_Rn_Lp_' fn]);
%V.dim(4)=16;
lui_spm_write_vol(V,1*((Rmsk<0)&(Lmsk>0)));

V.fname = fullfile(dr,['mask_Rn_Ln_' fn]);
%V.dim(4)=16;
lui_spm_write_vol(V,1*((Rmsk<0)&(Lmsk<0)));

V.fname = fullfile(dr,['Lmean_' fn]);
%V.dim(4)=16;
lui_spm_write_vol(V,Lmsk);

V.fname = fullfile(dr,['Rmean_' fn]);
%V.dim(4)=16;
lui_spm_write_vol(V,Rmsk);

%Why not just L-R??????
%null hypothesis is 0.....
%units are in delta %signal strength

%ranges between 0 and 1
%(2/pi)*atan(abs(L)./(abs(R)+eps))

%ranges between -1 and 1 (better because null hypothesis is 0?s
%(abs(L)-abs(R))./(abs(L)+abs(R)+eps)

% out = zeros(size(L));
% for j = 1:prod(size(L)),
%     l = L(j);r = R(j);
%     if ((l~=0)&(r~=0)),
%         out(j) = (abs(l)-abs(r))./(abs(l)+abs(r)+eps);
%         if ((l>0)&(r>0)),
%             out(j) = out(j)/2;
%         elseif ((l>0)&(r<0)),
%             out(j) = -0.5-0.5*abs(out(j));
%         elseif ((l<0)&(r<0)),
%             out(j) = -out(j)/2;
%         elseif ((l<0)&(r>0)),
%             out(j) = 0.5+0.5*abs(out(j));
%         end;
%     end;
% end;
