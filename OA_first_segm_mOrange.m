clear all
%-----this is code for indicating which cells were previously choosen------
% %---------------------------------------------------------
%   clear all
%   load('AD_whi5X3(only)_30min_met_may08_2010_pos_no_4_1')
%   %imshow(bwmorph(all_obj.cells(:,:,180)>0,'remove',inf))
%   c_num   ='c1';%phase
%   alpha_add_t=numbM-0;
%   [im_name] = get_image_name(prefix,pos_num,suffix,alpha_add_t, c_num,type,suffix2,numbM);
%   I=imread(im_name);
%   [x_s,y_s]=size(I);
%   Inew=uint8(zeros(x_s,y_s,3));
%   Inew(:,:,1)=I +255.* uint8(bwmorph(all_obj.cells(:,:,numbM)>0,'remove',inf));
%   Inew(:,:,2)=I;
%   Inew(:,:,3)=I;
%   figure(9)
%   imshow(Inew)
% % %---------------------------------------------------------
%1 2 3 4 20 21 23 24
% 1 3 10 11 17 19
for pos_here=[2 3 5 7 8 10:12] 
%clear all
prefix  ='OA_071015_OA042_refurb_3nM6min';
suffix2='3';
pos_num =pos_here;%1;
suffix  ='_';%_p2z2
numbM   =90; %the (last)timepoint
type    ='.tif';
c_num   ='c1';%phase
alpha_add_t=numbM-0;
[im_name] = get_image_name(prefix,pos_num,suffix,alpha_add_t, c_num,type,suffix2,numbM);
I=imread(im_name);

c_num   ='c3';%whi5-orange
[im_name] = get_image_name(prefix,pos_num,suffix,alpha_add_t, c_num,type,suffix2,numbM);
IG=imread(im_name);
c_num   ='c2';%gfp
[im_name] = get_image_name(prefix,pos_num,suffix,alpha_add_t, c_num,type,suffix2,numbM);
IY=imread(im_name);

%==============================
%I =0.33.*I(:,:,1)+0.33.*I(:,:,2)+0.33.*I(:,:,3);
%IG=IG(:,:,1);

disp(im_name)
bw=I<80;
bw2=bwmorph(bw,'thicken',1);
D = bwdist(~bw2);%original
D = -D;
D(~bw2) = -Inf; %note it is not bw2 but bw!
L = watershed(D);
%[Ifin2]=manually_get_labeled_image((L>0).*bw,I);
[Ifin2]= OA_manually_get_labeled_image((L>0).*bw,I,IG,IY); 
Lcells=bwlabel(Ifin2);
LcellsO=Lcells;%o for original
%name3=['tmpsave' num2str(pos_num)];
name2=['tmpsave' prefix '_pos' num2str(pos_num)];
save(name2, 'LcellsO')
close all
%---------------------------------------------------------
end


