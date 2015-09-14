function [all_obj]=get_fl_OA_1(all_obj,i,prefix,suffix2,pos_num,suffix,numbM,type,Lfinal,cell_exists,cell_no)

%this is a program for detecting (nuclear) fus1 & a weak mCherry (w?)
%nuclear localization okt12 2010
%better gausian nuclar fittign added mars 17 2011

%new version as of april 19th 2013
%cell_outline=(new_c_Image2>0

% figure(3);imagesc(Lfinal)

gauss_peak_cutoff=0.85; %where we cutoff the 2D-gaussian fit
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c_num   ='c2';%gfp <-use c4 for F1E strain, comment added nov21 2013
[im_name] = get_image_name(prefix,pos_num,suffix,i, c_num,type,suffix2,numbM);
IN=imread(im_name);
%==============================
%IN =0.5.*IN(:,:,1)+0.5.*IN(:,:,2);

[x_size,y_size]=size(IN);
all_obj.med_backgr_Far1(1,i) =median(median(IN));

all_obj.med_backgr_gfp(1,i)= double(median(median(IN)));

j=cell_no;
if cell_exists(j,1)==1
    put_mC=uint8(Lfinal==j).*IN;
    all_obj.max_nucl_int(j,i)=max(max(put_mC));
    all_obj.mean_Far1_int_per_area(j,i)=sum(sum(double(put_mC)))./sum(sum(double(put_mC>0)));
end

%--------------next kind of fit not using the nucleus....------------------
%------------------------------------------------------------------------
IG=IN;%tmp reassignment for simplicity
j=cell_no;
if cell_exists(j,1)==1
    put_whi5=uint8(Lfinal==j).*IG;
    %-----gaussian fit--------------------
    h = fspecial('average', 5);
    I_gf=  imfilter(IG,h,'replicate');
    Ignew=I_gf;
    nuc_tmp=Ignew.*uint8(Lfinal==j);
    
    Amp=max(max(nuc_tmp));
    s_f=regionprops(bwlabel(nuc_tmp==Amp),'Centroid');
    xtmpcenter=s_f(1).Centroid(1,1);
    ytmpcenter=s_f(1).Centroid(1,2);
    wind_size=10;
    xtmpcoordhere=max(1,xtmpcenter-wind_size):min(y_size,xtmpcenter+wind_size);
    ytmpcoordhere=max(1,ytmpcenter-wind_size):min(x_size,ytmpcenter+wind_size);
    [X, Y] = meshgrid(xtmpcoordhere,ytmpcoordhere);
    
    sX=size(X);
    sY=size(Y);
    ok_size=wind_size.*2+1;
    if sX(1,1)==ok_size && sX(1,2)==ok_size && sY(1,1)==ok_size && sY(1,2)==ok_size
        In_part=double(nuc_tmp(round(diag(Y)),round(diag(X))));
        edge_med=median([In_part(1,:) In_part(end,:) In_part(:,1)' In_part(:,end)']);
        x0=xtmpcenter;
        y0=ytmpcenter;
        mu=[x0,y0];
        best_fit_score=1e9;
        best_fit_no=[0 0];
        for jj=1:2:50
            for jj2=1:5:50
                limit_h=round(sqrt(jj.*jj2));
                for ii=-(limit_h-1):2:(limit_h-1)
                    Sigma = [jj -ii; -ii jj2];
                    F = mvnpdf([X(:) Y(:)],mu,Sigma);
                    F = reshape(F,length(X(:,1)),length(Y(:,1)));
                    Fmax=max(max(F));
                    Z=(( ((double(Amp)-edge_med)./Fmax)).*F+edge_med);
                    Z=Z.*(In_part>0);
                    tmp_score=sum(sum(abs((Z-In_part))));
                    if tmp_score<best_fit_score
                        best_fit_no=Sigma;%[ii jj];
                        best_fit_score=tmp_score;
                    end
                end
            end
        end
        %----------------- get the best solution------------------------
        %Sigma = [best_fit_no(1,2) -best_fit_no(1,1); -best_fit_no(1,1) best_fit_no(1,2)];
        %F = mvnpdf([X(:) Y(:)],mu,Sigma);
        F = mvnpdf([X(:) Y(:)],mu,best_fit_no);
        F = reshape(F,length(X(:,1)),length(Y(:,1)));
        Fmax=max(max(F));
        Z=(( ((double(Amp)-edge_med)./Fmax)).*F+edge_med);
        Z=Z.*(In_part>0);
        %---------------------------------------------------------------
        p_nuc=zeros(x_size,y_size);
        p_nuc(round(diag(Y)),round(diag(X)))=p_nuc(round(diag(Y)),round(diag(X)))+(Z>(Amp.*gauss_peak_cutoff));
        %---------------------------------------------------------------------
        p_cyt=bwmorph((Lfinal==j).*bwmorph(~p_nuc,'thicken',1),'erode',2);
        nuc_area=sum(sum(p_nuc));
        cyt_area=sum(sum(p_cyt));
        all_obj.tot_nucl_areaR(j,i)= nuc_area;
        all_obj.tot_cyt_areaR(j,i) = cyt_area;
        all_obj.cyt_Far1(j,i)      = sum(sum(double(p_cyt).*double(put_whi5)))./cyt_area;
        all_obj.nuc_Far1(j,i)      = sum(sum(double(p_nuc).*double(put_whi5)))./nuc_area;
        all_obj.nuclLFar1(:,:,i)   = all_obj.nuclLFar1(:,:,i) +p_nuc;%+j.* p_nuc;
        
        % ring around nucleus
        rad = 4;
        nuc=(all_obj.cells(:,:,i)==j).*(all_obj.nuclLFar1(:,:,i));
        OutI=((bwmorph(nuc,'thicken',rad).*~bwmorph(nuc,'thicken',rad-2)).*(all_obj.cells(:,:,i)==j));
        ring = double(OutI).*double(put_whi5);
        all_obj.med_nuc_GFP(j,i)= median(ring(ring>0));
        ring_sig = all_obj.med_nuc_GFP(j,i) - all_obj.med_backgr_gfp(1,i);
        cyt = double(p_cyt).*double(put_whi5);
        sig_cyt = cyt(cyt > 0) - all_obj.med_backgr_gfp(1,i);
        sig_cyt(sig_cyt > ring_sig) = ring_sig; % larger than nuclear ring
        sig_cyt(sig_cyt <= 0) = ring_sig/5; % smaller than background
        cytarea_fl = sum(sum(sig_cyt./ring_sig));
        all_obj.area_fl(j,i) = cytarea_fl + nuc_area;
        all_obj.volume_fl(j,i) = all_obj.area_fl(j,i)*all_obj.minAp_ax_mean(j,i);
        
        
        
        
    end%if for size
    
end
%------------------------end venus part---------------------------------

%---------------------start whi5-mcherry part---------------------------

c_num   ='c3';%mc used as whi5 here
[im_name] = get_image_name(prefix,pos_num,suffix,i, c_num,type,suffix2,numbM);
IG=imread(im_name);
%==============================
IG=IG(:,:,1);
[x_size,y_size]=size(IG);
all_obj.med_backgr_w5(1,i)=median(median(IG));

j=cell_no;
if cell_exists(j,1)==1
    put_nuc=uint8(Lfinal==j).*IG;
    all_obj.max_whi5_int(j,i)=max(max(put_nuc));
    
    p_nuc=bwmorph(bwmorph(bwareaopen(bwmorph((put_nuc>gauss_peak_cutoff.*all_obj.max_nucl_int(j,i)),'clean',inf),5,4),'close',1),'erode',0);
    %p_cyt=(Lfinal==j).*~p_nuc;
    p_cyt=bwmorph((Lfinal==j).*bwmorph(~p_nuc,'thicken',1),'erode',2);
    nuc_area=sum(sum(p_nuc));
    cyt_area=sum(sum(p_cyt));
    all_obj.mean_nuc_int_per_area(j,i)=sum(sum(double(p_nuc).*double(put_nuc)))./nuc_area;
    all_obj.tot_nucl_area(j,i)=nuc_area;
    all_obj.tot_cyt_area(j,i)=cyt_area;
    all_obj.cyt_whi5(j,i)=sum(sum(double(p_cyt).*double(put_whi5)))./cyt_area;
    all_obj.nuc_whi5(j,i)=sum(sum(double(p_nuc).*double(put_whi5)))./nuc_area;
    

end


%--------------------------------------------------------------------------
%--------------next kind of fit with whi5 gaussian....------------------
%------------------------------------------------------------------------

j=cell_no;
if cell_exists(j,1)==1
    %first identify the nucleus
    
    %put_nuc =uint8(Lfinal==j).*IN;
    put_whi5=uint8(Lfinal==j).*IG;
    
    %all_obj.max_nucl_int(j,i)=max(max(put_nuc));
    
    
    %-----gaussian fit--------------------
    h = fspecial('average', 5);
    I_gf=  imfilter(IG,h,'replicate');
    % im_min=min(min(I_gf));%(FLZ);
    % Ignew=I_gf-im_min;
    % im_max=double(max(max(Ignew)));
    % Max_f=255./im_max;
    % Ignew=uint8(Max_f.*double(Ignew));
    Ignew=I_gf;
    
    %nuc_tmp=IG.*uint8(Lcells);
    nuc_tmp=Ignew.*uint8(Lfinal==j);
    
    Amp=max(max(nuc_tmp));
    s_f=regionprops(bwlabel(nuc_tmp==Amp),'Centroid');
    xtmpcenter=s_f(1).Centroid(1,1);
    ytmpcenter=s_f(1).Centroid(1,2);
    wind_size=10;
    %[X, Y] = meshgrid(xtmpcenter-wind_size:1:xtmpcenter+wind_size, ytmpcenter-wind_size:1:ytmpcenter+wind_size);
    
    xtmpcoordhere=max(1,xtmpcenter-wind_size):min(y_size,xtmpcenter+wind_size);
    ytmpcoordhere=max(1,ytmpcenter-wind_size):min(x_size,ytmpcenter+wind_size);
    %if ~(min(xtmpcoordhere)==1 || max(xtmpcoordhere)==y_size || min(ytmpcoordhere)==1 || max(ytmpcoordhere)==x_size)
    [X, Y] = meshgrid(xtmpcoordhere,ytmpcoordhere);
    
    sX=size(X);
    sY=size(Y);
    ok_size=wind_size.*2+1;
    %[sX(1,1) sX(1,2) sY(1,1) sY(1,2)]
    if sX(1,1)==ok_size && sX(1,2)==ok_size && sY(1,1)==ok_size && sY(1,2)==ok_size
        %[X, Y] = meshgrid(xtmpcoordhere,ytmpcoordhere);
        
        
        In_part=double(nuc_tmp(round(diag(Y)),round(diag(X))));
        
        edge_med=median([In_part(1,:) In_part(end,:) In_part(:,1)' In_part(:,end)']);
        x0=xtmpcenter;
        y0=ytmpcenter;
        mu=[x0,y0];
        best_fit_score=1e9;
        best_fit_no=[0 0];
   
        for jj=1:2:50
            for jj2=1:5:50
                limit_h=round(sqrt(jj.*jj2));
                for ii=-(limit_h-1):2:(limit_h-1)
                    Sigma = [jj -ii; -ii jj2];
                    F = mvnpdf([X(:) Y(:)],mu,Sigma);
                    F = reshape(F,length(X(:,1)),length(Y(:,1)));
                    Fmax=max(max(F));
                    Z=(( ((double(Amp)-edge_med)./Fmax)).*F+edge_med);
                    Z=Z.*(In_part>0);
                    tmp_score=sum(sum(abs((Z-In_part))));
                    if tmp_score<best_fit_score
                        best_fit_no=Sigma;%[ii jj];
                        best_fit_score=tmp_score;
                    end
                end
            end
        end
        
        %----------------- get the best solution------------------------
        %Sigma = [best_fit_no(1,2) -best_fit_no(1,1); -best_fit_no(1,1) best_fit_no(1,2)];
        %F = mvnpdf([X(:) Y(:)],mu,Sigma);
        F = mvnpdf([X(:) Y(:)],mu,best_fit_no);
        F = reshape(F,length(X(:,1)),length(Y(:,1)));
        Fmax=max(max(F));
        Z=(( ((double(Amp)-edge_med)./Fmax)).*F+edge_med);
        Z=Z.*(In_part>0);
        %imshow(Z>(Amp.*gauss_peak_cutoff))
        %---------------------------------------------------------------
        p_nuc=zeros(x_size,y_size);
        p_nuc(round(diag(Y)),round(diag(X)))=p_nuc(round(diag(Y)),round(diag(X)))+(Z>(Amp.*gauss_peak_cutoff)); %was *0.75
        
        %---------------------------------------------------------------------
        p_cyt=bwmorph((Lfinal==j).*bwmorph(~p_nuc,'thicken',1),'erode',2);
        nuc_area=sum(sum(p_nuc));
        cyt_area=sum(sum(p_cyt));
        %  all_obj.mean_nuc_int_per_area_w5(j,i)=sum(sum(double(p_nuc).*double(put_nuc)))./nuc_area;
        all_obj.tot_nucl_area_w5(j,i)=nuc_area;
        all_obj.tot_cyt_area_w5(j,i) =cyt_area;
        all_obj.cyt_whi5R(j,i)=sum(sum(double(p_cyt).*double(put_whi5)))./cyt_area;
        all_obj.nuc_whi5R(j,i)=sum(sum(double(p_nuc).*double(put_whi5)))./nuc_area;
        all_obj.nuclLR(:,:,i)       = all_obj.nuclLR(:,:,i) +p_nuc;% +j.* p_nuc;
       
    end%if for size
    
    
    
    %--------------------------------------------------------------------------
    tmpX=double(all_obj.cells(:,:,i)==j).*double(IN);
    all_obj.abs_Far1(j,i)=sum(tmpX(:))-all_obj.med_backgr_Far1(i).*sum(sum(double((all_obj.cells(:,:,i)==j)>0)));
    if sum(sum(tmpX>0))>100
        all_obj.Far1_appr_conc(j,i)=all_obj.abs_Far1(j,i)./ (sum(sum(tmpX>0)).^(1.5));
    end
    all_obj.nuclear_area_F1(j,i)  =sum(sum((all_obj.cells(:,:,i)==j).*all_obj.nuclLFar1(:,:,i)));
    all_obj.nuclear_area_W5(j,i)  =sum(sum((all_obj.cells(:,:,i)==j).*all_obj.nuclLR(:,:,i)));
    all_obj.nuclear_area_FoW(j,i) =sum(sum((all_obj.cells(:,:,i)==j).*((all_obj.nuclLFar1(:,:,i)).*(all_obj.nuclLR(:,:,i)))));
   
   
    %total Far1 nuclear int accordingly
    [im_name] = get_image_name(prefix,pos_num,suffix,i,c_num,type,suffix2,numbM);INt=imread(im_name);
    tmpX=double(all_obj.cells(:,:,i)==j).*double(INt);
    
    cell_mask=uint8(all_obj.cells(:,:,i)==j);
  %  ntmp1 = INt.*cell_mask.*uint8(all_obj.nuclLFar1(:,:,i));
  %  ntmp2 = INt.*cell_mask;
  % all_obj.total_nuclear_Far1_F1(j,i)=sum(ntmp1(:))-all_obj.med_backgr_Far1(i).*uint8(sum(cell_mask(:)>0));
  %  all_obj.abs_Far1_wo_med(j,i)=sum(ntmp2(:));
   
    current_nuclear_mask_Far1=double(all_obj.nuclLFar1(:,:,i)).*double(all_obj.cells(:,:,i)==j);   
    current_nuclear_Far1=current_nuclear_mask_Far1.*tmpX;
   
   % sum(current_nuclear_mask_Far1(:)).* all_obj.med_backgr_Far1(1,i) %total background in nucleus
   % sum(current_nuclear_mask_Far1(:)).* all_obj.med_backgr_Far1(1,i) %total background in cell
   
    all_obj.abs_Far1_wo_med(j,i)=sum(tmpX(:))-(sum(cell_mask(:)).*all_obj.med_backgr_Far1(1,i));
  %  all_obj.total_nuclear_Far1_F1(j,i)=sum(current_nuclear_Far1(:));
    
  % all_obj.total_nuclear_Far1_F1(j,i)= sum(sum(tmpX.*double(all_obj.nuclLFar1(:,:,i))));
    all_obj.total_nuclear_Far1_W5(j,i)=sum(sum(tmpX.*double(all_obj.nuclLR(:,:,i))));
    all_obj.total_nuclear_Far1_B(j,i) =sum(sum(tmpX.*((all_obj.nuclLFar1(:,:,i)).*(all_obj.nuclLR(:,:,i)))));
   
    all_obj.total_nuclear_Far1_F1(j,i)=sum(current_nuclear_Far1(:))-(all_obj.nuclear_area_F1(j,i).*all_obj.med_backgr_Far1(1,i));
     
     
    
   
     
     
    
end

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% clear all
% close all
% 
% %program to fix annoying fucked up fl fcn, written june 26th 2013
% ptmp=[22 24:33];
% 
% for posh=ptmp
%     current_pos=['AD_18_5_240Gto0D_jun15_2013_pos_no_' num2str(posh) '_1'];
%     disp(current_pos)
%     load(current_pos)
% 
%     no_cells=length(cell_exists(:,1));
%     no_cells_here=no_cells;
%     for j=1:no_cells     
%     for i=numbM:-1:1
%         
%         if cell_exists(j,2)<i          
%                  new_c_Image2=double(all_obj.cells(:,:,i)==j);       
%                  [vol_c,area_c,maj_appr_axis,min_appr_axis]=get_area3(new_c_Image2,i,j); %i=cell number            
%                  all_obj.volume(j,i)    = vol_c;           
%                  all_obj.area(j,i)      =  area_c;    
%                  all_obj.proj_area(j,i) = sum(new_c_Image2(:)>0);    
%                  all_obj.majAp_ax(j,i)= maj_appr_axis;     
%                  all_obj.minAp_ax_med(j,i)= median(min_appr_axis(min_appr_axis>0));     
%                  all_obj.minAp_ax_sd(j,i)  =sqrt(var(min_appr_axis(min_appr_axis>0)));    
%                  all_obj.minAp_ax_len(j,i)  =length(min_appr_axis(min_appr_axis>0));    
%       
%                  %--------------------------------------------------------------------------------
%               all_obj.nuclear_area_F1(j,i)  =sum(sum((all_obj.cells(:,:,i)==j).*all_obj.nuclLFar1(:,:,i)));
%     all_obj.nuclear_area_W5(j,i)  =sum(sum((all_obj.cells(:,:,i)==j).*all_obj.nuclLR(:,:,i)));
%     all_obj.nuclear_area_FoW(j,i) =sum(sum((all_obj.cells(:,:,i)==j).*((all_obj.nuclLFar1(:,:,i)).*(all_obj.nuclLR(:,:,i)))));
%     %total Far1 nuclear int accordingly
%     [im_name] = get_image_name(prefix,pos_num,suffix,i,'c3',type,suffix2,numbM);INt=imread(im_name);
%     tmpX=double(all_obj.cells(:,:,i)==j).*double(INt);
%     
%     cell_mask=uint8(all_obj.cells(:,:,i)==j);
%     current_nuclear_mask_Far1=double(all_obj.nuclLFar1(:,:,i)).*double(all_obj.cells(:,:,i)==j);   
%     current_nuclear_Far1=current_nuclear_mask_Far1.*tmpX;
%     all_obj.abs_Far1_wo_med(j,i)=sum(tmpX(:))-(sum(cell_mask(:)).*all_obj.med_backgr_Far1(1,i));
%     
%     
%    % all_obj.total_nuclear_Far1_F1(j,i)=sum(current_nuclear_Far1(:));
%     all_obj.total_nuclear_Far1_F1(j,i)=sum(current_nuclear_Far1(:))-(all_obj.nuclear_area_F1(j,i).*all_obj.med_backgr_Far1(1,i));
%    
%     
%     all_obj.total_nuclear_Far1_W5(j,i)=sum(sum(tmpX.*double(all_obj.nuclLR(:,:,i))));
%     all_obj.total_nuclear_Far1_B(j,i) =sum(sum(tmpX.*((all_obj.nuclLFar1(:,:,i)).*(all_obj.nuclLR(:,:,i)))));
% 
% %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%   %   sum(current_nuclear_mask_Far1(:))
%             
%         end %cell exists loop
%     end %time
%     end %no cells
%         name1=[prefix '_pos_no_'  num2str(pos_num) '_X2'];
%         save(name1,    'all_obj',    'prefix',    'suffix2',    'pos_num',    'suffix',    'numbM',...
%             'type',    'LcellsO',    'cell_exists',    'GFP_subtraction_factor',     'higher_threshold',...
%             'lower_threshold',     'max_area_increase_per_tp',     'max_size_vs_largest_cell',...
%             'no_cells_here',     'phase_subtraction_factor',     'threshold_increase_factor')
%  
% end %no positions
% 

     
      
                    
               
             
