function [all_obj]=get_fl_N1(all_obj,i,prefix,suffix2,pos_num,suffix,numbM,type,Lfinal,cell_exists,cell_no)
%this version fixed for better gaussian nucl fitting

%cell_no=1
%i=10;
%Lfinal=Lcells;

% X all_obj.nuclL  =single(zeros(x_size,y_size,numbM));
%  all_obj.nuclLR  =single(zeros(x_size,y_size,numbM));
%  all_obj.nuclLG  =single(zeros(x_size,y_size,numbM));

%  c_num   ='c1';%phase
%  [im_name] = get_image_name(prefix,pos_num,suffix,i, c_num,type,suffix2,numbM);
%  I=imread(im_name);
c_num   ='c3';% orange whi5
[im_name] = get_image_name(prefix,pos_num,suffix,i, c_num,type,suffix2,numbM);
IN=imread(im_name);
c_num   ='c2';% gfp
[im_name] = get_image_name(prefix,pos_num,suffix,i, c_num,type,suffix2,numbM);
IG=imread(im_name);
[x_size,y_size]=size(IG);

all_obj.med_backgr_GFP(1,i)  = median(median(IG));
all_obj.med_backgr_mkok(1,i) = median(median(IN));
%for j=1:max(max(Lfinal));
j=cell_no;
if cell_exists(j,1)==1
    %first identify the nucleus
    % assume whi5 is the nuclear marker for this part

    put_far1 = uint8(Lfinal==j).*IG;
    put_nuc  = uint8(Lfinal==j).*IN;
    all_obj.max_nucl_int(j,i)=max(max(put_nuc));
    all_obj.max_far1_int(j,i)=max(max(put_far1));
    % areas in which more than 50% of the max. intensity of nuclear marker
    % is present is called nuclear area
    p_nuc=bwmorph(bwmorph(bwareaopen(bwmorph((put_nuc>0.5.*all_obj.max_nucl_int(j,i)),'clean',inf),5,4),'close',1),'erode',0);
    p_cyt=bwmorph((Lfinal==j).*bwmorph(~p_nuc,'thicken',1),'erode',2); % cytoplasmic if not nuclear
    nuc_area=sum(sum(p_nuc));
    cyt_area=sum(sum(p_cyt));
    all_obj.mean_nuc_int_per_area(j,i)=sum(sum(double(p_nuc).*double(put_nuc)))./nuc_area;
    all_obj.tot_nucl_area(j,i)=nuc_area;
    all_obj.tot_cyt_area(j,i)=cyt_area;
    all_obj.cyt_far1(j,i)=sum(sum(double(p_cyt).*double(put_far1)))./cyt_area;
    all_obj.nuc_far1(j,i)=sum(sum(double(p_nuc).*double(put_far1)))./nuc_area;
    all_obj.nuclL(:,:,i)  = all_obj.nuclL(:,:,i) + p_nuc;
    all_obj.mean_cyt_far1(j,i)=sum(sum(put_far1))./sum(sum(put_far1>0));
end

%------------------------------------------------------------------------
j=cell_no;
if cell_exists(j,1)==1
    %first identify the nucleus (assuming it is whi5)
    put_nuc =uint8(Lfinal==j).*IN;
    put_far1=uint8(Lfinal==j).*IG;
    all_obj.max_nucl_int(j,i)=max(max(put_nuc));
    %-----gaussian fit--------------------
    h = fspecial('average', 5);
    I_nf=  imfilter(IN,h,'replicate');
    im_min=min(min(I_nf));%(FLZ);
    Innew=I_nf-im_min;
    im_max=double(max(max(Innew)));
    Max_f=255./im_max;
    Innew=uint8(Max_f.*double(Innew));
    
    %nuc_tmp=IN.*uint8(Lcells);
    nuc_tmp=Innew.*uint8(Lfinal==j);
    
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
                    % imagesc(F);pause(0.01)
                    %  surf(Z-In_part);title(num2str(tmp_score));pause(0.001)
                    %  imagesc(abs(Z-In_part));pause(0.001)
                end
            end
        end
        
        %----------------- get the best solution------------------------
        %   Sigma = [best_fit_no(1,2) -best_fit_no(1,1); -best_fit_no(1,1) best_fit_no(1,2)];
        %   F = mvnpdf([X(:) Y(:)],mu,Sigma);
        F = mvnpdf([X(:) Y(:)],mu, best_fit_no);
        F = reshape(F,length(X(:,1)),length(Y(:,1)));
        Fmax=max(max(F));
        Z=(( ((double(Amp)-edge_med)./Fmax)).*F+edge_med);
        Z=Z.*(In_part>0);
        
        %figure(9);imagesc(double(Z)+double((Lfinal==j)>0));%imshow(Z>(Amp.*0.5))
        %---------------------------------------------------------------
        p_nuc=zeros(x_size,y_size);
        p_nuc(round(diag(Y)),round(diag(X)))=p_nuc(round(diag(Y)),round(diag(X)))+(Z>(Amp.*0.5));
        
        %---------------------------------------------------------------------
        p_cyt=bwmorph((Lfinal==j).*bwmorph(~p_nuc,'thicken',1),'erode',2);
        nuc_area=sum(sum(p_nuc));
        cyt_area=sum(sum(p_cyt));
        all_obj.mean_nuc_int_per_areaG(j,i)=sum(sum(double(p_nuc).*double(put_nuc)))./nuc_area;
        all_obj.tot_nucl_areaG(j,i)=nuc_area;
        all_obj.tot_cyt_areaG(j,i)=cyt_area;
        all_obj.cyt_far1G(j,i)=sum(sum(double(p_cyt).*double(put_far1)))./cyt_area;
        all_obj.nuc_far1G(j,i)=sum(sum(double(p_nuc).*double(put_far1)))./nuc_area;
        all_obj.nuclLG(:,:,i)  = all_obj.nuclLG(:,:,i) + p_nuc;%+j.* p_nuc;
    end%if for size
    
end

%--------------------------------------------------------------------------
%--------------next kind of fit not using the nucleus....------------------
%------------------------------------------------------------------------

j=cell_no;
if cell_exists(j,1)==1
    %first identify the nucleus
    
    %put_nuc =uint8(Lfinal==j).*IN;
    put_far1=uint8(Lfinal==j).*IG;
    
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
                    % imagesc(F);pause(0.01)
                    %  surf(Z-In_part);title(num2str(tmp_score));pause(0.001)
                    %  imagesc(abs(Z-In_part));pause(0.001)
                end
            end
        end
        %         for jj=1:1:50
        %             for ii=-(jj-1):1:(jj-1)
        %                 Sigma = [jj -ii; -ii jj];
        %                 F = mvnpdf([X(:) Y(:)],mu,Sigma);
        %                 F = reshape(F,length(X(:,1)),length(Y(:,1)));
        %                  Fmax=max(max(F));
        %                  Z=(( ((double(Amp)-edge_med)./Fmax)).*F+edge_med);
        %                  Z=Z.*(In_part>0);
        %                  tmp_score=sum(sum(abs((Z-In_part))));
        %                  %surf(Z-In_part);title(num2str(tmp_score));pause(0.001)
        %                  % imagesc(abs(Z-In_part));pause(0.001)
        %                  if tmp_score<best_fit_score
        %                      best_fit_no=[ii jj];
        %                      best_fit_score=tmp_score;
        %                  end
        %             end
        %         end
        %----------------- get the best solution------------------------
        % Sigma = [best_fit_no(1,2) -best_fit_no(1,1); -best_fit_no(1,1) best_fit_no(1,2)];
        % F = mvnpdf([X(:) Y(:)],mu,Sigma);
        F = mvnpdf([X(:) Y(:)],mu, best_fit_no);
        F = reshape(F,length(X(:,1)),length(Y(:,1)));
        Fmax=max(max(F));
        Z=(( ((double(Amp)-edge_med)./Fmax)).*F+edge_med);
        Z=Z.*(In_part>0);
        
        %imshow(Z>(Amp.*0.75))
        %---------------------------------------------------------------
        p_nuc=zeros(x_size,y_size);
        p_nuc(round(diag(Y)),round(diag(X)))=p_nuc(round(diag(Y)),round(diag(X)))+(Z>(Amp.*0.75));
        
        %---------------------------------------------------------------------
        p_cyt=bwmorph((Lfinal==j).*bwmorph(~p_nuc,'thicken',1),'erode',2);
        nuc_area=sum(sum(p_nuc));
        cyt_area=sum(sum(p_cyt));
        all_obj.mean_nuc_int_per_areaR(j,i)=sum(sum(double(p_nuc).*double(put_nuc)))./nuc_area;
        all_obj.tot_nucl_areaR(j,i)=nuc_area;
        all_obj.tot_cyt_areaR(j,i)=cyt_area;
        all_obj.cyt_far1R(j,i)=sum(sum(double(p_cyt).*double(put_far1)))./cyt_area;
        all_obj.nuc_far1R(j,i)=sum(sum(double(p_nuc).*double(put_far1)))./nuc_area;
        all_obj.nuclLR(:,:,i)  = all_obj.nuclLR(:,:,i) + p_nuc;%+j.* p_nuc;
    end%if for size
    
    
    
    %--------------------------------------------------------------------------
    % far1 concentration and other data
    tmpX=double(all_obj.cells(:,:,i)==j).*double(IN);
    all_obj.abs_far1(j,i)=sum(tmpX(:))-all_obj.med_backgr_GFP(i).*sum(sum(double((all_obj.cells(:,:,i)==j)>0)));
    if sum(sum(tmpX>0))>100
        all_obj.far1_appr_conc(j,i)=all_obj.abs_far1(j,i)./ (sum(sum(tmpX>0)).^(1.5));
    end
    all_obj.nuclear_area_far1(j,i)  =sum(sum((all_obj.cells(:,:,i)==j).*all_obj.nuclLR(:,:,i)));
    all_obj.nuclear_area_whi5(j,i)  =sum(sum((all_obj.cells(:,:,i)==j).*all_obj.nuclLG(:,:,i)));
    all_obj.nuclear_area_conv(j,i) =sum(sum((all_obj.cells(:,:,i)==j).*((all_obj.nuclLR(:,:,i)).*(all_obj.nuclLG(:,:,i)))));
   
   
    %total Far1 nuclear int accordingly
    [im_name] = get_image_name(prefix,pos_num,suffix,i,c_num,type,suffix2,numbM);INt=imread(im_name);
    tmpX=double(all_obj.cells(:,:,i)==j).*double(INt);
    
    cell_mask=uint8(all_obj.cells(:,:,i)==j);
  %  ntmp1 = INt.*cell_mask.*uint8(all_obj.nuclLFar1(:,:,i));
  %  ntmp2 = INt.*cell_mask;
  % all_obj.total_nuclear_Far1_F1(j,i)=sum(ntmp1(:))-all_obj.med_backgr_Far1(i).*uint8(sum(cell_mask(:)>0));
  %  all_obj.abs_Far1_wo_med(j,i)=sum(ntmp2(:));
   
    current_nuclear_mask_far1=double(all_obj.nuclLR(:,:,i)).*double(all_obj.cells(:,:,i)==j);   
    current_nuclear_far1=current_nuclear_mask_far1.*tmpX;
   
   % sum(current_nuclear_mask_Far1(:)).* all_obj.med_backgr_Far1(1,i) %total background in nucleus
   % sum(current_nuclear_mask_Far1(:)).* all_obj.med_backgr_Far1(1,i) %total background in cell
   
    all_obj.abs_Far1_wo_med(j,i)=sum(tmpX(:))-(sum(cell_mask(:)).*all_obj.med_backgr_GFP(1,i));
  %  all_obj.total_nuclear_Far1_F1(j,i)=sum(current_nuclear_Far1(:));
    
  % all_obj.total_nuclear_Far1_F1(j,i)= sum(sum(tmpX.*double(all_obj.nuclLFar1(:,:,i))));
    all_obj.total_nuclear_far1_G(j,i)=sum(sum(tmpX.*double(all_obj.nuclLG(:,:,i))));
    all_obj.total_nuclear_far1_conv(j,i) =sum(sum(tmpX.*((all_obj.nuclLR(:,:,i)).*(all_obj.nuclLG(:,:,i)))));
   
    all_obj.total_nuclear_far1_R(j,i)=sum(current_nuclear_far1(:))-(all_obj.nuclear_area_far1(j,i).*all_obj.med_backgr_GFP(1,i));
     
     
end



% for jj=1:2:50
%     for jj2=1:5:50
%         limit_h=round(sqrt(jj.*jj2));
%         for ii=-(limit_h-1):2:(limit_h-1)
%
%             Sigma = [jj -ii; -ii jj2];
%
%             F = mvnpdf([X(:) Y(:)],mu,Sigma);
%             F = reshape(F,length(X(:,1)),length(Y(:,1)));
%
%             % imagesc(F);pause(0.01)
%             %  surf(Z-In_part);title(num2str(tmp_score));pause(0.001)
%             %  imagesc(abs(Z-In_part));pause(0.001)
%         end
%     end
% end