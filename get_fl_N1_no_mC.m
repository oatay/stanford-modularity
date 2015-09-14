function [all_obj]=get_fl_N1_no_mC(all_obj,i,prefix,suffix2,pos_num,suffix,numbM,type,Lfinal,cell_exists,cell_no)
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
gauss_peak_cutoff = 0.85;
c_num   ='c4';%mC-nucl
[im_name] = get_image_name(prefix,pos_num,suffix,i, c_num,type,suffix2,numbM);
IN=imread(im_name);
c_num   ='c4';%GFP w5
[im_name] = get_image_name(prefix,pos_num,suffix,i, c_num,type,suffix2,numbM);
IG=imread(im_name);
[x_size,y_size]=size(IG);

all_obj.med_backgr_gfp(1,i)=median(median(IG));
%for j=1:max(max(Lfinal));
j=cell_no;
if cell_exists(j,1)==1
    %first identify the nucleus
    put_nuc =uint8(Lfinal==j).*IN;
    put_whi5=uint8(Lfinal==j).*IG;
    all_obj.max_nucl_int(j,i)=max(max(put_nuc));
    all_obj.max_whi5_int(j,i)=max(max(put_whi5));
    p_nuc=bwmorph(bwmorph(bwareaopen(bwmorph((put_nuc>gauss_peak_cutoff.*all_obj.max_nucl_int(j,i)),'clean',inf),5,4),'close',1),'erode',0);
    %  p_cyt=(Lfinal==j).*~p_nuc;
    p_cyt=bwmorph((Lfinal==j).*bwmorph(~p_nuc,'thicken',1),'erode',2);
    nuc_area=sum(sum(p_nuc));
    cyt_area=sum(sum(p_cyt));
    all_obj.mean_nuc_int_per_area(j,i)=sum(sum(double(p_nuc).*double(put_nuc)))./nuc_area;
    all_obj.tot_nucl_area(j,i)=nuc_area;
    all_obj.tot_cyt_area(j,i)=cyt_area;
    all_obj.cyt_whi5(j,i)=sum(sum(double(p_cyt).*double(put_whi5)))./cyt_area;
    all_obj.nuc_whi5(j,i)=sum(sum(double(p_nuc).*double(put_whi5)))./nuc_area;
    all_obj.nuclL(:,:,i)  = all_obj.nuclL(:,:,i) + p_nuc;%+j.* p_nuc;
    %new code as of nov7 2011
    all_obj.mean_cytGFP(j,i)=sum(sum(put_whi5))./sum(sum(put_whi5>0));  
    all_obj.tot_cytGFP(j,i)= sum(sum(put_whi5));   
end

%------------------------------------------------------------------------
j=cell_no;
if cell_exists(j,1)==1
    %first identify the nucleus
    put_nuc =uint8(Lfinal==j).*IN;
    put_whi5=uint8(Lfinal==j).*IG;
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
        p_nuc(round(diag(Y)),round(diag(X)))=p_nuc(round(diag(Y)),round(diag(X)))+(Z>(Amp.*gauss_peak_cutoff));
        
        %---------------------------------------------------------------------
        p_cyt=bwmorph((Lfinal==j).*bwmorph(~p_nuc,'thicken',1),'erode',2);
        nuc_area=sum(sum(p_nuc));
        cyt_area=sum(sum(p_cyt));
        all_obj.mean_nuc_int_per_areaG(j,i)=sum(sum(double(p_nuc).*double(put_nuc)))./nuc_area;
        all_obj.tot_nucl_areaG(j,i)=nuc_area;
        all_obj.tot_cyt_areaG(j,i)=cyt_area;
        all_obj.cyt_whi5G(j,i)=sum(sum(double(p_cyt).*double(put_whi5)))./cyt_area;
        all_obj.nuc_whi5G(j,i)=sum(sum(double(p_nuc).*double(put_whi5)))./nuc_area;
        all_obj.nuclLG(:,:,i)  = all_obj.nuclLG(:,:,i) + p_nuc;%+j.* p_nuc;
        
        % ring around nucleus
        rad = 4;
        nuc=(all_obj.cells(:,:,i)==j).*(all_obj.nuclLG(:,:,i));
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

%--------------------------------------------------------------------------
%--------------next kind of fit not using the nucleus....------------------
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
        all_obj.cyt_whi5R(j,i)=sum(sum(double(p_cyt).*double(put_whi5)))./cyt_area;
        all_obj.nuc_whi5R(j,i)=sum(sum(double(p_nuc).*double(put_whi5)))./nuc_area;
        all_obj.nuclLR(:,:,i)  = all_obj.nuclLR(:,:,i) + p_nuc;%+j.* p_nuc;
    end%if for size
    
    all_obj.med_backgr_Cln3(1,i) =median(median(IN));
    tmpX=double(all_obj.cells(:,:,i)==j).*double(IN);
    all_obj.abs_Cln3(j,i)=sum(tmpX(:))-all_obj.med_backgr_Cln3(i).*sum(sum(double((all_obj.cells(:,:,i)==j)>0)));
    if sum(sum(tmpX>0))>100
        all_obj.Cln3_appr_conc(j,i)=all_obj.abs_Cln3(j,i)./ (sum(sum(tmpX>0)).^(1.5));
    end
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