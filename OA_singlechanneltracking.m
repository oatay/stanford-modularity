

%-------program which segments and tracks-------------------------------
clear all
pos_num_here=[2 3 4 5 7 8 10 11];%[90 91 95:97 99 101:107 109:112];

%---set all relevant thresholds---
max_size_vs_largest_cell =3; %how much larger a cell we allow than the largest cell [fraction]
max_area_increase_per_tp =0.15;%maximal allowed area increase before switching to higher threshold (was 0.075) % 0.03
higher_threshold         =0.5;
lower_threshold          =0.2; % was 0.2 % -0.3
threshold_increase_factor=0.15;
phase_subtraction_factor =4; % was 3.0 % 12
GFP_subtraction_factor   =0.15; 

for all_pos=pos_num_here
    %--------------------------------------------------------------------------
    %clear all
    tic
    prefix  ='OA_070815_11Abar1d_refurb_3nM6min';
    suffix2='1';
    
    pos_num =all_pos;
    suffix  ='_';%_p2z2
    numbM    = 140; %the (last)timepoint
    type    ='.tif';
    c_num   ='c1';%phase
    alpha_add_t=numbM-0;
    [im_name] = get_image_name(prefix,pos_num,suffix,alpha_add_t, c_num,type,suffix2,numbM);
    I=imread(im_name);
    c_num   ='c4';%mCitrine (all whi5 values are actually Cln3)
    [im_name] = get_image_name(prefix,pos_num,suffix,alpha_add_t, c_num,type,suffix2,numbM);
    IG=imread(im_name);
    h = fspecial('average', 5);
    Igf=  imfilter(IG,h,'replicate');
    im_min=min(min(Igf));%(FLZ);
    Ignew=Igf-im_min;
    im_max=double(max(max(Ignew)));
    Max_f=255./im_max;
    Ignew=uint8(Max_f.*double(Ignew));
    [x_size,y_size]=size(I);
    %----------------------------------------------------
    
    
    load(['tmpsave' prefix '_pos' num2str(pos_num)]);
    %-load('tmpsaveAD_whi52X2_w_136_30min_met_may12_2010_pos4')
    Ifin2=LcellsO>0;
    %-----------------------------------------------------
    close all
    Lcells=bwlabel(Ifin2);
    LcellsO=Lcells;%o for original
    cell_exists=ones((max(max(Lcells))),2);
    %----------------assign things----------------------------
    no_obj=max(max(Lcells));
    all_obj.med_backgr_gfp       =zeros(1,numbM);
    
    all_obj.max_nucl_int         =zeros(no_obj,numbM);
    all_obj.max_whi5_int         =zeros(no_obj,numbM);
    all_obj.mean_nuc_int_per_area=zeros(no_obj,numbM);
    all_obj.tot_nucl_area        =zeros(no_obj,numbM);
    all_obj.tot_cyt_area         =zeros(no_obj,numbM);
    all_obj.nuc_whi5absv         =zeros(no_obj,numbM);
    all_obj.nuc_whi5frac         =zeros(no_obj,numbM);
    all_obj.cyt_whi5             =zeros(no_obj,numbM);
    all_obj.nuc_whi5             =zeros(no_obj,numbM);
    all_obj.cells  =single(zeros(x_size,y_size,numbM));
    all_obj.nuclL  =logical(zeros(x_size,y_size,numbM));
    %for the gaussian fits...
    all_obj.mean_nuc_int_per_areaG=zeros(no_obj,numbM);
    all_obj.tot_nucl_areaG        =zeros(no_obj,numbM);
    all_obj.tot_cyt_areaG         =zeros(no_obj,numbM);
    all_obj.cyt_whi5G             =zeros(no_obj,numbM);
    all_obj.nuc_whi5G             =zeros(no_obj,numbM);
    all_obj.nuclLG  =logical(zeros(x_size,y_size,numbM));
    %for the gaussian GFP-fit...
    all_obj.mean_nuc_int_per_areaR=zeros(no_obj,numbM);
    all_obj.tot_nucl_areaR        =zeros(no_obj,numbM);
    all_obj.tot_cyt_areaR         =zeros(no_obj,numbM);
    all_obj.cyt_whi5R             =zeros(no_obj,numbM);
    all_obj.nuc_whi5R             =zeros(no_obj,numbM);
    all_obj.nuclLR  =logical(zeros(x_size,y_size,numbM));
    
    all_obj.med_nuc_GFP         = zeros(no_obj,numbM);
    all_obj.area_fl             = zeros(no_obj,numbM);
    all_obj.volume_fl           = zeros(no_obj,numbM);
    
    all_obj.volume_ax           = zeros(no_obj,numbM);
    all_obj.area_ax             = zeros(no_obj,numbM);
    all_obj.proj_area_ax        = zeros(no_obj,numbM);%projected area
    all_obj.majAp_ax            = zeros(no_obj,numbM);
    all_obj.minAp_ax_med        = zeros(no_obj,numbM);
    all_obj.minAp_ax_mean       = zeros(no_obj,numbM);
    all_obj.minAp_ax_sd         = zeros(no_obj,numbM);
    all_obj.minAp_ax_len        = zeros(no_obj,numbM);
    
    %---new parameter as of nov7 2011---
    all_obj.mean_cytGFP           =zeros(no_obj,numbM);
    all_obj.area                  =zeros(no_obj,numbM);
    all_obj.Cln3_appr_conc        =zeros(no_obj,numbM);
    all_obj.med_backgr_Cln3       =zeros(no_obj,numbM);
    all_obj.abs_Cln3              =zeros(no_obj,numbM);
    %---------------------------------------------


    Lcells=LcellsO;
    min_cell_size=10;%in pixels
    tmp_sizes=zeros(1,max(max((LcellsO))));for i5=1:max(max((LcellsO)));tmp_sizes(i5)=sum(sum(LcellsO==i5));end
    max_allowed_cell_size=max_size_vs_largest_cell.*max(tmp_sizes);
    cell_margin=12;
    [x_size,y_size]=size(I);
    %--------------------------------------
    
    %--------------------------------------
    for i=1:max(max(LcellsO))
        Lcells=LcellsO;%bwlabel(LcellsO==i);
        for c_time=alpha_add_t:-1:1
            if cell_exists(i,1)==1
                disp([i c_time])
                c_num   ='c1';%phase
                [im_name] = get_image_name(prefix,pos_num,suffix,c_time, c_num,type,suffix2,numbM);
                I=imread(im_name);
                %-------------------------------------------
                c_num   ='c4';%GFP
                [im_name] = get_image_name(prefix,pos_num,suffix,c_time, c_num,type,suffix2,numbM);
                IG=imread(im_name);
                h = fspecial('average', 5);
                Igf=  imfilter(IG,h,'replicate');
                im_min=min(min(Igf));%(FLZ);
                Ignew=Igf-im_min;
                im_max=double(max(max(Ignew)));
                Max_f=255./im_max;
                Ignew=uint8(Max_f.*double(Ignew));
                %--------------------------------------------
                s_f=regionprops(bwlabel(Lcells==i),'Centroid','BoundingBox');
                new_c_Image =zeros(x_size,y_size);
                new_c_Image2=zeros(x_size,y_size);
                bbox=round(s_f(1).BoundingBox);
                lower_x_limit=max(1,bbox(1)-cell_margin);
                upper_x_limit=min(y_size,bbox(1)+bbox(3)+cell_margin);
                lower_y_limit=max(1,bbox(2)-cell_margin);
                upper_y_limit=min(x_size,bbox(2)+bbox(4)+cell_margin);
                x_cn=lower_x_limit:upper_x_limit;
                y_cn=lower_y_limit:upper_y_limit;
                tmpI1=(I(y_cn,x_cn));
                tmpI2=(Lcells==i);%was i
                tmpI2=tmpI2(y_cn,x_cn);
                [Ires] = rescale_im09(tmpI1,-.5,3);%
                %------------------------------------------------
                newI2=tmpI1<-100;
 
                exI=uint8(newI2);
                exI(1,1:end)=1;exI(end,1:end)=1;exI(1:end,1)=1;exI(1:end,end)=1; %added march 26th 2012
                for j=1:256
                    bw=tmpI1<(j); % was < !!!!!
                    bw2=bwmorph(bw,'thicken',1);
                    D = bwdist(~bw2);%original
                    D = -D;
                    D(~bw2) = -Inf; %note it is not bw2 but bw!
                    L = watershed(D);
                %    imagesc(L);pause
                %    intr_segm=setxor([0:max(max(L))],[unique(L.*exI)']);%added march 26th 2012
                    intr_segm=setxor(uint8([0:max(max(L))]),[unique(uint8(L).*uint8(exI))']);
                    %----add all relevant parts of L to the new object---------
                    newI=tmpI1<-100;
                    
                    for j2=intr_segm %1:max(max(L))
                        if sum(sum((L==j2).*tmpI2))./sum(sum((L==j2)))>0.4 ...
                                && sum(sum(((L==j2).*~bw)==1))./sum(sum((L==j2)))<0.33
                            newI=newI+(L==j2);
                        end
                    end
                    newI=bwmorph(imfill(bwmorph(newI,'close',1),'holes'),'thicken',1);
                    newI2=newI2+double(newI);
                end
                local_IG_mean=sum(sum(Ignew(y_cn,x_cn).*uint8(tmpI2)))./sum(sum(tmpI2));
                maxI2= max(max(newI2));%first maxima to scale against
                lin_im=newI2(:)';
                [q,w]=hist(lin_im,50);
                [qq,ww]=min(abs(0.90.*maxI2-w));
                maxI2=w(ww);
                %%--------- change of local IGmean factor as of apr13 2010----
                %newI2=newI2-(maxI2./256).*double(Ires) -0.15.*maxI2.*double(Ignew(y_cn,x_cn)<(0.65.*local_IG_mean));
                %%------------------------------------------------------------
                %--------- change of local weight for Ires as of oct4 2010----
                newI2=newI2-phase_subtraction_factor.*(maxI2./256).*double(Ires) -GFP_subtraction_factor.*maxI2.*double(Ignew(y_cn,x_cn)<(0.65.*local_IG_mean));
                %------------------------------------------------------------
              
                newI2=newI2-bwdist(tmpI2);
                maxI2= max(max(newI2));
                lin_im=newI2(:)';
                [q,w]=hist(lin_im,50);
                [qq,ww]=min(abs(0.90.*maxI2-w));
                maxI2=w(ww);
                
                t_mod=0;
                %         [high_thr,low_thr,t_mod,newI] = get_im_thr2(newI2,tmpI2,maxI2,numbM,c_time,t_mod);
              
                high_thr=higher_threshold;%0.55;%was 0.6
                low_thr =lower_threshold;%0.35; %was 0.4
                pcell=(newI2-tmpI2.*1)>(maxI2.*high_thr);
                prem =((newI2-tmpI2.*1)<(maxI2.*high_thr)).*((newI2-tmpI2.*1)>(maxI2.*low_thr));%was 0.2
                newI=(  pcell+prem.*tmpI2);
                newI=(imfill( bwmorph( bwareaopen(newI,10,4),'close',1),'holes'));
                
                
                %-------------------------------------------------
              
                area_increase=(sum(sum(newI))-sum(sum(tmpI2)))./sum(sum(tmpI2));
                if area_increase>max_area_increase_per_tp
                    disp(['higher thr enabled incr. frac. = ' num2str((sum(sum(newI))-sum(sum(tmpI2)))./sum(sum(tmpI2))) ''])
                    high_thr = higher_threshold + threshold_increase_factor;
                    low_thr  = lower_threshold  + threshold_increase_factor;
                    pcell=(newI2-tmpI2.*1)>(maxI2.*high_thr);
                    prem =((newI2-tmpI2.*1)<(maxI2.*high_thr)).*((newI2-tmpI2.*1)>(maxI2.*low_thr));%was 0.2
                    newI=(  pcell+prem.*tmpI2);
                    newI=(imfill( bwmorph( bwareaopen(newI,10,4),'close',1),'holes'));
                end
                
                %-new code by Oguzhan - 16/04/12
%                 area_decrease =((sum(sum(newI))-sum(sum(tmpI2)))./sum(sum(tmpI2)));disp(area_decrease)
%                 if abs(area_decrease)>max_area_increase_per_tp
%                     disp(['lower thr enabled incr. frac. = ' num2str((sum(sum(newI))-sum(sum(tmpI2)))./sum(sum(tmpI2))) ''])
%                     high_thr = higher_threshold - threshold_increase_factor;
%                     low_thr  = lower_threshold  - threshold_increase_factor;
%                     pcell=(newI2-tmpI2.*1)>(maxI2.*high_thr);
%                     prem =((newI2-tmpI2.*1)<(maxI2.*high_thr)).*((newI2-tmpI2.*1)>(maxI2.*low_thr));%was 0.2
%                     newI=(  pcell+prem.*tmpI2);
%                     newI=(imfill( bwmorph( bwareaopen(newI,10,4),'close',1),'holes'));
%                 end
                %-end new code by Oguzhan - 16/04/12-------------------------------------------------
                
                %------------------------------------------------------
%                                             figure(1);subplot(2,2,1);imagesc(newI2.*double(newI));title(c_time)
%                                             subplot(2,2,2);imagesc(newI2);title(c_time)
%                                             subplot(2,2,3);imagesc(pcell);title(c_time)
%                                             subplot(2,2,4);imagesc(prem);title(c_time);drawnow ;
                %      %------------------------------------------------
                %--------------remove multiple parts----------------------------
                
                
                Lny=bwlabel(newI);
                Lny2=newI<-1000;
                
                if sum(sum(newI))>min_cell_size
                    if max(max(Lny))==1  %all ok
                        Lny2=Lny2+(Lny==1);
                    else
                        maxindex=0;
                        tmparea=0;
                        for j2=1:max(max(Lny))
                            if sum(sum(Lny==j2))>tmparea
                                tmparea=sum(sum(Lny==j2));
                                maxindex=j2;
                            end
                        end
                        Lny2=Lny2+(Lny==maxindex);
                    end
                end
                %-----------------------------------------------------------
                new_c_Image(y_cn,x_cn)=new_c_Image(y_cn,x_cn)+i.*Lny2;
                %remove places where two cells share the same area:
                Ltmp=bwlabel(new_c_Image==i);
                s_f4=regionprops(Ltmp,'Area');
                if length(s_f4)>1
                    tmparea=0;
                    maindex=0;
                    for j4=1:length(s_f4)
                        if s_f4(j4).Area>tmparea
                            tmparea=s_f4(j4).Area;
                            maindex=j4;
                        end
                    end
                    new_c_Image2=new_c_Image2+i.*(Ltmp==maindex);
                else
                    new_c_Image2=new_c_Image2+i.*(new_c_Image==i);
                end
                
                if sum(sum(new_c_Image2==i))<min_cell_size ...
                        || sum(sum(uint8(newI).*Ignew(y_cn,x_cn)))./sum(sum(uint8(newI)))<mean(mean(Ignew)) ...
                        || sum(sum(new_c_Image2>0))>max_allowed_cell_size ...
                        || area_increase>0.3
                    cell_exists(i,1)=0;
                    cell_exists(i,2)=c_time;
                end
                all_obj.cells(:,:,c_time)=all_obj.cells(:,:,c_time)+new_c_Image2;

                [vol_c,area_c,maj_appr_axis,min_appr_axis]=get_area2(new_c_Image2,c_time,i); %i=cell number
                all_obj.volume_ax(i,c_time)    = vol_c;
                all_obj.area_ax(i,c_time)      =  area_c;
                all_obj.proj_area(i,c_time) = sum(new_c_Image2(:)>0);
                all_obj.majAp_ax(i,c_time)= maj_appr_axis;
                all_obj.minAp_ax_med(i,c_time)= median(min_appr_axis(min_appr_axis>0));
                all_obj.minAp_ax_mean(i,c_time)= mean(min_appr_axis(min_appr_axis>0));
                all_obj.minAp_ax_sd(i,c_time)  =sqrt(var(min_appr_axis(min_appr_axis>0)));
                all_obj.minAp_ax_len(i,c_time)  =length(min_appr_axis(min_appr_axis>0));
                

                all_obj.area(i,c_time)=sum(sum(double(new_c_Image2>0)));
                %------------------------------------------
                new_c_Image3=I<-1000;
                j5=i;
                if cell_exists(j5,1)==1
                    new_c_Image3=new_c_Image3+bwmorph(new_c_Image2==j5,'remove',inf);
                end
                %--------------------new code as of nov23 2009----------------------
                 Inew_here=uint8(zeros(x_size,y_size,3));
                 Inew_here(:,:,1)=0.5.*I+   255.*uint8(bwmorph(new_c_Image3,'remove',inf));
                 Inew_here(:,:,2)=0.5.*I+IG+  0.*uint8(bwmorph(new_c_Image3,'remove',inf));
                 Inew_here(:,:,3)=0.5.*I+     0.*uint8(bwmorph(new_c_Image3,'remove',inf));
                 if all_pos == pos_num_here(1)
                    figure(4);imshow(Inew_here);title(c_time);drawnow
                 end
                %--------------------end of new code--------------------------------
                if (cell_exists(i,1))>0
                    im_here=new_c_Image2==i;
                    s_ftmp=regionprops(bwlabel(im_here),'Centroid');
                    center_here=round(s_ftmp(1).Centroid);
                    text(center_here(1,1)+20, center_here(1,2),[num2str(i)],'BackgroundColor',[.7 .9 .7]);
                end
                hold off
                %--------------------------------------------
%                 %-------------------new code as of nov23 2009----------------------
%                 Inew_here=uint8(zeros(x_size,y_size,3));
%                 
%                 Inew_here(:,:,1)=0.5.*I+   255.*uint8(bwmorph((Lcells>0),'remove',inf));
%                 Inew_here(:,:,2)=0.5.*I+IG+  0.*uint8(bwmorph((Lcells>0),'remove',inf));
%                 Inew_here(:,:,3)=0.5.*I+     0.*uint8(bwmorph((Lcells>0),'remove',inf));
%                 figure(6);imshow(Inew_here);title(c_time);drawnow
                %--------------------end of new code----------------------
                %pause  
             
                
                Lcells=new_c_Image2;
                
                
                
                 % [all_obj]= get_fl_N1(all_obj,c_time,prefix,suffix2,pos_num,suffix,numbM,type,Lcells,cell_exists,i);
              %------------- for when no mC added nov2013   ------------------------------------------------------------------------
               disp('only single channel')
               %[all_obj]=get_fl_N1_no_mC(all_obj,c_time,prefix,suffix2,pos_num,suffix,numbM,type,Lcells,cell_exists,i);
               [all_obj]=get_fl_N1_no_mC(all_obj,c_time,prefix,suffix2,pos_num,suffix,numbM,type,Lcells,cell_exists,i);
               %-------------------------------------------------------------------------------------
               %   [all_obj]=get_fl_N1(all_obj,c_time,prefix,suffix2,pos_num,suffix,numbM,type,Lcells,cell_exists,i);
               if all_pos == pos_num_here(1)
                   figure(5);
                   plot(all_obj.nuc_whi5(i,:)-all_obj.cyt_whi5(i,:));hold on
                   plot(all_obj.nuc_whi5G(i,:)-all_obj.cyt_whi5G(i,:),'g');
                   plot(all_obj.nuc_whi5R(i,:)-all_obj.cyt_whi5R(i,:),'r');
                   plot(all_obj.mean_cytGFP(i,:),'m')
                   hold off
                   drawnow
               end
            end %end of cell exists loop
        end
    end %end of cell loop
    
    ttime=toc;
    disp(['total time: ' num2str(ttime./60) 'min'])
    name1=[prefix '_pos_no_'  num2str(pos_num) '_re_exp_vol'];
    save(name1, 'all_obj', 'prefix', 'suffix2', 'pos_num', 'suffix', 'numbM', 'type', 'LcellsO' ...
        ,'cell_exists','max_size_vs_largest_cell','max_area_increase_per_tp','higher_threshold','lower_threshold' ...
    ,'threshold_increase_factor','phase_subtraction_factor','GFP_subtraction_factor');   
    close all
    
    %--------------------------------------------------------------------------
end

%99999999999999999999999999999999999999999999999999999999999999999999999999




















