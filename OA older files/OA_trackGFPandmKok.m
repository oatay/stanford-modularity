
%---------------------------------------------------------
clear all

%--- setting of the constants ---
max_size_vs_largest_cell =3.0; %how much larger (fraction) pixel area a cell can be compared with the final image
max_area_increase_per_tp =0.15;%maximal allowed area increase before switching to higher threshold (was 0.075) % 0.03
higher_threshold         =0.4;  %segmentation threshold for regions that wasnt part of the cell in the previous time point
lower_threshold          =0.1;  %segmentation threshold for regions that was part of the cell in the previous time point
threshold_increase_factor=0.15;  %change of thresholds in case of large change in area
phase_subtraction_factor =4;  %factor that determines how much of the phase image we deduct from the segmentation
GFP_subtraction_factor   =0.15;  %(originally 0.5) factor determining how much of gfp-image we deduct, set to zero (0) if this feature is not desired
min_cell_size            =10;   %minimal cell size allowed in pixels
cell_margin              =12;   %margin around the cell from the previous segmentaiton, for current segmentation

pos_num_here=[18];
for all_pos=pos_num_here
    %--------------------------------------------------------------------------
    %clear all
    tic
    prefix  ='OA_123113_OA037andOA038_120nMmetexpt';
    suffix2='7';
    pos_num =all_pos;
    suffix  ='_';%_p2z2
    numbM   =99; %the (last)timepoint
    type    ='.tif';
    c_num   ='c1';%phase
    alpha_add_t=numbM-0;
    %read in phase, I
    [im_name] = get_image_name(prefix,pos_num,suffix,alpha_add_t, c_num,type,suffix2,numbM);
    I=imread(im_name);
    [x_size,y_size]=size(I);
    %----------------------------------------------------
    load(['tmpsave' prefix '_pos' num2str(pos_num)]);
    Ifin2=LcellsO>0;
    %-----------------------------------------------------
    close all
    %algorithm for finding the individual cells "no matter what"
    Lcells=bwlabel(Ifin2);
    LcellsO=Lcells;%o for original
    cell_exists=ones((max(max(Lcells))),2);
    %--------------------------------------------
    no_obj=max(max(Lcells));
    
    % all_obj data storage
    all_obj.cells                  = single(zeros(x_size,y_size,numbM)); % store segmented cells
    all_obj.uniquecells            = single(zeros(x_size,y_size,no_obj,numbM));
    
    % assuming whi5 as a nuclear marker (only true during arrest)
    all_obj.med_backgr_GFP         = zeros(1,numbM); % med background GFP (c2)
    all_obj.max_nucl_int           = zeros(no_obj,numbM); % max whi5-mkok intensity
    all_obj.mean_GFP_int_per_area  = zeros(no_obj,numbM); %
    % areas in which more than 50% of the max. intensity of nuclear marker is present is called nuclear area
    all_obj.nuclL                  = logical(zeros(x_size,y_size,numbM)); % nuclear areas based on whi5
    all_obj.mean_nuc_int_per_area  = zeros(no_obj,numbM); % whi5
    all_obj.tot_nucl_area          = zeros(no_obj,numbM); % whi5 area
    all_obj.tot_cyt_area           = zeros(no_obj,numbM); % non-whi5 area
    all_obj.cyt_far1               = zeros(no_obj,numbM); % total sum of intensities in non-whi5-area far1 divided by non-whi5 area
    all_obj.nuc_far1               = zeros(no_obj,numbM); % total sum of intensities in whi5-area far1 divided by whi5 area
    all_obj.med_backgr_mkok        = zeros(1,numbM); % med background mkok (c3)
    all_obj.max_far1_int           = zeros(no_obj,numbM); % max far1-gfp intensity

    
    % gaussian fitting whi5 to identify it as nucleus
    all_obj.mean_nuc_int_per_areaG = zeros(no_obj,numbM); % whi5
    all_obj.tot_nucl_areaG         = zeros(no_obj,numbM); % whi5 area
    all_obj.tot_cyt_areaG          = zeros(no_obj,numbM); % non-whi5 area
    all_obj.cyt_far1G              = zeros(no_obj,numbM); % total sum of intensities in non-whi5-area far1 divided by non-whi5 area
    all_obj.nuc_far1G              = zeros(no_obj,numbM); % total sum of intensities in whi5-area far1 divided by whi5 area
    all_obj.nuclLG                  = logical(zeros(x_size,y_size,numbM));
   
    % fit without using the nucleus (only fits to far1-mkok) -- may not be
    % very useful if the signal is not in nucleus
    all_obj.tot_nucl_areaR         = zeros(no_obj,numbM);
    all_obj.tot_cyt_areaR          = zeros(no_obj,numbM);
    all_obj.cyt_far1R                = zeros(no_obj,numbM);
    all_obj.nuc_far1R                = zeros(no_obj,numbM);
    all_obj.nuclLR                   = logical(zeros(x_size,y_size,numbM));

    % additional data for far1
    all_obj.abs_far1          =  zeros(no_obj,numbM); % median subtracted far1
    all_obj.abs_far1_wo_med   = zeros(no_obj,numbM);
    all_obj.far1_appr_conc    =  zeros(no_obj,numbM);
    
    % different far1 nuclear areas and total far1
    all_obj.nuclear_area_far1        = zeros(no_obj,numbM);
    all_obj.nuclear_area_whi5        = zeros(no_obj,numbM);
    all_obj.nuclear_area_conv        = zeros(no_obj,numbM);
    all_obj.total_nuclear_far1_R     = zeros(no_obj,numbM);
    all_obj.total_nuclear_far1_G     = zeros(no_obj,numbM);
    all_obj.total_nuclear_far1_conv  = zeros(no_obj,numbM);

    
    % data from get_area function
    all_obj.volume              = zeros(no_obj,numbM);
    all_obj.area                = zeros(no_obj,numbM);
    all_obj.proj_area           = zeros(no_obj,numbM);
    all_obj.majAp_ax            = zeros(no_obj,numbM);
    all_obj.minAp_ax            = zeros(no_obj,numbM);
    %---------------------------------------------
    
    %backward in time
    Lcells=LcellsO;
    
    %find largest cell to calculate max allowed cell size---
    tmp_sizes=zeros(1,max(max((LcellsO))));for i5=1:max(max((LcellsO)));tmp_sizes(i5)=sum(sum(LcellsO==i5));end
    max_allowed_cell_size=max_size_vs_largest_cell.*max(tmp_sizes);
    
    [x_size,y_size]=size(I);
    no_cells_here=max(max(LcellsO));
    %--------------------------------------
    for i=1:no_cells_here
        Lcells=LcellsO;%bwlabel(LcellsO==i);
        for c_time=alpha_add_t:-1:1
            if cell_exists(i,1)==1
                disp([i c_time])
                c_num   ='c1';%phase
                [im_name] = get_image_name(prefix,pos_num,suffix,c_time, c_num,type,suffix2,numbM);
                I=imread(im_name);
                %-------------------------------------------
                c_num   ='c3'; % orange
                [im_name] = get_image_name(prefix,pos_num,suffix,c_time, c_num,type,suffix2,numbM);
                IY=imread(im_name);
                %---------------------------
                c_num   ='c2';% GFP read
                [im_name] = get_image_name(prefix,pos_num,suffix,c_time, c_num,type,suffix2,numbM);
                IG=imread(im_name);
                % XX filtering
                h = fspecial('average', 5);
                Igf=  imfilter(IG,h,'replicate');
                im_min=min(min(Igf));
                Ignew=Igf-im_min;
                im_max=double(max(max(Ignew)));
                Max_f=255./im_max;
                Ignew=uint8(Max_f.*double(Ignew));
                %--------------------------------------------
                s_f=regionprops(bwlabel(Lcells==i),'Centroid','BoundingBox','Area'); %identify cell, label cell, calculate characteristics
                new_c_Image =zeros(x_size,y_size); %define new images
                new_c_Image2=zeros(x_size,y_size);
                
                %---get the cell and not an intersection
                tmpareas=zeros(1,length(s_f));
                for i5=1:length(s_f);tmpareas(i5)=s_f(i5).Area;end
                [val,pos]= max(tmpareas); ok_obj_no=pos;
                %---
                
                %retrieve the bounding box + margins for the current cell
                bbox=round(s_f(ok_obj_no).BoundingBox);
                lower_x_limit=max(1,bbox(1)-cell_margin);
                upper_x_limit=min(y_size,bbox(1)+bbox(3)+cell_margin);
                lower_y_limit=max(1,bbox(2)-cell_margin);
                upper_y_limit=min(x_size,bbox(2)+bbox(4)+cell_margin);
                x_cn=lower_x_limit:upper_x_limit;
                y_cn=lower_y_limit:upper_y_limit;
                
                %open images around the current cell phase and binary
                tmpI1=(I(y_cn,x_cn));
                tmpI2=(Lcells==i);%was i
                tmpI2=tmpI2(y_cn,x_cn);
                [Ires] = rescale_im09(tmpI1,-.5,3);%rescales the phase image
                %------------------------------------------------
                
                newI2=tmpI1<-100;  %make new image for iterative segmentation
                %calculate edges of image (to be excluded from segmentation)
                exI=uint8(newI2);exI(1,1:end)=1;exI(end,1:end)=1;exI(1:end,1)=1;exI(1:end,end)=1; 
                for j=1:256 %for all possible image intensities
                    %again watershed algorithm
                    bw=tmpI1<(j); 
                    bw2=bwmorph(bw,'thicken',1);
                    D = bwdist(~bw2);%original
                    D = -D;
                    D(~bw2) = -Inf; 
                    L = watershed(D);
                    %--end watershed
                    
                    intr_segm=setxor([0:max(max(L))],[unique(uint8(L).*uint8(exI))']);
                    %----add all relevant parts of L to the new object---------
                    newI=tmpI1<-100;
                    for j2=intr_segm %only go over parts that are not included in edges
                        if sum(sum((L==j2).*tmpI2))./sum(sum((L==j2)))>0.4 ...         %select pieces that are over 40% cell and less than 33% non-cell
                                && sum(sum(((L==j2).*~bw)==1))./sum(sum((L==j2)))<0.33
                            newI=newI+(L==j2);                                         %add all selected pieces to the same image
                        end
                    end
                    newI=bwmorph(imfill(bwmorph(newI,'close',1),'holes'),'thicken',1); %removes holes and such
                    newI2=newI2+double(newI);                                          %add to selected images from previous image                  
                % figure(3);imagesc(L);title(j);pause(0.05)
                end
                   
                local_IG_mean=sum(sum(Ignew(y_cn,x_cn).*uint8(tmpI2)))./sum(sum(tmpI2));%was off
                maxI2= max(max(newI2));%first maxima to scale against
                lin_im=newI2(:)';
                [q,w]=hist(lin_im,50);
                [qq,ww]=min(abs(0.90.*maxI2-w));
                maxI2=w(ww);
                
                newI2=newI2-phase_subtraction_factor.*(maxI2./256).*double(Ires) -GFP_subtraction_factor.*maxI2.*double(Ignew(y_cn,x_cn)<(0.85.*local_IG_mean));
                newI2=newI2-bwdist(tmpI2);
                
                
                maxI2= max(max(newI2));%first maxima to scale against
                lin_im=newI2(:)';
                [q,w]=hist(lin_im,50);
                [qq,ww]=min(abs(0.90.*maxI2-w));
                maxI2=w(ww);
    
                high_thr=higher_threshold; %threshold for regions not-previously designated as "cell"
                low_thr =lower_threshold;  %threshold for regions     previously designated as "cell"                
                %calculate which putative cell regions that corresponds to previously segmented areas and which that doesnt
                pcell=(newI2-tmpI2.*1)>(maxI2.*high_thr);
                prem =((newI2-tmpI2.*1)<(maxI2.*high_thr)).*((newI2-tmpI2.*1)>(maxI2.*low_thr));%was 0.2                
                %new segmented image
                newI=(  pcell+prem.*tmpI2);
                newI=(imfill( bwmorph( bwareaopen(newI,10,4),'close',1),'holes'));
                %----change segmentation if too large increase-decrease in area (if needed)---------------------------------------------                
                % if the area increases more than some specified fraction (see "thresholds" above) 
                % then increase the segmentation thresholds and recalculate the segmentated image
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
                %same as above but for sudden decreases in cell area
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
                %----end change segmentation if too large increase-decrease in area (if needed)--------------------------------------------- 
                
             
                for z1=1 %cosmetic loop
                    %dont change this code
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
                %-----------------end remove multiple parts------------------------------------------
                new_c_Image(y_cn,x_cn)=new_c_Image(y_cn,x_cn)+i.*Lny2;
                
                %----- remove places where two cells share the same area:
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
                %----- end remove places where two cells share the same area:
                end %z1 cosmetic loop
                %determine if the segmentation is ok, ie if it fulfills the following conditions: 
                %1. larger than minimal cell size
                %2. not containing all lower than median mCitr (if used)
                %3. not over max size
                %4. no drastic increase of more than 50% cell area
                if sum(sum(new_c_Image2==i))<min_cell_size ...
                        || sum(sum(uint8(newI).*Ignew(y_cn,x_cn)))./sum(sum(uint8(newI)))<0.99.*mean(mean(Ignew)) ...
                        || sum(sum(new_c_Image2>0))>max_allowed_cell_size ...
                        || area_increase>0.5
                    cell_exists(i,1)=0;
                    cell_exists(i,2)=c_time;
                end
                
                %store final segmented image in structure "all.obj"
                all_obj.cell_area(i,c_time)=sum(sum(new_c_Image2>0));
                all_obj.cells(:,:,c_time)=all_obj.cells(:,:,c_time)+new_c_Image2;
                all_obj.uniquecells(:,:,i, c_time) = new_c_Image2;
                Lcells=new_c_Image2; %make final labeled image
             
                %-------this code is for image display only-----------------------------------
                new_c_Image3=I<-1000;
                j5=i;
                if cell_exists(j5,1)==1
                    new_c_Image3=new_c_Image3+bwmorph(new_c_Image2==j5,'remove',inf);
                end
                
                Inew_here=uint8(zeros(x_size,y_size,3));
                Inew_here(:,:,1)=0.5.*I+     255.*uint8(bwmorph(new_c_Image3,'remove',inf));
                Inew_here(:,:,2)=0.5.*I+ 0.5.*IG+   0.*uint8(bwmorph(new_c_Image3,'remove',inf));
                Inew_here(:,:,3)=0.5.*I+       0.*uint8(bwmorph(new_c_Image3,'remove',inf));
                figure(4);imshow(Inew_here);title(c_time);drawnow
                
                if (cell_exists(i,1))>0
                    im_here=new_c_Image2==i;
                    s_ftmp=regionprops(bwlabel(im_here),'Centroid');
                    center_here=round(s_ftmp(1).Centroid);
                    text(center_here(1,1)+20, center_here(1,2),[num2str(i)],'BackgroundColor',[.7 .9 .7]);
                end
                hold off
                %--------------------end image display code--------------------------------
                 
                % get  data
                [all_obj]=OA_get_fl_GFPandmKok(all_obj,c_time,prefix,suffix2,pos_num,suffix,numbM,type,Lcells,cell_exists,i);
                [vol_c,area_c,maj_appr_axis,min_appr_axis]=OA_get_area(new_c_Image2,c_time,i); %i=cell number                         
                all_obj.volume(i,c_time)    = vol_c;           
                all_obj.area(i,c_time)      =  area_c;    
                all_obj.proj_area(i,c_time) = sum(new_c_Image2(:)>0);    
                
                % display fluorescent data
                figure(5);
                
                plot(all_obj.nuc_far1G(i,:),'k');hold on
                %plot(all_obj.far1_appr_conc(i,:),'r');
                hold off
                drawnow
                
            end %end of cell exists loop
        end
    end %end of cell loop
    
    ttime=toc;
    disp(['total time: ' num2str(ttime./60) 'min'])
    name1=[prefix '_pos_no_'  num2str(pos_num) '_2'];
    save(name1, 'all_obj', 'prefix', 'suffix2', 'pos_num', 'suffix', 'numbM', 'type', 'LcellsO' ...
        ,'cell_exists','max_size_vs_largest_cell','max_area_increase_per_tp','higher_threshold','lower_threshold' ...
        ,'threshold_increase_factor','phase_subtraction_factor','no_obj');
    close all
    %---------------------------------------------------------
end






