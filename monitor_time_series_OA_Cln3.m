function [cell_data]=monitor_time_series_OA_Cln3(flc,initial_tp,no_cell_of_interest,...
    prefix,pos_num,suffix,type,suffix2,numbM,all_obj,lim_h)



%middle button exits program
%left and right button as well as left and right arrow moves time +-1
%"1" and "3" on the numpad will move it 10 timeopints -+

no_tp   =numbM;%length(all_3obj.max_nucl_int(1,:));
cell_data=-ones(1,10);
no_cells=length(all_obj.max_nucl_int(:,1));
%cell_data=0;%time of whi5-exit if cell is ok
disp(['cell no: ' num2str(no_cell_of_interest) '/' num2str(no_cells) ', left click for mother - right click for daughter at time of bud emergence'])

% intr_cell=0;
%initial_tp=10;
%no_cell_of_interest=5;
c_num   ='c1';%phase
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);I=imread(im_name);
c_num   =flc;%'c3';%here venus
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);IN=imread(im_name);
IG=IN;

area = all_obj.tot_nucl_areaG + all_obj.tot_cyt_areaG;
background1 = (89.29-area(area < 1440).*0.0023);%%%background subt.
j2 = logical((area >= 1440).*(area <1765));
background2 = (97.91-area(j2).*0.0086);
background3 = (83.41+area(area >= 1765).*0.000);

all_obj.nuc_whi5G(area < 1440) = all_obj.nuc_whi5G(area < 1440) - background1;
all_obj.nuc_whi5G(j2) = all_obj.nuc_whi5G(j2)-background2; %%%%%%%%%%%%
all_obj.nuc_whi5G(area >= 1765) = all_obj.nuc_whi5G(area >= 1765)-background3; %%%%%%%%%%%%
all_obj.cyt_whi5G(area < 1440) = all_obj.cyt_whi5G(area < 1440) - background1;
all_obj.cyt_whi5G(j2) = all_obj.cyt_whi5G(j2)-background2; %%%%%%%%%%%%
all_obj.cyt_whi5G(area >= 1765) = all_obj.cyt_whi5G(area >= 1765)-background3; 
% all_obj.nuc_whi5G = all_obj.nuc_whi5G-background_sp;
% all_obj.cyt_whi5G = all_obj.cyt_whi5G-background_sp;

volume = all_obj.volume_ax;

curr_plot_tot_nucl_area=[all_obj.tot_nucl_areaG(no_cell_of_interest,:) zeros(1,5)];
curr_plot_area=[area(no_cell_of_interest,:) zeros(1,5)];
curr_plot_area(find( (isinf(curr_plot_area)+isnan(curr_plot_area))>0) )=0;
curr_plot_volume=[volume(no_cell_of_interest,:) zeros(1,5)];
curr_plot_volume(find( (isinf(curr_plot_volume)+isnan(curr_plot_volume))>0) )=0;

%[curr_plot_Cln3N]=fix_N_trace(curr_plot_Cln3N);
totCln3 = all_obj.tot_nucl_areaG.*all_obj.nuc_whi5G + all_obj.tot_cyt_areaG.*all_obj.cyt_whi5G;
curr_plot_Cln3T =[totCln3(no_cell_of_interest,:) zeros(1,5)];
curr_plot_Cln3C = curr_plot_Cln3T./(curr_plot_volume);
curr_plot_Cln3C(find( (isinf(curr_plot_Cln3C)+isnan(curr_plot_Cln3C))>0) )=0;
curr_plot_Cln3N =[all_obj.nuc_whi5G(no_cell_of_interest,:)-all_obj.cyt_whi5G(no_cell_of_interest,:) zeros(1,5)];
curr_plot_Cln3NC = curr_plot_Cln3N./(curr_plot_volume);
%-----------------------
%curr_plot_Cln3C=1.*all_obj.nuc_Far1(no_cell_of_interest,:)-1.*all_obj.cyt_Far1(no_cell_of_interest,:);


scale_factor1=max(curr_plot_Cln3T)./max(curr_plot_Cln3N);
curr_plot_Cln3N = 1.2*scale_factor1.*curr_plot_Cln3N;
scale_factor2=max(curr_plot_Cln3T)./max(curr_plot_Cln3C);
curr_plot_Cln3C = 0.7*scale_factor2.*curr_plot_Cln3C;
scale_factor3=max(curr_plot_Cln3T)./max((curr_plot_volume));
curr_plot_volume=1.0.*scale_factor3.*curr_plot_volume;
scale_factor4=max(curr_plot_Cln3T)./max(curr_plot_Cln3NC);
curr_plot_Cln3NC = 0.3*scale_factor4.*curr_plot_Cln3NC;

%p1 = 7.493; p2 = 47.27;
%zdim = all_obj.volume(no_cell_of_interest,:)./all_obj.area(no_cell_of_interest,:);
%autofl_nucfar1 = p1.*zdim + p2;
%plot(curr_plot_area);pause
max_v=10+max(curr_plot_Cln3T);

figure(1)
subplot(2,2,1)
plot(curr_plot_Cln3T,'color',[0.0 0.5 0.0]);hold on; % green
plot(curr_plot_Cln3C,'color',[0.0 0.2 1.0]); % blue
plot(curr_plot_Cln3N,'color',[0.2 0.0 0.0]); % black
plot(curr_plot_Cln3NC, 'color',[1.0 0.0 0.0]); % red
plot(curr_plot_volume,'color',[0.733 0 0]); % dark red
title(['cell no = ' num2str(no_cell_of_interest)]);
for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end

% plot([all_obj.nuc_whi5R(no_cell_of_interest,:)-all_obj.cyt_whi5R(no_cell_of_interest,:) zeros(1,5)], 'color',[1.0 0.0 0.0]);hold on
% plot(curr_plot_Cln3C,'color',[0.0 0.2 1.0]);
 

 
for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 30],'k--','linewidth',1); end

subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,initial_tp)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,initial_tp),'remove')));
Nuc_outl=255.*uint8(bwmorph((all_obj.cells(:,:,initial_tp)==no_cell_of_interest).*all_obj.nuclLR(:,:,initial_tp),'remove',inf));
imshow(outline+Nuc_outl+I);title('phase')
text(10,10,['time = ' num2str(initial_tp)],'color',[0 1 0])
subplot(2,2,3)
imshow(IN+Nuc_outl);title('Cln3-mCitrine')
text(10,10,['time = ' num2str(initial_tp)],'color',[0 1 0])




checkbutton=1;
current_time=initial_tp;
while checkbutton~=2
    [x1,y1,button] = ginput(1);
    nosatpixel = sum(sum((IG.*uint8(all_obj.cells(:,:,current_time)==no_cell_of_interest)) > 254));
    nosatpixelN = sum(sum((IN.*uint8(all_obj.cells(:,:,current_time)==no_cell_of_interest)) > 254));
    if (nosatpixel > 9)
    disp(strcat('saturation' ,' time:' , num2str(current_time) , ' cell:' , num2str(no_cell_of_interest), ' no of sat. pixels:', num2str(nosatpixel)))
    end
    
    if (nosatpixelN > 9)
    disp(strcat('saturation' ,' time:' , num2str(current_time) , ' cell:' , num2str(no_cell_of_interest), ' no of sat. pixels (IN):', num2str(nosatpixel)))
    end
    
    if (button==29) %forward in time
        current_time= current_time+1;
        if current_time>no_tp
            disp('end of time series')
            current_time= current_time-1;
        else
            figure(1)
            subplot(2,2,1);plot(curr_plot_Cln3T,'color',[0.0 0.5 0.0]);
            hold on; plot(curr_plot_Cln3C,'color',[0.0 0.2 1.0]);
            plot(curr_plot_Cln3N,'color',[0.2 0.0 0.0]);
            plot(curr_plot_Cln3NC, 'color',[1.0 0.0 0.0]);
            plot(curr_plot_volume,'color',[0.7 0 0]); 
            title(['cell no = ' num2str(no_cell_of_interest)]);
            for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end
            plot(current_time,curr_plot_Cln3T(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.5 0]);
            plot(current_time,curr_plot_Cln3C(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.2 1]);
            plot(current_time,curr_plot_Cln3N(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.2 0.0 0]);
            plot(current_time,curr_plot_volume(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.7 0 0]);
            plot(current_time,curr_plot_Cln3NC(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 0]);hold off
            
            c_num   ='c1';%phase
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            I=imread(im_name);subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,current_time)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,current_time),'remove')));
            Nuc_outl=255.*uint8(bwmorph((all_obj.cells(:,:,current_time)==no_cell_of_interest).*all_obj.nuclLR(:,:,current_time),'remove',inf));
            imshow(outline+I+Nuc_outl);title('phase')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   =flc;%'c3';%Far1
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IN=imread(im_name);subplot(2,2,3)
            imshow(IN+Nuc_outl);title('Cln3-mCitrine')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])

        end
        %---------------------first +-10 step------------------------------
    elseif button==51 %forward in time
        current_time= current_time+10;
        if current_time>no_tp
            current_time=no_tp;
            disp('end of time series')
        else
            figure(1)
            subplot(2,2,1);plot(curr_plot_Cln3T,'color',[0.0 0.5 0.0]);hold on; plot(curr_plot_Cln3C,'color',[0.0 0.2 1.0]);plot(curr_plot_Cln3N,'color',[0.2 0.0 0.0]);
            plot(curr_plot_Cln3NC, 'color',[1.0 0.0 0.0]);plot(curr_plot_volume,'color',[0.7 0 0]); title(['cell no = ' num2str(no_cell_of_interest)]);
            for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end
            plot(current_time,curr_plot_Cln3T(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.5 0]);
            plot(current_time,curr_plot_Cln3C(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.2 1]);
            plot(current_time,curr_plot_Cln3N(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.2 0.0 0]);
            plot(current_time,curr_plot_volume(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.7 0 0]);
            plot(current_time,curr_plot_Cln3NC(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 0]);hold off
            
            c_num   ='c1';%phase
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            I=imread(im_name);subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,current_time)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,current_time),'remove')));
            Nuc_outl=255.*uint8(bwmorph((all_obj.cells(:,:,current_time)==no_cell_of_interest).*all_obj.nuclLR(:,:,current_time),'remove',inf));
            imshow(outline+I+Nuc_outl);title('phase')
            
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   =flc;%'c3';%Far1
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IN=imread(im_name);subplot(2,2,3)
            imshow(Nuc_outl+IN+outline);title('Cln3-mCitrine')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])

            
            
        end
        %---------------------end first +-10 step--------------------------
    elseif  button==30
        checkbutton=2;
    elseif (button==28)
        current_time= current_time-1;
        if current_time<1
            disp('beginning of time series reached')
        else
            figure(1)
            subplot(2,2,1);plot(curr_plot_Cln3T,'color',[0.0 0.5 0.0]);hold on; plot(curr_plot_Cln3C,'color',[0.0 0.2 1.0]);plot(curr_plot_Cln3N,'color',[0.2 0.0 0.0]);
            plot(curr_plot_Cln3NC, 'color',[1.0 0.0 0.0]);plot(curr_plot_volume,'color',[0.7 0 0]); title(['cell no = ' num2str(no_cell_of_interest)]);
            for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end
            plot(current_time,curr_plot_Cln3T(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.5 0]);
            plot(current_time,curr_plot_Cln3C(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.2 1]);
            plot(current_time,curr_plot_Cln3N(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.2 0.0 0]);
            plot(current_time,curr_plot_volume(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.7 0 0]);
            plot(current_time,curr_plot_Cln3NC(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 0]);hold off
            
            c_num   ='c1';%phase
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            I=imread(im_name);subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,current_time)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,current_time),'remove')));
            Nuc_outl=255.*uint8(bwmorph((all_obj.cells(:,:,current_time)==no_cell_of_interest).*all_obj.nuclLR(:,:,current_time),'remove',inf));
            imshow(outline+I+Nuc_outl);title('phase')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   =flc;%'c3';%Far1
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IN=imread(im_name);subplot(2,2,3)
            imshow(Nuc_outl+IN+outline);title('Cln3-mCitrine')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])

            
        end
        %---------------------for the second +-10step---------------------
    elseif button==49
        current_time= current_time-10;
        if current_time<1
            current_time=1;
            disp('beginning of time series reached')
        else
            figure(1)
            subplot(2,2,1);plot(curr_plot_Cln3T,'color',[0.0 0.5 0.0]);hold on; plot(curr_plot_Cln3C,'color',[0.0 0.2 1.0]);plot(curr_plot_Cln3N,'color',[0.2 0.0 0.0]);
            plot(curr_plot_Cln3NC, 'color',[1.0 0.0 0.0]);plot(curr_plot_volume,'color',[0.7 0 0]); title(['cell no = ' num2str(no_cell_of_interest)]);
            for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end
            plot(current_time,curr_plot_Cln3T(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.5 0]);
            plot(current_time,curr_plot_Cln3C(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.2 1]);
            plot(current_time,curr_plot_Cln3N(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.2 0.0 0]);
            plot(current_time,curr_plot_volume(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.7 0 0]);
            plot(current_time,curr_plot_Cln3NC(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 0]);hold off
            
            c_num   ='c1';%phase
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            I=imread(im_name);subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,current_time)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,current_time),'remove')));
            Nuc_outl=255.*uint8(bwmorph((all_obj.cells(:,:,current_time)==no_cell_of_interest).*all_obj.nuclLR(:,:,current_time),'remove',inf));
            imshow(outline+I+Nuc_outl);title('phase')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   =flc;%'c3';%Far1
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IN=imread(im_name);subplot(2,2,3)
            imshow(Nuc_outl+IN+outline);title('Cln3-mCitrine')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])

            
        end
        %---------------end -+10 step-------------------------------------
    elseif button==1 %= left buttom for mother
        %    tmpX=min([numbM round(x1)]);tmpX2=min([numbM x1]);
        %  cell_data=[0 sum(sum((all_obj.cells(:,:,tmpX)==no_cell_of_interest)>0)) sum(sum((all_obj.cells(:,:,lim_h(1))==no_cell_of_interest)>0)) tmpX2 no_cell_of_interest pos_num];
        %0 for mother 1 for daughter
        cell_data(1:4)=[no_cell_of_interest pos_num 0 x1];
        checkbutton=2;
    elseif button==3 %= right button for daughter
        %   tmpX=min([numbM round(x1)]);tmpX2=min([numbM x1]);
        %  cell_data=[1 sum(sum((all_obj.cells(:,:,tmpX)==no_cell_of_interest)>0)) sum(sum((all_obj.cells(:,:,lim_h(1))==no_cell_of_interest)>0)) tmpX2 no_cell_of_interest pos_num];
        %0 for mother 1 for daughter
        cell_data(1:4)=[no_cell_of_interest pos_num 1 x1];
        checkbutton=2;
    else
        disp('odd error')
    end
end

close all

% %interesting cell
if button ==1 || button ==3
    for j=1:2
        plot(curr_plot_Cln3T,'color',[0.0 0.5 0.0]);
        hold on; 
        plot(curr_plot_Cln3C,'color',[0.0 0.2 1.0]);
        plot(curr_plot_Cln3N,'color',[0.2 0.0 0.0]); 
        plot(curr_plot_Cln3NC, 'color',[1.0 0.0 0.0]);
        plot(curr_plot_volume,'color',[0.7 0 0]); 
        title(['cell no = ' num2str(no_cell_of_interest)]);
        for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end

%         if j==1
%             xlabel('pre-arrest-begin')
%             [x1,y1,btmp]=ginput(1);
%             plot([x1 x1],[0 max_v],'color',[1 0 0],'linewidth',2)
%             cell_data(1,5)=min([numbM x1]);
%         elseif j==2
%             xlabel('pre-arrest-end')
%             [x1,y1,btmp]=ginput(1);
%             plot([x1 x1],[0 max_v],'color',[1 0 0],'linewidth',2)
%             cell_data(1,6)=min([numbM x1]);
        if j==1
            xlabel('arrest-begin')
            [x1,y1,btmp]=ginput(1);
            plot([x1 x1],[0 max_v],'color',[1 1 0],'linewidth',2)
            cell_data(1,7)=min([numbM x1]);
        elseif j==2
            xlabel('arrest-end')
            [x1,y1,btmp]=ginput(1);
            plot([x1 x1],[0 max_v],'color',[1 1 0],'linewidth',2)
            cell_data(1,8)=min([numbM x1]);
         end
    end
    
end










