function [cell_data]=monitor_time_series_OA_Whi5(flc,initial_tp,no_cell_of_interest,...
    prefix,pos_num,suffix,type,suffix2,numbM,all_obj,lim_h)

%Whi5-mKOK 

no_tp   =numbM;%length(all_obj.max_nucl_int(1,:));
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
c_num   ='c2';%mcherry
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);IG=imread(im_name);

volume = all_obj.volume_fl;
area= all_obj.area;

background1 = median(all_obj.med_backgr_w5)/29.*(30.7743+area(area < 3215).*0.0010);%%%background subt.
j2 = logical((area >= 3215).*(area <4180));
background2 = median(all_obj.med_backgr_w5)/29.*(35.0318-area(j2).*0.0004);
background3 = median(all_obj.med_backgr_w5)/29.*(30.8450+area(area >= 4180).*0.0006);
% background_sp = all_obj.med_backgr_w5)/73.*(ppval(pp,area))

all_obj.nuc_whi5R(area < 3215) = all_obj.nuc_whi5R(area < 3215) - background1;
all_obj.nuc_whi5R(j2) = all_obj.nuc_whi5R(j2)-background2; %%%%%%%%%%%%
all_obj.nuc_whi5R(area >= 4180) = all_obj.nuc_whi5R(area >= 4180)-background3; %%%%%%%%%%%%
all_obj.cyt_whi5R(area < 3215) = all_obj.cyt_whi5R(area < 3215) - background1;
all_obj.cyt_whi5R(j2) = all_obj.cyt_whi5R(j2)-background2; %%%%%%%%%%%%
all_obj.cyt_whi5R(area >= 4180) = all_obj.cyt_whi5R(area >= 4180)-background3;

curr_plot_Far1T=[all_obj.nuc_Far1(no_cell_of_interest,:) zeros(1,5)];

only = size(all_obj.volume,1);
curr_plot_Whi5 = all_obj.tot_nucl_area_w5.*all_obj.nuc_whi5R + all_obj.tot_cyt_area_w5(1:only,:).*all_obj.cyt_whi5R;
curr_plot_Whi5 =[curr_plot_Whi5(no_cell_of_interest,:) zeros(1,5)];

curr_plot_Whi5N = all_obj.nuc_whi5R-all_obj.cyt_whi5R;
curr_plot_Whi5N = 10.*[curr_plot_Whi5N(no_cell_of_interest,:) zeros(1,5)];
%-----------------------
%curr_plot_Whi5C=1.*all_obj.nuc_Far1(no_cell_of_interest,:)-1.*all_obj.cyt_Far1(no_cell_of_interest,:);

 


curr_plot_volume = 0.01.*[volume(no_cell_of_interest,:) zeros(1,5)];
curr_plot_volume(find( (isinf(curr_plot_volume)+isnan(curr_plot_volume))>0) )=0;
curr_plot_Whi5C=[curr_plot_Whi5./curr_plot_volume zeros(1,5)];
curr_plot_Whi5C(find( (isinf(curr_plot_Whi5C)+isnan(curr_plot_Whi5C))>0) )=0;

curr_plot_Far1T = 25*all_obj.Far1_appr_conc(no_cell_of_interest,:);
curr_plot_Whi5= (5.*max(curr_plot_Far1T)./max(curr_plot_Whi5)).*curr_plot_Whi5;

%plot(curr_plot_area);pause
max_v=10+max(curr_plot_Far1T);

% figure(2)
% hold on
% plot(6:6:(84*6),3.*curr_plot_Whi5(1:84), 'color',[1.0 0.0 0.0],'LineWidth',3); % red
% plot(6:6:(84*6),3.*curr_plot_volume(1:84),'color',[0.7 0 0],'LineWidth',3); % dark red
% plot(6:6:(84*6),curr_plot_Whi5C(1:84),'color',[0.0 0.2 1.0],'LineWidth',3); % blue
% hold off
figure(1)
subplot(2,2,1)
plot(curr_plot_Far1T,'color',[0.0 0.5 0.0]); % green
hold on; 
plot(curr_plot_Whi5C,'color',[0.0 0.2 1.0]); % blue
plot(curr_plot_Whi5N,'color',[0.2 0.0 0.0]); % black
plot(curr_plot_Whi5, 'color',[1.0 0.0 0.0]); % red
plot(curr_plot_volume,'color',[0.7 0 0]); % dark red
title(['cell no = ' num2str(no_cell_of_interest)]);


for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end

% plot([all_obj.nuc_whi5R(no_cell_of_interest,:)-all_obj.cyt_whi5R(no_cell_of_interest,:) zeros(1,5)], 'color',[1.0 0.0 0.0]);hold on
% plot(curr_plot_Whi5C,'color',[0.0 0.2 1.0]);
 
 
 
for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 30],'k--','linewidth',1); end

subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,initial_tp)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,initial_tp),'remove')));
Nuc_outl=255.*uint8(bwmorph((all_obj.cells(:,:,initial_tp)==no_cell_of_interest).*all_obj.nuclLFar1(:,:,initial_tp),'remove',inf));
imshow(outline+Nuc_outl+I);title('phase')
text(10,10,['time = ' num2str(initial_tp)],'color',[0 1 0])
subplot(2,2,3)
imshow(outline+IN);title('Whi5-mkok')
text(10,10,['time = ' num2str(initial_tp)],'color',[0 1 0])
subplot(2,2,4)
imshow(outline+IG);title('Far1-GFP')
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
            subplot(2,2,1);plot(curr_plot_Far1T,'color',[0.0 0.5 0.0]);hold on; plot(curr_plot_Whi5C,'color',[0.0 0.2 1.0]);plot(curr_plot_Whi5N,'color',[0.2 0.0 0.0]);
            plot(curr_plot_Whi5, 'color',[1.0 0.0 0.0]);plot(curr_plot_volume,'color',[0.7 0 0]); title(['cell no = ' num2str(no_cell_of_interest)]);
            for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end
            plot(current_time,curr_plot_Far1T(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.5 0]);
            plot(current_time,curr_plot_Whi5C(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.2 1]);
            plot(current_time,curr_plot_Whi5N(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.2 0.0 0]);
            plot(current_time,curr_plot_volume(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.7 0 0]);
            plot(current_time,curr_plot_Whi5(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 0]);hold off
            c_num   ='c1';%phase
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            I=imread(im_name);subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,current_time)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,current_time),'remove')));
            Nuc_outl=255.*uint8(bwmorph((all_obj.cells(:,:,current_time)==no_cell_of_interest).*all_obj.nuclLFar1(:,:,current_time),'remove',inf));
            imshow(outline+I+Nuc_outl);title('phase')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   =flc;%'c3';%Far1
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IN=imread(im_name);subplot(2,2,3)
            imshow(outline+IN);title('Whi5-mKok')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   ='c2';%whi5
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IG=imread(im_name);subplot(2,2,4)
            imshow(outline+IG);title('Far1-GFP')
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
            subplot(2,2,1);plot(curr_plot_Far1T,'color',[0.0 0.5 0.0]);hold on; plot(curr_plot_Whi5C,'color',[0.0 0.2 1.0]);plot(curr_plot_Whi5N,'color',[0.2 0.0 0.0]);
            plot(curr_plot_Whi5, 'color',[1.0 0.0 0.0]);plot(curr_plot_volume,'color',[0.7 0 0]); title(['cell no = ' num2str(no_cell_of_interest)]);
            for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end
            plot(current_time,curr_plot_Far1T(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.5 0]);
            plot(current_time,curr_plot_Whi5C(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.2 1]);
            plot(current_time,curr_plot_Whi5N(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.2 0.0 0]);
            plot(current_time,curr_plot_volume(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.7 0 0]);
            plot(current_time,curr_plot_Whi5(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 0]);hold off
            
            c_num   ='c1';%phase
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            I=imread(im_name);subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,current_time)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,current_time),'remove')));
            Nuc_outl=255.*uint8(bwmorph((all_obj.cells(:,:,current_time)==no_cell_of_interest).*all_obj.nuclLFar1(:,:,current_time),'remove',inf));
            imshow(outline+I+Nuc_outl);title('phase')
            
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   =flc;%'c3';%Far1
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IN=imread(im_name);subplot(2,2,3)
            imshow(outline+IN);title('Whi5-mKok')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   ='c2';%whi5
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IG=imread(im_name);subplot(2,2,4)
            imshow(outline+IG);title('Far1-GFP')
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
            subplot(2,2,1);plot(curr_plot_Far1T,'color',[0.0 0.5 0.0]);hold on; plot(curr_plot_Whi5C,'color',[0.0 0.2 1.0]);plot(curr_plot_Whi5N,'color',[0.2 0.0 0.0]);
            plot(curr_plot_Whi5, 'color',[1.0 0.0 0.0]);plot(curr_plot_volume,'color',[0.7 0 0]); title(['cell no = ' num2str(no_cell_of_interest)]);
            for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end
            plot(current_time,curr_plot_Far1T(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.5 0]);
            plot(current_time,curr_plot_Whi5C(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.2 1]);
            plot(current_time,curr_plot_Whi5N(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.2 0.0 0]);
            plot(current_time,curr_plot_volume(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.7 0 0]);
            plot(current_time,curr_plot_Whi5(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 0]);hold off
            
            c_num   ='c1';%phase
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            I=imread(im_name);subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,current_time)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,current_time),'remove')));
            Nuc_outl=255.*uint8(bwmorph((all_obj.cells(:,:,current_time)==no_cell_of_interest).*all_obj.nuclLFar1(:,:,current_time),'remove',inf));
            imshow(outline+I+Nuc_outl);title('phase')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   =flc;%'c3';%Far1
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IN=imread(im_name);subplot(2,2,3)
            imshow(outline+IN);title('Whi5-mKok')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   ='c2';%whi5
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IG=imread(im_name);subplot(2,2,4)
            imshow(outline+IG);title('Far1-GFP')
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
            subplot(2,2,1);plot(curr_plot_Far1T,'color',[0.0 0.5 0.0]);hold on; plot(curr_plot_Whi5C,'color',[0.0 0.2 1.0]);plot(curr_plot_Whi5N,'color',[0.2 0.0 0.0]);
            plot(curr_plot_Whi5, 'color',[1.0 0.0 0.0]);plot(curr_plot_volume,'color',[0.7 0 0]); title(['cell no = ' num2str(no_cell_of_interest)]);
            for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end
            plot(current_time,curr_plot_Far1T(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.5 0]);
            plot(current_time,curr_plot_Whi5C(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.2 1]);
            plot(current_time,curr_plot_Whi5N(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.2 0.0 0]);
            plot(current_time,curr_plot_volume(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.7 0 0]);
            plot(current_time,curr_plot_Whi5(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 0]);hold off
            
            c_num   ='c1';%phase
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            I=imread(im_name);subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,current_time)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,current_time),'remove')));
            Nuc_outl=255.*uint8(bwmorph((all_obj.cells(:,:,current_time)==no_cell_of_interest).*all_obj.nuclLFar1(:,:,current_time),'remove',inf));
            imshow(outline+I+Nuc_outl);title('phase')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   =flc;%'c3';%Far1
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IN=imread(im_name);subplot(2,2,3)
            imshow(outline+IN);title('Whi5-mKok')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   ='c2';%whi5
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IG=imread(im_name);subplot(2,2,4)
            imshow(outline+IG);title('Far1-GFP')
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
        plot(curr_plot_Far1T,'color',[0.0 0.5 0.0]);hold on; 
        plot(curr_plot_Whi5C,'color',[0.0 0.2 1.0]);
        plot(curr_plot_Whi5N,'color',[0.2 0.0 0.0]); 
        plot(curr_plot_Whi5, 'color',[1.0 0.0 0.0]);
        plot(curr_plot_volume,'color',[0.7 0 0]); 
        title(['cell no = ' num2str(no_cell_of_interest)]);
        for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end

        if j==1
            xlabel('mark arrest (begin)')
            [x1,y1,btmp]=ginput(1);
            plot([x1 x1],[0 max_v],'color',[1 0 0],'linewidth',2)
            cell_data(1,5)=min([numbM x1]);
        elseif j==2
            xlabel('mark arrest (end)')
            [x1,y1,btmp]=ginput(1);
            plot([x1 x1],[0 max_v],'color',[1 0 0],'linewidth',2)
            cell_data(1,6)=min([numbM x1]);
%         elseif j==3
%             xlabel('mark (begin)-G1 for daughters')
%             [x1,y1,btmp]=ginput(1);
%             plot([x1 x1],[0 max_v],'color',[1 1 0],'linewidth',2)
%             cell_data(1,7)=min([numbM x1]);
%         elseif j==4
%             xlabel('mark (end)-G1 for daughters')
%             [x1,y1,btmp]=ginput(1);
%             plot([x1 x1],[0 max_v],'color',[1 1 0],'linewidth',2)
%             cell_data(1,8)=min([numbM x1]);
        end
    end
    
end