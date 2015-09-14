function [cell_data]=monitor_time_series_s1w5_oa1(flc,initial_tp,no_cell_of_interest,...
    prefix,pos_num,suffix,type,suffix2,numbM,all_obj,lim_h,no_cells)


%for the Far1 venus + whi5 mcherry strain

%middle button exits program
%left and right button as well as left and right arrow moves time +-1
%"1" and "3" on the numpad will move it 10 timeopints -+


no_tp   =numbM;%length(all_obj.max_nucl_int(1,:));
cell_data=-ones(1,23);
%cell_data=0;%time of whi5-exit if cell is ok
disp(['cell no: ' num2str(no_cell_of_interest) '/' num2str(no_cells) ', left click for mother - right click for daughter at time of cc-arrest'])


% intr_cell=0;
%initial_tp=10;
%no_cell_of_interest=5;
c_num   ='c1';%phase
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);I=imread(im_name);
c_num   =flc;%'c3';%here venus
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);IN=imread(im_name);
c_num   ='c2';%mcherry
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);IG=imread(im_name);


curr_plot_Far1 =[all_obj.nuc_sic1G(no_cell_of_interest,:)-all_obj.cyt_sic1G(no_cell_of_interest,:) zeros(1,5)];

curr_plot_Whi5 =[all_obj.nuc_whi5R(no_cell_of_interest,:)-all_obj.cyt_whi5R(no_cell_of_interest,:) zeros(1,5)];

curr_plot_Far1(curr_plot_Far1>(3.*max(curr_plot_Whi5)))=3.*max(curr_plot_Whi5);

curr_plot_area =zeros(1,length(curr_plot_Far1));
curr_plot_Far1C=zeros(1,length(curr_plot_Far1));
curr_plot_Far1T=zeros(1,length(curr_plot_Far1));
curr_plot_Far1N=zeros(1,length(curr_plot_Far1));

% curr_plot_Far1C=[all_obj.Far1_appr_conc(no_cell_of_interest,:) zeros(1,5)];
% curr_plot_area=[(all_obj.abs_Far1(no_cell_of_interest,:)./all_obj.Far1_appr_conc(no_cell_of_interest,:)).^(2./3) zeros(1,5)];
% curr_plot_area(find( (isinf(curr_plot_area)+isnan(curr_plot_area))>0) )=0;
% scale_factor1= max(curr_plot_Far1)./max(curr_plot_Far1C);
% curr_plot_Far1C=0.9.*scale_factor1.*curr_plot_Far1C;
% scale_factor2=(max(curr_plot_Far1C)./max(all_obj.abs_Far1(no_cell_of_interest,:)));
% curr_plot_Far1T=scale_factor2.*[all_obj.abs_Far1(no_cell_of_interest,:) zeros(1,5)];%abs Far1 scaled
% scale_factor3=max(curr_plot_Far1)./max(curr_plot_area);
% curr_plot_area=0.15.*scale_factor3.*curr_plot_area;
% 
% curr_plot_area2=[(all_obj.abs_Far1(no_cell_of_interest,:)./all_obj.Far1_appr_conc(no_cell_of_interest,:)).^(2./3) zeros(1,5)];
% curr_plot_area2(find( (isinf(curr_plot_area2)+isnan(curr_plot_area2))>0) )=0;
% 
% %plot(curr_plot_area);pause
% curr_plot_Far1C(1,41)=0.5.*(curr_plot_Far1C(1,40)+curr_plot_Far1C(1,42));
% curr_plot_Far1T(1,41)=0.5.*(curr_plot_Far1T(1,40)+curr_plot_Far1T(1,42));


figure(1)
subplot(2,2,1)
plot(curr_plot_Far1,'b');hold on;plot(curr_plot_area.^(3./2),'color',[0.3 0 0]);plot(curr_plot_Whi5,'r');plot(curr_plot_Far1C,'m');plot(curr_plot_Far1T,'color',[0 0.5 0]);plot(curr_plot_Far1N,'k');title(['cell no = ' num2str(no_cell_of_interest)])

max_v=10+max(curr_plot_Far1N);
for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2);end
subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,initial_tp)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,initial_tp),'remove')));
imshow(outline+I);title('phase')
text(10,10,['time = ' num2str(initial_tp)],'color',[0 1 0])
subplot(2,2,3)
imshow(IN);title('Far1-venus')
text(10,10,['time = ' num2str(initial_tp)],'color',[0 1 0])
subplot(2,2,4)
imshow(outline+IG);title('Whi5-mCherry')
text(10,10,['time = ' num2str(initial_tp)],'color',[0 1 0])

checkbutton=1;
current_time=initial_tp;
while checkbutton~=2
    [x1,y1,button] = ginput(1);
    if (button==29) %forward in time
        current_time= current_time+1;
        if current_time>no_tp
            disp('end of time series')
        else
            figure(1)
            subplot(2,2,1);plot(curr_plot_Far1,'b');hold on;plot(curr_plot_area.^(3./2),'color',[0.3 0 0]);plot(curr_plot_Whi5,'r');plot(curr_plot_Far1C,'m');plot(curr_plot_Far1T,'color',[0 0.5 0]);plot(curr_plot_Far1N,'k');title(['cell no = ' num2str(no_cell_of_interest)])
            for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end
            plot(current_time,curr_plot_Far1(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0 1]);
            plot(current_time,curr_plot_Far1N(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0 0]);
            
            plot(current_time,curr_plot_Far1T(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.5 0]);
            plot(current_time,curr_plot_Far1C(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 1]);
            plot(current_time,curr_plot_area(current_time).^(3./2),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.3 0 0]);
            
            plot(current_time,curr_plot_Whi5(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 0]);hold off
            c_num   ='c1';%phase
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            I=imread(im_name);subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,current_time)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,current_time),'remove')));
            imshow(outline+I);title('phase')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   =flc;%'c3';%Far1
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IN=imread(im_name);subplot(2,2,3)
            imshow(IN);title('Far1-venus')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   ='c2';%whi5
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IG=imread(im_name);subplot(2,2,4)
            imshow(outline+IG);title('Whi5-mCherry')
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
            subplot(2,2,1);plot(curr_plot_Far1,'b');hold on;plot(curr_plot_area.^(3./2),'color',[0.3 0 0]);plot(curr_plot_Far1N,'k');plot(curr_plot_Far1C,'m');plot(curr_plot_Far1T,'color',[0 0.5 0]);plot(curr_plot_Whi5,'r');title(['cell no = ' num2str(no_cell_of_interest)])
            for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end
            plot(current_time,curr_plot_Far1(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0 1]);
            plot(current_time,curr_plot_Far1N(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0 0]);
            
            plot(current_time,curr_plot_Far1T(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.5 0]);
            plot(current_time,curr_plot_Far1C(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 1]);
            plot(current_time,curr_plot_area(current_time).^(3./2),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.3 0 0]);
            
            plot(current_time,curr_plot_Whi5(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 0]);hold off
            c_num   ='c1';%phase
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            I=imread(im_name);subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,current_time)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,current_time),'remove')));
            imshow(outline+I);title('phase')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   =flc;%'c3';%Far1
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IN=imread(im_name);subplot(2,2,3)
            imshow(IN);title('Far1-venus')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   ='c2';%whi5
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IG=imread(im_name);subplot(2,2,4)
            imshow(outline+IG);title('Whi5-mCherry')
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
            subplot(2,2,1);plot(curr_plot_Far1,'b');hold on;plot(curr_plot_area.^(3./2),'color',[0.3 0 0]);plot(curr_plot_Far1N,'k');plot(curr_plot_Far1C,'m');plot(curr_plot_Far1T,'color',[0 0.5 0]);plot(curr_plot_Whi5,'r');title(['cell no = ' num2str(no_cell_of_interest)])
            for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end
            plot(current_time,curr_plot_Far1(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0 1]);
            plot(current_time,curr_plot_Far1N(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0 0]);
            
            plot(current_time,curr_plot_Far1T(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.5 0]);
            plot(current_time,curr_plot_Far1C(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 1]);
            plot(current_time,curr_plot_area(current_time).^(3./2),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.3 0 0]);
            
            plot(current_time,curr_plot_Whi5(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 0]);hold off
            c_num   ='c1';%phase
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            I=imread(im_name);subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,current_time)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,current_time),'remove')));
            imshow(outline+I);title('phase')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   =flc;%'c3';%Far1
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IN=imread(im_name);subplot(2,2,3)
            imshow(IN);title('Far1-venus')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   ='c2';%whi5
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IG=imread(im_name);subplot(2,2,4)
            imshow(outline+IG);title('Whi5-mCherry')
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
            subplot(2,2,1);plot(curr_plot_Far1,'b');hold on;plot(curr_plot_area.^(3./2),'color',[0.3 0 0]);plot(curr_plot_Far1N,'k');plot(curr_plot_Far1C,'m');plot(curr_plot_Far1T,'color',[0 0.5 0]);plot(curr_plot_Whi5,'r');title(['cell no = ' num2str(no_cell_of_interest)])
            for ii=1:length(lim_h);plot([lim_h(ii) lim_h(ii)],[0 max_v],'k','linewidth',2); end
            plot(current_time,curr_plot_Far1(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0 1]);
            plot(current_time,curr_plot_Far1N(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0 0]);
            
            plot(current_time,curr_plot_Far1T(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0 0.5 0]);
            plot(current_time,curr_plot_Far1C(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 1]);
            plot(current_time,curr_plot_area(current_time).^(3./2),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[0.3 0 0]);
            
            plot(current_time,curr_plot_Whi5(current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',[1 0 0]);hold off
            c_num   ='c1';%phase
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            I=imread(im_name);subplot(2,2,2);outline=(255.*uint8(bwmorph(all_obj.cells(:,:,current_time)==no_cell_of_interest,'remove'))+0.*uint8(bwmorph(all_obj.cells(:,:,current_time),'remove')));
            imshow(outline+I);title('phase')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   =flc;%'c3';%Far1
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IN=imread(im_name);subplot(2,2,3)
            imshow(IN);title('Far1-venus')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])
            c_num   ='c2';%whi5
            [im_name] = get_image_name(prefix,pos_num,suffix,current_time, c_num,type,suffix2,numbM);
            IG=imread(im_name);subplot(2,2,4)
            imshow(outline+IG);title('Whi5-mCherry')
            text(10,10,['time = ' num2str(current_time)],'color',[0 1 0])            
        end
        %---------------end -+10 step-------------------------------------
    elseif (button==1) || (button==3) %1 = left buttom for mother, 3 = right button for daughter      
        cell_data(1:8)=[no_cell_of_interest pos_num x1 0 0 0 0 (button-1)>0 ];
        checkbutton=2;
    else
        disp('odd error')
        button
    end  
end
close all

%-----------------------------------------------------------------------------------

% 1.cell no
% 2.fov no
% 3.arrest time (det by whi5)
% 4.whi5 exit time
% 5.sic1 exit time
% 6.dw5/dt at cell cycle reentry
% 7.ds1/dt at cell cycle reentry
% 8.cell type (0 = m) or (1 = d)
% 9.sic1 midexit from fitpoints   
%10.cell size at beginning of sic1 drop
%11.cell size at end of drop
%12.hill coefficient for fit
%13.'K' for fit
%14. rsquare for fit
%15. sse for fit
%16. time for full sic1 drop;
%17. start time (vs exp time) of sic1 drop
%18. end time (vs exp time) of sic1 drop
%19.hill coefficient for area fit
%20.'K' for area fit
%21. rsquare for area fit
%22. sse for area fit
%23. ...as 19 for smoothened area


% sic1 drop in reentering cells:
% ---example cell and its trace (as a function of time)
% ---delta time for sic1 drop (histogram for a population of reentering cells)
% ---cell size at the beginning and ending of sic1 drop (for many cells)-- Ultimately we would like to get a histogram of hill coefficients: For instance, 
%    if we have this time data, we can plot normalized sic1 (from 1 to 0) to normalized change in cell size for many cells for this time period, fit hill 
%    coefficients to these plots, and have a histogram of hill coefficients for sic1 drop (we should need to think about how best to do this to avoid some 
%    common problems that others make in use of hill coefficients).
% ---median time that the cell stays arrested until deciding to reenter (to point out that sic1 does not decrease for x hours while cln3 is increasing in this same period)

% K=2;
% 
% for K=1:10
% X=0:0.001:10;
% n=5;
% plot(X,1-1./(1+(K./X).^n));hold on
% pause
% end




% % interesting cell
 if button ==1 || button ==3

     sic1_window=15;
     plot(curr_plot_Far1,'b');hold on;plot(curr_plot_Whi5,'r');hold on
     
     % [pks,loc] = findpeaks(X);

     
% %     plot(curr_plot_area.^(3./2),'color',[0.3 0 0]);hold on;plot(curr_plot_Whi5,'r');
% %     title(['cell no = ' num2str(no_cell_of_interest)])
     xlabel('mark time of cell cycle reentry (Whi5 exit - red)')
     [x1,y1,btmp]=ginput(1);
     plot([x1,x1],[0 max(curr_plot_Far1)],'k--')
     cell_data(1,4)=x1;
     xlabel('mark time of Sic1 (blue line) nucl exit at cc reentry')
     [x1,y1,btmp]=ginput(1);
     cell_data(1,5)=x1;
     
     
     tmpWL=1;
     while tmpWL~=0
     
     close all
     
     plot(curr_plot_Far1((-sic1_window+round(cell_data(1,5)):(sic1_window+round(cell_data(1,5))))));hold on;
     xlabel('select start point for sic1 fit')
     [x1,y1,tmpb]=ginput(1);
     tmpS=[x1,y1];
     plot(x1,y1,'o');
     xlabel('select endpoint for sic1 fit')
     [x1,y1,tmpb]=ginput(1);
     tmpE=[x1,y1];
     
     tmpCh=round(cell_data(1,5))-sic1_window:round(cell_data(1,5))+sic1_window;
     
     
     cell_data(1,7)=(tmpS(1,2)-tmpE(1,2))./(tmpE(1,1)-tmpS(1,1));   
     cell_data(1,9)= cell_data(1,5) + tmpS(1,1)+0.5.*(tmpE(1,1)-tmpS(1,1)) - (sic1_window +1);
     
     
     SP=round(cell_data(1,5))+round(tmpS(1,1))-(1+sic1_window);%start point
     EP=round(cell_data(1,5))+round(tmpE(1,1))-(1+sic1_window);%end point
     
     
     
    % SP=round(cell_data(1,9)+round(tmpS(1,1))-(1+sic1_window));%start point
    % EP=round(cell_data(1,9)+round(tmpE(1,1))-(1+sic1_window));%end point
     
     
     
     CH=all_obj.cells(:,:,SP)==no_cell_of_interest;
     cell_data(1,10)=sum(CH(:)>0);
     CH=all_obj.cells(:,:,EP)==no_cell_of_interest;
     cell_data(1,11)=sum(CH(:)>0);
    
     to_fit=curr_plot_Far1(SP:EP);
     bfaftl=5;
     tf2= [ (to_fit-(min(to_fit)))./(max(to_fit-(min(to_fit)))) zeros(1,bfaftl)];
     %tf2= [ones(1,bfaftl) (to_fit-(min(to_fit)))./(max(to_fit-(min(to_fit)))) zeros(1,bfaftl)]; 
% --- hill coeff fit for normalized sic1 ---
     to_fit_conc=tf2;
     x_conc=1:length(to_fit_conc); 
     xi=linspace(x_conc(1),x_conc(end),100);Y=tf2; yi = interp1(Y,xi);to_fit_conc=yi;x_conc=xi;
     
     fcn_h2=fittype('1-(x^n/(K+x^n))');options = fitoptions(fcn_h2);options.StartPoint = [to_fit_conc(1,1) to_fit_conc(1,end)];  
     [fitobject2,stats2] = fit(x_conc',to_fit_conc',fcn_h2,options);
     figure(3);plot(x_conc,1 -( x_conc.^(fitobject2.n)./((fitobject2.K)+ x_conc.^(fitobject2.n))),x_conc,to_fit_conc);title('sic1 decay hill coeff fit') 
      
     disp(['n = ' num2str(fitobject2.n)])

cell_data(1,12)=fitobject2.n;
cell_data(1,13)=fitobject2.K;
cell_data(1,14)=stats2.rsquare;
cell_data(1,15)=stats2.sse;
cell_data(1,16)=tmpS(1,1)-tmpE(1,1);



%========================================================================================
% fit the area to sic1 hill coeff
cell_data(1,17)= tmpCh(round(tmpS(1,1)));
cell_data(1,18) =tmpCh(round(tmpE(1,1)));

%area_tmp=zeros(1,length(1,round(tmpE(1,1))-round(tmpS(1,1))));
area_tmp=1;

 %area_tmp=zeros(1,2.*ws+1);
      % area_tmp2=1;
      % sic1_tmp=zeros(1,2.*ws+1);
      % sic1_tmp=all_obj.nuc_sic1G(cell_no,round(sd_time)-ws:round(sd_time)+ws)-all_obj.cyt_sic1G(cell_no,round(sd_time)-ws:round(sd_time)+ws);
       
       index=1;
       for j=cell_data(1,17):(cell_data(1,18)+bfaftl)
           CH=all_obj.cells(:,:,j)==no_cell_of_interest;
           area_tmp(index)=sum(CH(:)>0);
           index=index+1;         
       end
       
       %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
           CH=all_obj.cells(:,:,round(cell_data(1,3)))==no_cell_of_interest;
           area_tmp2(1)=sum(CH(:)>0);
        
       
       
       %===================================================================================
       
     to_fit_conc=tf2;
     %x_conc=1:length(to_fit_conc); 
     %x_conc=area_tmp./area_tmp(1); 
     %x_conc=(1+area_tmp)-area_tmp(1);
     
     x_conc=(area_tmp./area_tmp2)-(area_tmp(1)./area_tmp2-1); 
     

     %xi=linspace(x_conc(1),x_conc(end),100);Y=tf2; yi = interp1(Y,xi);to_fit_conc=yi;x_conc=xi;
     
     fcn_h2=fittype('1-(x^n/(K+x^n))');options = fitoptions(fcn_h2);options.StartPoint = [to_fit_conc(1,1) to_fit_conc(1,end)];  
     [fitobject2,stats2] = fit(x_conc',to_fit_conc',fcn_h2,options);
     figure(4);plot(x_conc,1 -( x_conc.^(fitobject2.n)./((fitobject2.K)+ x_conc.^(fitobject2.n))),x_conc,to_fit_conc);title('sic1 decay hill coeff fit') 
      
     disp(['n = ' num2str(fitobject2.n)])

cell_data(1,19)=fitobject2.n;
cell_data(1,20)=fitobject2.K;
cell_data(1,21)=stats2.rsquare;
cell_data(1,22)=stats2.sse;

xN=linspace(x_conc(1),x_conc(end),100);
figure(4);hold on;plot(xN,1 -( xN.^(fitobject2.n)./((fitobject2.K)+ xN.^(fitobject2.n))),'r')

 
index=1;
area_tmp3=zeros(1,numbM);
       for j=1:numbM
           CH=all_obj.cells(:,:,j)==no_cell_of_interest;
           area_tmp3(index)=sum(CH(:)>0);
           index=index+1;         
       end
       area_tmp3=(area_tmp3./area_tmp2)-(area_tmp3(round(cell_data(3)))./area_tmp2-1); 
       
%figure(6);plot(area_tmp3);

cpS1=all_obj.nuc_sic1G(no_cell_of_interest,:)-all_obj.cyt_sic1G(no_cell_of_interest,:);
w5t= all_obj.nuc_whi5R(no_cell_of_interest,:)-all_obj.cyt_whi5R(no_cell_of_interest,:); 

figure(7);rh=124:153;plot(area_tmp3(rh),cpS1(rh),'-o')
figure(6);plot(area_tmp3(rh));
figure(8);plot(rh,w5t(rh),rh,cpS1(rh))
%figure(6);plot(area_tmp3(1:end),cpS1,'-o')


  %===================================================================================
     area_tmp=smooth(area_tmp)';  
     to_fit_conc=tf2;
   
     x_conc=(area_tmp./area_tmp2)-(area_tmp(1)./area_tmp2-1); 
    % xi=linspace(x_conc(1),x_conc(end),100);Y=tf2; yi = interp1(Y,xi);to_fit_conc=yi;x_conc=xi;
     
     fcn_h2=fittype('1-(x^n/(K+x^n))');options = fitoptions(fcn_h2);options.StartPoint = [to_fit_conc(1,1) to_fit_conc(1,end)];  
     [fitobject2,stats2] = fit(x_conc',to_fit_conc',fcn_h2,options);
     figure(5);plot(x_conc,1 -( x_conc.^(fitobject2.n)./((fitobject2.K)+ x_conc.^(fitobject2.n))),x_conc,to_fit_conc);title('sic1 decay hill coeff fit') 
      
     disp(['n = ' num2str(fitobject2.n)])

cell_data(1,23)=fitobject2.n;


xlabel('left button = ok, other button = redo')
[xx,yy,tmpB]=ginput(1)
if tmpB==1
    tmpWL=0;
%else;tmpWL=1;
end


%xlabel('paused - click to continue')
%pause
     end %while loop

 end % interesting cell


%