clear all
close all

add_data=0;
if add_data==1

%load('OA_040815_JS21211b-9nMrepeat_pos_no_11_1');lim_h=[42 numbM];
load('OA_042015_RV200_6nM_pos_no_5_1');lim_h=[42 numbM];

all_cell_data=-1.*ones(1,23); %[-1 -1 -1 -1 -1];

no_cells=length(all_obj.mean_cytGFP(:,1));
intr_cells_here=1:no_cells;
disp(['no of cells here = ' num2str(no_cells)])
index=1;

for no_cell_of_interest=1:no_cells
    %   if cell_exists(no_cell_of_interest,2)<32
    initial_tp=lim_h(end);
    intr_cell= monitor_time_series_s1w5_oa1('c4',initial_tp,no_cell_of_interest,prefix,pos_num,suffix,type,suffix2,numbM,all_obj,lim_h,no_cells);
    if intr_cell(1,1)~=-1
        all_cell_data(index,:)=intr_cell;
        disp(all_cell_data)
        index=index+1;
    end
    %   end
end

if all_cell_data(1,1)~=-1
  %  save([prefix '_' num2str(pos_num) '_s1w5_oa_2'], 'all_cell_data','lim_h','prefix','pos_num','suffix','type','suffix2','numbM');
else
    disp('nothing saved for this FOV')
end
close all


elseif add_data==0
 
     load('OA_040815_JS21211b-9nMrepeat_1_s1w5_oa_2');cell_data_OA=all_cell_data;tmp_p=[2 4 5 8:11];
     for tmppos=tmp_p;load(['OA_040815_JS21211b-9nMrepeat_' num2str(tmppos) '_s1w5_oa_2']); cell_data_OA=[cell_data_OA; all_cell_data];end

    tmp_p=[5 6 8 12 14 15]; for tmppos=tmp_p;load(['OA_042015_RV200_6nM_' num2str(tmppos) '_s1w5_oa_2']); cell_data_OA=[cell_data_OA; all_cell_data];end
     
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
%[h,p]=kstest2(cell_data_OA(1:23,12),cell_data_OA(24:end,12))

% disp('----- time -----')
%      [x1,y1]=AD_cdf1(cell_data_OA(:,12));figure(1);plot(x1,y1);title('hill coeff vs time')
%     [ mean(x1) sqrt(var(x1)) sqrt(var(x1))./sqrt(length(x1)) length(x1) ]
%   disp('----- area -----')
%   [x1,y1]=AD_cdf1(cell_data_OA(:,19));figure(2);plot(x1,y1);title('hill coeff vs area')
%     [ mean(x1) sqrt(var(x1)) sqrt(var(x1))./sqrt(length(x1)) length(x1) ]
%   disp('----- smoothened area -----')
%   [x1,y1]=AD_cdf1(cell_data_OA(:,23));figure(3);plot(x1,y1);title('hill coeff vs area')
%     [ mean(x1) sqrt(var(x1)) sqrt(var(x1))./sqrt(length(x1)) length(x1) ]
%     
    
 %     [x2,y2]=AD_cdf1(cell_data_FK(:,7));
%     plot(x1,y1,x2,y2)


% to add new hill coeff fit vs size
% add_new_hc=0;
% if add_new_hc==1
%     
%  
%     for i=1:31 % first exp
%        if i<24
%         load(['OA_040815_JS21211b-9nMrepeat_pos_no_' num2str(cell_data_OA(i,2)) '_1']);
%        else
%            load(['OA_042015_RV200_6nM_pos_no_' num2str(cell_data_OA(i,2)) '_1']);lim_h=[42 numbM];
%        end
%        cell_no=cell_data_OA(i,1);
%        st_size=cell_data_OA(i,10);
%        ed_size=cell_data_OA(i,11);
%        
%        sd_time=cell_data_OA(i,9);
%        ws=5;
%        area_tmp=zeros(1,2.*ws+1);
%        area_tmp2=zeros(1,numbM);
%       % sic1_tmp=zeros(1,2.*ws+1);
%        sic1_tmp=all_obj.nuc_sic1G(cell_no,round(sd_time)-ws:round(sd_time)+ws)-all_obj.cyt_sic1G(cell_no,round(sd_time)-ws:round(sd_time)+ws);
%        
%        index=1;
%        for j=round(sd_time)-ws:round(sd_time)+ws
%            CH=all_obj.cells(:,:,j)==cell_no;
%            area_tmp(index)=sum(CH(:)>0);
%            index=index+1;         
%        end
%        
%        %============================
%        index=1;
%        for j=1:numbM
%            CH=all_obj.cells(:,:,j)==cell_no;
%            area_tmp2(index)=sum(CH(:)>0);
%            index=index+1;         
%        end
%        %============================
%        area_tmp2=smooth(area_tmp2);
%        area_tmp=area_tmp2(round(sd_time)-ws:round(sd_time)+ws);
%        
%        
%        figure(1);plot(area_tmp,'-o');hold on;
%        plot([1 2.*ws+1],[st_size st_size],'k--')
%        plot([1 2.*ws+1],[ed_size ed_size],'k--')
% %      
% 
%      
% 
%       [q1,w1]=min(abs(area_tmp(1:ws)-st_size));
%       [q2,w2]=min(abs(area_tmp(ws+1:end)-ed_size));
%       
%       area_tmp(w1:ws+w2)
%       sic1_tmp(w1:ws+w2)
%       
%       area_tmpO=area_tmp(w1:ws+w2);     
%       area_tmp=area_tmp(w1:ws+w2)./area_tmp(w1);
%       sic1_tmp=(sic1_tmp(w1:ws+w2)-sic1_tmp(ws+w2))./(sic1_tmp(w1)-sic1_tmp(ws+w2));
%       
%       rc=round(sd_time)-ws:round(sd_time)+ws;
%       figure(3);plot(area_tmp2);hold on;plot(rc(w1:ws+w2),area_tmpO,'k--','linewidth',5);hold off
%      %figure(3);plot(area_tmp2);hold on;plot(round(sd_time)-ws:round(sd_time)+ws,area_tmp,'k--','linewidth',5);hold off
%       
%       figure(4);plot(area_tmp,sic1_tmp,'-o');title(i)   
%       
%       figure(5);plot(sic1_tmp,'r-o')
%       
%       drawnow 
%       
%       pause
% %       figure(2);plot(sic1_tmp);hold on;
% %     %  plot([1 2.*ws+1],[st_size st_size],'k--')
% %       plot([round(1.*ws+1) round(1.*ws+1)],[min(sic1_tmp) max(sic1_tmp)],'k--')
% %        
% %       
%       
%       
%     end
%   
% end
end