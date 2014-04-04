close all


far1ratio = all_obj.abs_Far1./all_obj.volume;
totwhi5 = all_obj.cyt_whi5R.*all_obj.tot_cyt_areaR + all_obj.tot_nucl_areaR.* all_obj.nuc_whi5R;
w5ratio = totwhi5./all_obj.volume;
w5ratio(isinf(w5ratio)) = NaN;
mean(nanmean(w5ratio))
figure(1)
hold all
for i = 1:size(w5ratio,1)
    plot(w5ratio(i,:),'x')
end
hold off
figure(2)
far1ratio = all_obj.abs_Far1./all_obj.volume;
far1ratio(isinf(far1ratio)) = NaN;
mean(nanmean(far1ratio))

hold all
figure(2)
for i = 1:size(far1ratio,1)
    plot(far1ratio(i,:),'o')
end

hold off

%%
totwhi5 = all_obj.cyt_whi5R.*all_obj.tot_cyt_areaR + all_obj.tot_nucl_areaR.* all_obj.nuc_whi5R;
for i = 1:size(all_obj.Far1_appr_conc,1)
%     for j = 1:99
%         figure(1)        
%         image(all_obj.uniquecells(:,:,i,j));
%         figure(2)
%         image(all_obj.cells(:,:,j)) 
%         pause(0.005)
%     end
    figure(3) 
    hold all
    plot(all_obj.volume(i,:),all_obj.Far1_appr_conc(i,:))
    hold off
    hold all
    x = [0 1e4 2e4 3e4 4e4 5e4];
    y = 4.*x;
    plot(all_obj.volume(i,:),all_obj.abs_Far1(i,:))
    plot(x,y,'LineWidth',4)
    hold off
    figure(4)
    hold all
    plot(all_obj.volume(i,:),totwhi5(i,:))
    hold off
end
%% 
close all
load('OA_012514_OA037and38_12nMaFmetCln3pulses_pos_no_1_re_exp')
figure(1)
i = 2;             
hold all
plot(all_obj.nuc_whi5R(i,:)-all_obj.cyt_whi5R(i,:),'b')
nucfar1 = all_obj.cyt_Far1(i,:);
plot(nucfar1,'r')
hold off

figure(2)
plot(all_obj.Far1_appr_conc(i,:))
%%
no_tp   =numbM;%length(all_obj.max_nucl_int(1,:));
cell_data=-ones(1,49);
no_cells=length(all_obj.max_nucl_int(:,1));
%cell_data=0;%time of whi5-exit if cell is ok
disp(['cell no: ' num2str(no_cell_of_interest) '/' num2str(no_cells) ', left click for mother - right click for daughter at time of bud emergence'])

% intr_cell=0;
initial_tp=80;
%no_cell_of_interest=5;
c_num   ='c1';%phase
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);I=imread(im_name);
c_num   ='c3';%'c3';%here venus
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);IN=imread(im_name);
c_num   ='c2';%mcherry
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);IG=imread(im_name);

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
%%
% find saturated images

finaltp = length(all_obj.max_nucl_int(1,:));
threshold = 9;
no_cells=length(all_obj.max_nucl_int(:,1));
saturatedcells = zeros(0,0);
for initial_tp = 1:finaltp
    %no_cell_of_interest=5;
    c_num   ='c1';%phase
    [im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);I=imread(im_name);
    c_num   ='c3';%whi5-mkok
    [im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);IN=imread(im_name);
    c_num   ='c2';%far1-gfp
    [im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);IG=imread(im_name);
    
    for no_cell_of_interest = 1:no_cells
    
        nosatpixel = sum(sum((IG.*uint8(all_obj.cells(:,:,initial_tp)==no_cell_of_interest)) > 254));
        nosatpixelN = sum(sum((IN.*uint8(all_obj.cells(:,:,initial_tp)==no_cell_of_interest)) > 254));

        if ((nosatpixel > threshold) && ~ismember(no_cell_of_interest, saturatedcells))
        disp(strcat('saturation' ,' time:' , num2str(initial_tp) , ' cell:' , num2str(no_cell_of_interest), ' no of sat. pixels:', num2str(nosatpixel)))
        saturatedcells = [saturatedcells no_cell_of_interest];
        end

        if ((nosatpixelN > threshold) && ~ismember(no_cell_of_interest, saturatedcells))
        disp(strcat('saturation' ,' time:' , num2str(initial_tp) , ' cell:' , num2str(no_cell_of_interest), ' no of sat. pixels (IN):', num2str(nosatpixel)))
        saturatedcells = [saturatedcells no_cell_of_interest];
        end
    end
end
saturatedcells
    
%% background fluorescence
close all
figure(1)
hold on
for i = 1:size(data,1)
    whi5_entry = round(data(i,5));
    whi5_exit = round(data(i,6));
    cellno = interesting_cells(i);
    for j = whi5_entry:whi5_exit
        plot(all_obj.volume(cellno, j), all_obj.Far1_appr_conc(cellno,j),'x')
    end
end
hold off

figure(2)
hold on
a = 1;
for i = 1:size(data,1)
    whi5_entry = round(data(i,5));
    whi5_exit = round(data(i,6));
    cellno = interesting_cells(i);
    nucfar1bckg = zeros(1,1);
    for j = whi5_entry:whi5_exit
        
        volumes(a) = all_obj.volume(cellno, j);
        length(a) = all_obj.volume(cellno,j)/all_obj.area(cellno,j);
        nucfar1s(a) = all_obj.nuc_Far1(cellno,j);
        far1conc(a) = all_obj.Far1_appr_conc(cellno,j);
        plot(length(a), all_obj.nuc_Far1(cellno,j),'x')
        a = a + 1;
    end
end
hold off
perdiff = 0;
perdiffcyt = 0;
for i = 1:size(data,1)
    far1_st = round(data(i,7));
    far1_end = round(data(i,8));
    cellno = interesting_cells(i);
    diff = all_obj.nuc_Far1(cellno,far1_st) - all_obj.nuc_Far1(cellno,far1_end);
    perdiff = perdiff + diff / all_obj.nuc_Far1(cellno,far1_st);
    
    diffcyt = all_obj.cyt_Far1(cellno,far1_st) - all_obj.cyt_Far1(cellno,far1_end);
    perdiffcyt = perdiffcyt + diffcyt/ all_obj.cyt_Far1(cellno,far1_st);
end
perdiff/i
perdiffcyt/i
    

%% percentage decrease in nuclear far1
perdiff =0;
j = 0;
for i = 1:size(data,1)
    far1_st = round(data(i,9));
    far1_end = round(data(i,10));
    cellno = interesting_cells(i);
    backgr_st = round(data(i,5));
    backgr_end = round(data(i,6));
    backgr_nucFar1 = mean(all_obj.nuc_Far1(cellno,backgr_st:backgr_end));
    hold on
    if far1_end < 100
        p1 = 5; % p(1)
        p2 = backgr_nucFar1;%47.27;
        zdim = all_obj.volume(cellno,:)./all_obj.area(cellno,:);
        % eliminate NaNs.
        vec_bckgr = backgr_st:backgr_end;
        init_zdim = mean(zdim(vec_bckgr(~(isnan(zdim(vec_bckgr))))));
        autofl_nucfar1 = p1.*(zdim-init_zdim) + p2;
        nucFar1 = all_obj.nuc_Far1(cellno,:) - autofl_nucfar1;
        plot(nucFar1)
        diff = nucFar1(far1_st) - nucFar1(far1_end);
        diff / nucFar1(far1_st)
        perdiff = perdiff + diff / nucFar1(far1_st);
        j = j +1;
    end
end

perdiff/j
%%


