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
load('OA_030214_29re253-re263_12nM_pos_no_20_re_exp')


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
        disp(strcat('saturation' ,' time:' , num2str(initial_tp) , ' cell:' , num2str(no_cell_of_interest), ' no of sat. pixels (IN):', num2str(nosatpixelN)))
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
clc; close all; clear all;
perdifftot = 0;
j = 0;
per_decreases = zeros(1,1);
abs_decreases = zeros(1,1);
for pos = [1 14:19 21] %
    filename = 'OA_030714_OA039L381-L392-Longpulse_pos_no_';
    load([filename int2str(pos) '_re_exp_analysis']);
    load([filename int2str(pos) '_re_exp'])


    for i = 1:size(data,1)
        far1_st = round(data(i,9));
        far1_end = round(data(i,10));
        cellno = interesting_cells(i);
        backgr_st = round(data(i,5));
        backgr_end = round(data(i,6));
        vec_bckgr = backgr_st:backgr_end;
        if isempty(vec_bckgr)
            vec_bckgr = backgr_end:backgr_st;
        end
        
        backgr_nucFar1 = mean(all_obj.nuc_Far1(cellno,vec_bckgr));
        hold on

        if far1_end < 100
            j = j +1; 
            p1 = 7; % p(1) % this must be calculated properly using autofluorescence.
            p2 = backgr_nucFar1;%47.27;
            zdim = all_obj.volume(cellno,:)./all_obj.area(cellno,:);
            % eliminate NaNs.
            init_zdim = mean(zdim(vec_bckgr(~(isnan(zdim(vec_bckgr))))));
            autofl_nucfar1 = p1.*(zdim-init_zdim) + p2;
            nucFar1 = all_obj.nuc_Far1(cellno,:) - autofl_nucfar1;
            avefar1cur = sumfar1/j;
            nucfar1add = nucFar1;
            nucfar1add(isnan(nucfar1add)) = 0;
            nucfar1add = isnan(nucFar1).*avefar1cur + ~isnan(nucFar1).*nucfar1add;
            sumfar1 = sumfar1 + nucfar1add;    
            plot(nucFar1)
            diff = nucFar1(far1_st) - nucFar1(far1_end);
            perdiff = diff / nucFar1(far1_st);
            perdifftot = perdifftot + diff / nucFar1(far1_st);
                       
            abs_decreases(j) = diff;
            per_decreases(j) = perdiff;
        end
    end
end
per_decreases_pos = per_decreases(per_decreases < 1);
per_decreases_pos = per_decreases_pos(per_decreases_pos > -1);
abs_decreases_pos = abs_decreases(per_decreases < 1);
abs_decreases_pos = abs_decreases_pos(per_decreases_pos > -1);

avefar1 = sumfar1./size(per_decreases_pos,2);
plot(avefar1, 'r','LineWidth',3)
title('high copy nuc Far1-L92P in response to a met pulse in arrested cells','FontSize',12)
xlabel('frame no (/6 minutes)','FontSize',12)
ylabel('corrected nuclear Far1 concentration (a.u.)', 'FontSize',12)



hold off

%%

load('OA_030214_29re253-re263_12nM_pos_no_1_re2_exp')
cellno = 1;

p1 = 7; % p(1) % this must be calculated properly using autofluorescence.
p2 = 47.27;
zdim = all_obj.volume(cellno,:)./all_obj.area(cellno,:);
autofl_nucfar1 = p1.*(zdim-init_zdim) + p2;
nucFar1 = all_obj.nuc_Far1(cellno,:) - autofl_nucfar1;
plot(all_obj.Far1_appr_conc(cellno,:))


%%
clear all; clc;
% y1 = [7.07 32.94]; 
% e1 = [0.48 2.88];
% h = barwitherr(e1,y1);
y2 = [24 41];
e2 = [6 1];
h2 = barwitherr(e2,y2);

set(gca,'XTickLabel',{'wt Far1','Far1-L92P'})
ylabel('% Far1 Degradation (a.u.)')
set(h2(1),'FaceColor','k');
%%
clc; close all; clear all;
pulse_end = 26;
daughters_wt = zeros(1,1);
mothers_wt = zeros(1,1);
j1 = 1;
k1 = 1;
for pos = [2:14]
    filename = 'OA_122513_29255and29262_120nMpulseto3nM_pos_no_';
    load([filename int2str(pos) '_re_exp_analysis']);
    load([filename int2str(pos) '_re_exp'])


    for i = 1:size(data,1)
        whi5_exit = data(i,8);
        arresttime = 6*(whi5_exit-pulse_end);
        if data(i,3) == 0
            mothers_wt(j1) = arresttime;
            j1 = j1 + 1;
        else
            daughters_wt(k1) = arresttime;
            k1 = k1 + 1;
        end
    end
end

daughters_L92P = zeros(1,1);
mothers_L92P = zeros(1,1);
j2 = 1;
k2 = 1;

for pos = [1 15:16 18:20 22:28]
    filename = 'OA_122513_29255and29262_120nMpulseto3nM_pos_no_';
    load([filename int2str(pos) '_re_exp_analysis']);
    load([filename int2str(pos) '_re_exp'])


    for i = 1:size(data,1)
        whi5_exit = data(i,8);
        arresttime = 6*(whi5_exit-pulse_end);
        if data(i,3) == 0
            mothers_L92P(j2) = arresttime;
            j2 = j2 + 1;
        else
            daughters_L92P(k2) = arresttime;
            k2 = k2 + 1;
        end
    end
end        

mean(daughters_wt(daughters_wt))
mean(daughters_L92P(daughters_L92P))

%%
% run through all data

% for each data point for far1; use abs_Far1 as total far1

% take median of last 20 points as zero

% subtract that median from the far1

% fit the curve
fitobj = fit(time',far1','exp1');
% get coefficients
coeffs = coeffvalues(fitobj);
% inverse of second coefficient
halflife = -1/a(2);
%%

hold all
perdevpertime = zeros(1,1);
j = 1;

medianperdev = zeros(1,1);
k = 1;
for pos = [2 3 5]

    filename = 'OA_030214_29re253-re263_12nM_pos_no_';
    load([filename int2str(pos) '_re2_exp_analysis']);
    load([filename int2str(pos) '_re2_exp'])
    posdevs = zeros(1,1);
    for i = 1:size(data,1)

        far1_st = round(data(i,9));
        far1_end = round(data(i,10));

        backgr_st = round(data(i,5));
        backgr_end = round(data(i,6));
        vec_bckgr = backgr_st:backgr_end;
        if isempty(vec_bckgr)
            vec_bckgr = backgr_end:backgr_st;
        end
        cellno = interesting_cells(i);
        backgrfar1 = median(all_obj.Far1_appr_conc(cellno, vec_bckgr));

        far1curve = all_obj.Far1_appr_conc(cellno, far1_st:far1_end) - backgrfar1;
        far1curve = far1curve(far1curve > 0);
        plot(far1curve)
        expect = far1curve(1)*length(far1curve);
        deviate = expect - trapz(far1curve,2);
        perdevpertime(j) = (deviate/length(far1curve))/far1curve(1);
        posdevs(i) = perdevpertime(j);
        j = j +1;
    end
    if length(posdevs) > 3
        medianperdev(k) = median(posdevs);
        k = k + 1;
    end
    
end

hold off
%% no pulse
mothers_L92P = zeros(1,1);

k1 = 1;

for pos = [12:20]
    filename = 'OA_011514_29255and29262_3nMcontrol_pos_no_';
    load([filename int2str(pos) '_re_exp_analysis_2']);
    load([filename int2str(pos) '_re_exp'])


    for i = 1:size(data,1)
        whi5_exit = data(i,8);
        whi5_entry = data(i,7);
        arresttime = 6*(whi5_exit-whi5_entry);
        if data(i,3) == 0
            mothers_L92P(k1) = arresttime;
            k1 = k1 + 1;
        end
    end
end        

mothers_L92P = mothers_L92P(mothers_L92P > 30);