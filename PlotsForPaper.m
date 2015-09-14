
% cln3 oscillations--concentrations (background subtracted values)
clear all;
load arrestCln3_oscillations
G1high = mean(cln3oscillationdata(:,3));
beforearrest = mean(cln3oscillationdata(:,4));
arresthigh = mean(cln3oscillationdata(:,5));
cellno = length(cln3oscillationdata(:,3));
G1highse = std(cln3oscillationdata(:,3))/sqrt(cellno);
beforearrestse = std(cln3oscillationdata(:,4))/sqrt(cellno);
arresthighse = std(cln3oscillationdata(:,5))/sqrt(cellno);
yvals = [G1high beforearrest arresthigh];
evals = [G1highse beforearrestse arresthighse];
ylim([0 0.4])
errorbar(1:3,yvals,evals,'x','LineWidth',3);
% n = 24

%% Cln3 cell for example trace
clear all;clc;
filename = 'OA_041315_OA064_11Abar1del_6nM_pos_no_5_re_exp_vol';
load(filename)
cellno = 3; 
%pos=4;cellno=18;%116;%63; %55; example mother
%pos=5;cellno=3;180;70; example daughter
initial_tp=180; %180;70
c_num   ='c1';%phase
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);
I=imread(im_name);
c_num   ='c4';
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);
IN=imread(im_name);

outline=(255.*uint8(bwmorph(all_obj.cells(:,:,initial_tp)==cellno,'remove')));
figure(1)
imshow(I+outline);title('phase')
figure(2)
imshow(outline+IN);title('mCitrine')


%% Cln3 hill coefficients
clear all; clc;
load arrestCln3_hills

prevpos = chosencellpos(1);
file = 'OA_041315_OA064_11Abar1del_6nM_pos_no_';
load([file int2str(prevpos) '_re_exp_vol_analysis']);
m = 1;
d = 1;
for i = 1:length(hills);
    pos = chosencellpos(i);
    if pos ~= prevpos % load correct pos
        load([file int2str(pos) '_re_exp_vol_analysis']);
        prevpos = pos;
    end
    cellno = chosencellno(i); % which cell in that pos;
    for j = 1:size(data,1) % scan through all cells in that pos
        if cellno == data(j,1) % found the right cell
            if data(j,3) == 0 % mother cell
                hillsm(m) = hills(i);
                m = m +1;
            else
                hillsd(d) = hills(i);
                d = d + 1;
            end
        end
    end
end
median(hillsm) %2.0794; n = 29
median(hillsd) %2.1827; n = 32
median(hills) %2.0794; n = 61
hist(hills,[1,3,5,7,9,11])
med = @(x)median(x);    
hillboot = bootstrp(1000,med,hills);
hillbootm = bootstrp(1000,med,hillsm);
hillbootd = bootstrp(1000,med,hillsd);
stdofhill = std(hillboot); %0.3225
stdofhillm = std(hillbootm); %0.3837
stdofhilld = std(hillbootd); %0.5753

%% Clb5 hill coefficients
%assuming ~20 min GFP maturation time
clear all;close all;clc
load arrestClb5_hills
median(hills) % 24.7068 % n = 32
median(hillsd) % 28.8400 % n =16
median(hillsm) % 19.6759 % n =16
med = @(x)median(x);    
hillboot = bootstrp(1000,med,hills); 
hillbootm = bootstrp(1000,med,hillsm);
hillbootd = bootstrp(1000,med,hillsd); 
%922 out of 1000 bigger
% p > 0.05 for difference (p = 0.08)
stdofhill = std(hillboot);  % 3.31
stdofhillm = std(hillbootm); % 4.26
stdofhilld = std(hillbootd); % 4.65

hist(hills(hills<50),[1 6 11 16 21 26 31 36 41 46 51])

%% Clb5 example cell and trace
clear all;clc;
filename = 'AD_wtc5p_30_17c_oct31_2013_pos_no_5_vol';
load(filename)
cellno = 1; 
%pos=5; cellno =136;159
initial_tp=159; % h = 16.0
c_num   ='c1';%phase
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);
I=imread(im_name);
c_num   ='c4';
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);
IN=imread(im_name);

outline=(255.*uint8(bwmorph(all_obj.cells(:,:,initial_tp)==cellno,'remove')));
figure(1)
imshow(I+outline);title('phase')
figure(2)
imshow(outline+IN);title('GFP')
%% whi5 plots are in arrestWhi5 and arrestWhi5at3nM m files
% here, only volume is plotted
clc;clear all;close all;
load AD_wtc5p_30_17c_oct31_2013_pos_no_arrestWhi5at3nM

daugterdouble3nM = median(log(2)./ppdexpb); 
motherdouble3nM = median(log(2)./ppmexpb); 
med = @(x)median(x);    
bootd3nM = bootstrp(1000,med,log(2)./ppdexpb); 
bootm3nM = bootstrp(1000,med,log(2)./ppmexpb);

stdofdaug3nM = std(bootd3nM); 
stdofmoth3nM = std(bootm3nM);

b = median(ppdexpb);
a = median(ppdexpa);
bootd3nMa = bootstrp(1000,med,ppdexpa);
bootd3nMb = bootstrp(1000,med,ppdexpb);
stda = std(bootd3nMa);
stdb = std(bootd3nMb);

xvals = 0:1:200;
yvals = a.*exp(b.*xvals);
highyvals = (a+stda).*exp((b+stdb).*xvals);
upper_er = highyvals - yvals;
lowyvals = (a-stda).*exp((b-stdb).*xvals);
lower_er = yvals - lowyvals;
figure(1)
errorbar(xvals,yvals,lower_er,upper_er)
xlim([0 200])

load OA_060514_OA040and42_4-5nM_pos_no_arrestWhi5_OA040_at4_5nM

daugterdouble45nM = median(log(2)./ppdexpb); 
motherdouble45nM = median(log(2)./ppmexpb); 
med = @(x)median(x);    
bootd45nM = bootstrp(1000,med,log(2)./ppdexpb); 
bootm45nM = bootstrp(1000,med,log(2)./ppmexpb);

stdofdaug45nM = std(bootd45nM); 
stdofmoth45nM = std(bootm45nM);
figure(2)
errorbar(1:2,[daugterdouble3nM daugterdouble45nM],[stdofdaug3nM stdofdaug45nM],'x','LineWidth',3);
ylim([0 180])
figure(3)
errorbar(1:2,[motherdouble3nM motherdouble45nM],[stdofmoth3nM stdofmoth45nM],'x','LineWidth',3);
ylim([0 180])
%% example cell traces for whi5

clear all;clc;
filename = 'OA_060514_OA040and42_4-5nM_pos_no_2_re_exp';
load(filename)
cellno = 1; 
%pos=2; cellno =1; t 25;50;75;
initial_tp=75; % t=25
c_num   ='c1';%phase
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);
I=imread(im_name);
c_num   ='c3';
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);
IN=imread(im_name);

outline=(255.*uint8(bwmorph(all_obj.cells(:,:,initial_tp)==cellno,'remove')));
figure(1)
imshow(I+outline);title('phase')
figure(2)
imshow(outline+IN);title('mKok')

clear all;clc;close all;
filename = 'AD_wtc5p_30_17c_oct31_2013_pos_no_7_vol';
load(filename)
cellno = 17; 
initial_tp=118; %118;70
c_num   ='c1';%phase
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);
I=imread(im_name);
c_num   ='c4';
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);
IN=imread(im_name);

outline=(255.*uint8(bwmorph(all_obj.cells(:,:,initial_tp)==cellno,'remove')));
figure(3)
imshow(I+outline);title('phase')
figure(4)
imshow(outline+IN);title('GFP')


volume = all_obj.volume;
area= all_obj.area;
median(all_obj.med_backgr_gfp)
background1 = median(all_obj.med_backgr_gfp)/45.*(74.3126+area(area < 3215).*0.0072);%%%background subt.
background3 = median(all_obj.med_backgr_gfp)/45.*(100.0599-area(area >= 3215).*0.0011);
only = size(all_obj.volume,1);

all_obj.nuc_whi5R(area < 3215) = all_obj.nuc_whi5R(area < 3215) - background1;
all_obj.nuc_whi5R(area >= 3215) = all_obj.nuc_whi5R(area >= 3215)-background3; %%%%%%%%%%%%
all_obj.cyt_whi5R(area < 3215) = all_obj.cyt_whi5R(area < 3215) - background1;
all_obj.cyt_whi5R(area >= 3215) = all_obj.cyt_whi5R(area >= 3215)-background3;

curr_plot_Whi5 = all_obj.tot_nucl_areaR.*all_obj.nuc_whi5R + all_obj.tot_cyt_areaR(1:only,:).*all_obj.cyt_whi5R;
curr_plot_Whi5 =0.001.*[curr_plot_Whi5(cellno,:) zeros(1,5)];
curr_plot_volume = 0.01.*[volume(cellno,:) zeros(1,5)];
curr_plot_volume(find( (isinf(curr_plot_volume)+isnan(curr_plot_volume))>0) )=0;
curr_plot_Whi5C=1000.*[curr_plot_Whi5./curr_plot_volume zeros(1,5)];
curr_plot_Whi5C(find( (isinf(curr_plot_Whi5C)+isnan(curr_plot_Whi5C) + curr_plot_Whi5C>2000)>0) )=0;


figure(5)
hold on
plot(3:3:(120*3),12.*curr_plot_Whi5(1:120), 'color',[1.0 0.0 0.0],'LineWidth',3); % red
plot(3:3:(120*3),3.*curr_plot_volume(1:120),'color',[0.7 0 0],'LineWidth',3); % dark red
plot(3:3:(120*3),curr_plot_Whi5C(1:120),'color',[0.0 0.2 1.0],'LineWidth',3); % blue
hold off



