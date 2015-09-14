cln3risetimes = [102.7930, 11.5866];
far1decaytimes = [6.4394,  1.1997];

barwitherr([11.5866 1.1997], [102.7930 6.4394])

%%
clear all; clc; close all;
backgroundfile = 'OA_070815_11Abar1d_refurb_3nM6min_pos_no_';
fnorm = 0.10; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');

position = 4;
cell = 15;
for pos = position
    load([backgroundfile int2str(pos) '_re_exp_vol']);
    load([backgroundfile int2str(pos) '_re_exp_vol_Cln3increase']);
    
    area = all_obj.tot_nucl_areaG + all_obj.tot_cyt_areaG;
    background1 = (89.29-area(area < 1440).*0.0023);%%%background subt.
    j2 = logical((area >= 1440).*(area <1765));
    background2 = (97.91-area(j2).*0.0086);
    background3 = (83.41+area(area >= 1765).*0.000);
    
    all_obj.nuc_whi5G(area < 1440) = all_obj.nuc_whi5G(area < 1440) - background1;
    all_obj.nuc_whi5G(j2) = all_obj.nuc_whi5G(j2)-background2; 
    all_obj.nuc_whi5G(area >= 1765) = all_obj.nuc_whi5G(area >= 1765)-background3; %%%%%%%%%%%%
    all_obj.cyt_whi5G(area < 1440) = all_obj.cyt_whi5G(area < 1440) - background1;
    all_obj.cyt_whi5G(j2) = all_obj.cyt_whi5G(j2)-background2; 
    all_obj.cyt_whi5G(area >= 1765) = all_obj.cyt_whi5G(area >= 1765)-background3;
    totCln3 = all_obj.tot_nucl_areaG.*all_obj.nuc_whi5G + all_obj.tot_cyt_areaG.*all_obj.cyt_whi5G;
    cln3conc = totCln3./all_obj.volume_ax;
    
    for i = 1:size(data,1)
        

        startp = data(i,7);
        endp = data(i,8);
        cellno = data(i,1);
        if cellno == cell;
            cell
            from = floor(startp)-5;
            to = ceil(endp)+5;

            cln3data = cln3conc(cellno,from:to);
            smoothcln3 = filtfilt(b,a,cln3data);
            hold all
            plot((from:to).*6, cln3data)
            plot((from:to).*6, smoothcln3)
            hold off

        end
    end
end
%%
backgroundfile = 'OA_061015_OA040_3nM6min_pos_no_';
fnorm = 0.10; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');
position = 13;
cell = 2;
for pos = position
    load([backgroundfile int2str(pos) '_re_exp_Whi5analysis']);
    load([backgroundfile int2str(pos) '_re_exp']);
    only = size(all_obj.volume,1);
    % background subtraction
    volume = all_obj.volume;
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
    
    
    % Whi5 signal
    volume(find( (isinf(volume)+isnan(volume))>0) )=0;
    
    Whi5nuc = all_obj.nuc_whi5R-all_obj.cyt_whi5R;
    Whi5tot = all_obj.tot_nucl_area_w5.*all_obj.nuc_whi5R + all_obj.tot_cyt_area_w5(1:only,:).*all_obj.cyt_whi5R;
    Whi5conc = Whi5tot./volume;
    Whi5conc(find( (isinf(Whi5conc)+isnan(Whi5conc))>0) )=0;
    
    
     for i = 1:size(data,1)
        
        startp = round(data(i,5));
        endp = round(data(i,6));
        cellno = data(i,1);
        % reentering cells
        if endp < 135
            if data(i,3) == 1 && cellno == cell %daughter cell
                hold all
                Whi5conccell = Whi5conc(cellno,:);
                Whi5nuccell = 0.035.*Whi5nuc(cellno,:);
                Whi5conccell = 1.3.*max(Whi5nuccell)/max(Whi5conccell).*Whi5conccell;
                smoothWhi5 = filtfilt(b,a,Whi5conccell);
                plot((1:133).*6,Whi5conccell(1:133),'LineWidth',2)
                plot((1:133).*6,Whi5nuccell(1:133),'LineWidth',2)
                plot((1:133).*6,smoothWhi5(1:133),'LineWidth',2)
            end
        end
     end
end
%%
clear all;clc;close all;
load OA_071615_OA049_50GalFar13nM_pos_no_decaytimes
decayhalflifed = log(0.5)*decaytimedaughters;
meddecay = median(decayhalflifed);
stddecay = std(bootstrp(100,@median,decayhalflifed));


load OA_061015_OA040_3nM6min_pos_no_arrestWhi5_OA040_at3nM
volumehalflifed= log(2)./ppdexpb;
medvol = median(volumehalflifed)
stdvol = std(bootstrp(1000,@median,volumehalflifed))

barwitherr([stddecay stdvol], [meddecay medvol])
ylim([0 150])

%%
clear all;clc;close all;

load arrestClb5pr3nM_results

indices = cumsum(darrestduration./6);
indices = [0 indices];
xre = zeros(1,1);
yre = dClb5duringarrest;
for i = 1:length(darrestduration)
    startp = indices(i);
    endp = indices(i+1);
    xre = [xre (0:endp-startp).*6];
end
xre = xre(2:end);
xre = xre(xre<400);
yre = yre(xre<400);


figure(1)
hold on
[bin binmeds binstds] = makebins(xre,yre,0,400,20);
binmeds = binmeds - min(binmeds);
ciplot(binmeds-binstds, binmeds+binstds, bin,'r')
plot(bin,binmeds,'LineWidth',3)
ylim([-0.4 0.4])
xlim([0 400])

%% Clb5 example trace
clc; close all; clear all;
file = 'OA_060415_AD3017c_3nM6min_pos_no_';
load arrestClb5pr3nM_results


fnorm = 0.95; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');


pos = 3; cellno = 23; i =4;
load([file int2str(pos) '_re_exp_analysis']);
load([file int2str(pos) '_re_exp']);

% background subtraction
all_obj.med_backgr_gfp = all_obj.med_backgr_Far1;
volume = all_obj.volume;

area= all_obj.tot_nucl_areaR+all_obj.tot_cyt_areaR; % check area % check 45 % correct backgr. subtraction
background1 = median(all_obj.med_backgr_gfp)-45+(74.3126+area(area < 3215).*0.0072);%%%background subt.
background3 = median(all_obj.med_backgr_gfp)-45+(100.0599-area(area >= 3215).*0.0011);

all_obj.nuc_Far1(area < 3215) = all_obj.nuc_Far1(area < 3215) - background1;
all_obj.nuc_Far1(area >= 3215) = all_obj.nuc_Far1(area >= 3215)-background3; %%%%%%%%%%%%
all_obj.cyt_Far1(area < 3215) = all_obj.cyt_Far1(area < 3215) - background1;
all_obj.cyt_Far1(area >= 3215) = all_obj.cyt_Far1(area >= 3215)-background3;

volume(find( (isinf(volume)+isnan(volume))>0) )=0;
startp = round(data(i,5));
reentry = round(data(i,6));
endp = round(data(i,7));

GFPtot = all_obj.tot_nucl_areaR.*all_obj.nuc_Far1 + all_obj.tot_cyt_areaR.*all_obj.cyt_Far1;
GFPconc = GFPtot./volume;
GFPconc(find( (isinf(GFPconc)+isnan(GFPconc))>0) )=0;
hold on
ylim([0 1.5])
filteredGFP = filtfilt(b,a,GFPconc(cellno,startp-1:endp+5));
plot((startp-1:endp+5).*6, filteredGFP)
%%
clear all;clc;close all
load arrestClb5pr3nM_results

medhalfmax = median(dClb5halfmaxduration);
stdhalfmax = std(bootstrp(1000,@median,dClb5halfmaxduration));
barwitherr([11.5866 stdhalfmax 1.1997], [102.7930 medhalfmax 6.4394])

%%

volume1_2 = [75.2655 3.1258];%**********
cln3risetimes = [102.7930, 11.5866];
clb5prtimes = [23.5282, 5.0682];
sic1d = [5.9185 1.2369];
far1decaytimes = [6.4394,  1.1997];


barwitherr([3.1258 11.5866 5.0682 1.2369 1.1997],[75.2655 102.7930 23.5282 5.9185 6.4394])
hold all
plot([75.2655 102.7930 23.5282 5.9185 6.4394],'o')

%%
cln3 = 0:0.01:1;
Km = 0.5;
n = 1;
b = 1;
hold all
for n = [1 2 4 8 16]
    clb5 = (b.*cln3.^n)./(Km^n+cln3.^n);
    plot(cln3,clb5,'LineWidth',1.5)
end
hold off
%%
devperapert =[0.5008    0.3692    0.2578    0.2344         
    0.4933    0.3516    0.2208    0.1831         
    0.4871    0.3382    0.1867    0.1323         
    0.4778    0.3537    0.2463    0.1422         
    0.4717    0.3812    0.2765    0.1452         
    0.4630    0.3953    0.2885    0.1424         
    0.4571    0.4012    0.2925    0.1379         
    0.4497    0.4028    0.2861    0.1284         
    0.4440    0.4067    0.2836    0.1215         
    0.4373    0.4086    0.2862    0.1148         
    0.4330    0.4085    0.2805    0.1077         
    0.4268    0.4064    0.2760    0.1008 ];

barwitherr(std(devperapert,1), mean(devperapert,1))
%% percentage decrease in nuclear far1
clc; close all; clear all;
perdifftot = 0;
j = 0;
per_decreases = zeros(1,1);
abs_decreases = zeros(1,1);
sumfar1 = zeros(1,1);

td = zeros(1,1);
xd = zeros(1,1);
yd = zeros(1,1);
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
            p1 = 7.49;
            p2 = backgr_nucFar1;%47.27;
            volume = all_obj.volume(cellno,:);
            zdim = all_obj.volume(cellno,:)./all_obj.area(cellno,:);
            % eliminate NaNs.
            init_zdim = mean(zdim(vec_bckgr(~(isnan(zdim(vec_bckgr))))));
            autofl_nucfar1 = p1.*(zdim-init_zdim) + p2;
            nucFar1 = all_obj.nuc_Far1(cellno,:) - autofl_nucfar1;
            avefar1cur = sumfar1/j;
            nucfar1add = nucFar1;
            nucFar1(isnan(nucFar1)) = 0;
            
            xd = [xd volume(31:56)];
            yd = [yd nucFar1(31:56)];
            td = [td 6.*(30:55)];
            
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
td = td(2:end);
yd = yd(2:end);
xd = xd(2:end);

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

figure(2)
hold on
[binl92p binmedsl92p binstdsl92p] = makebins(td,yd,30*6,55*6,30);
normf = 1; % max(binmedsl92p);
binmedsl92p = binmedsl92p./normf;
binstdsl92p = binstdsl92p./normf;
ciplot(binmedsl92p-binstdsl92p, binmedsl92p+binstdsl92p, binl92p,'g')
plot(binl92p,binmedsl92p,'b')
ylim([0 100])
xlim([30*6 55*6])
plot([240 240],[0 max(binmedsl92p)],'LineWidth',2)
plot([270 270],[0 min(binmedsl92p)],'LineWidth',2)

%%
clc; close all; clear all;
perdifftot = 0;
j = 0;
per_decreases = zeros(1,1);
abs_decreases = zeros(1,1);
sumfar1 = zeros(1,1);

td = zeros(1,1);
xd = zeros(1,1);
yd = zeros(1,1);
for pos = [2 3 4 5 8] %
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
            p1 = 7.49;
            p2 = backgr_nucFar1;%47.27;
            volume = all_obj.volume(cellno,:);
            zdim = all_obj.volume(cellno,:)./all_obj.area(cellno,:);
            % eliminate NaNs.
            init_zdim = mean(zdim(vec_bckgr(~(isnan(zdim(vec_bckgr))))));
            autofl_nucfar1 = p1.*(zdim-init_zdim) + p2;
            nucFar1 = all_obj.nuc_Far1(cellno,:) - autofl_nucfar1;
            avefar1cur = sumfar1/j;
            nucfar1add = nucFar1;
            nucFar1(isnan(nucFar1)) = 0;
            
            xd = [xd volume(31:56)];
            yd = [yd nucFar1(31:56)];
            td = [td 6.*(30:55)];
            
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
td = td(2:end);
yd = yd(2:end);
xd = xd(2:end);

per_decreases_pos = per_decreases(per_decreases < 1);
per_decreases_pos = per_decreases_pos(per_decreases_pos > -1);
abs_decreases_pos = abs_decreases(per_decreases < 1);
abs_decreases_pos = abs_decreases_pos(per_decreases_pos > -1);

avefar1 = sumfar1./size(per_decreases_pos,2);
plot(avefar1, 'r','LineWidth',3)
title('nuc Far1 in response to a met pulse in arrested cells','FontSize',12)
xlabel('frame no (/6 minutes)','FontSize',12)
ylabel('corrected nuclear Far1 concentration (a.u.)', 'FontSize',12)



hold off

figure(2)
hold on
[binwt binmedswt binstdswt] = makebins(td,yd,30*6,55*6,30);
normf = 1; % max(binmedswt);
binmedswt = binmedswt./normf;
binstdswt = binstdswt./normf;
ciplot(binmedswt-binstdswt, binmedswt+binstdswt, binwt,'r')
plot(binwt,binmedswt,'LineWidth',3)
ylim([0 100])
xlim([30*6 55*6])
plot([240 240],[0 max(binmedswt)],'LineWidth',2)
plot([270 270],[0 min(binmedswt)],'LineWidth',2)
%%
barwitherr([26.0 5.2],[42.7 -49.5])
hold all
plot([42.7 -49.5],'o')
ylim([-70 70])
%%
clear all;close all;clc;
backgroundfile = 'OA_071015_OA042_refurb_3nM6min_pos_no_';
for pos = 2
    load([backgroundfile int2str(pos) '_re_exp_Far1L92Pnucentry']);
    load([backgroundfile int2str(pos) '_re_exp']);
    cellno = data(i,1);
    if cellno == 13
        hold all
        
        

        curr_plot_Far1nuc = all_obj.nuc_Far1(cellno,:);
        p1 = 7.493; p2 = 47.27;
        zdim = all_obj.volume(cellno,:)./all_obj.area(cellno,:);
        autofl_nucfar1 = p1.*zdim + p2;
        far1nuc = curr_plot_Far1nuc-autofl_nucfar1;
        far1nuc(isnan(far1nuc)) = 0;
        whi5nuc = (all_obj.nuc_whi5R(cellno,:)-all_obj.cyt_whi5R(cellno,:))./zdim;

        far1nucnorm = far1nuc./zdim;
        plot(35*6:6:74*6,far1nucnorm(36:75))
        plot(35*6:6:74*6, whi5nuc(36:75),'r')
        plot([data(i,5)*6 data(i,5)*6], [0 max(far1nucnorm)])
        plot([data(i,6)*6 data(i,6)*6], [0 max(far1nucnorm)])
        xlim([35*6 75*6])
        ylim([0 30])
    end
end
%%
clear all;close all;clc;
backgroundfile = 'OA_061015_OA040_3nM6min_pos_no_';
pos = 6;

load([backgroundfile int2str(pos) '_re_exp_Far1analysis']);
load([backgroundfile int2str(pos) '_re_exp']);

i=3;
cellno = data(i,1);
if cellno == 17
    hold all
    
    
    
    curr_plot_Far1nuc = all_obj.nuc_Far1(cellno,:);
    p1 = 7.493; p2 = 47.27;
    zdim = all_obj.volume(cellno,:)./all_obj.area(cellno,:);
    autofl_nucfar1 = p1.*zdim + p2;
    far1nuc = curr_plot_Far1nuc-autofl_nucfar1;
    far1nuc(isnan(far1nuc)) = 0;
    whi5nuc = (all_obj.nuc_whi5R(cellno,:)-all_obj.cyt_whi5R(cellno,:))./zdim;
    
    far1nucnorm = far1nuc./zdim;
    plot(0*6:6:139*6,far1nucnorm(1:140))
    plot(0*6:6:139*6, whi5nuc(1:140),'r')
    plot([data(i,5)*6 data(i,5)*6], [0 max(far1nucnorm)])
    xlim([41*6 110*6])
    ylim([-1.1 4.9])
end

