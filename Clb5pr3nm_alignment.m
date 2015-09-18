close all;clc;clear all;
file = 'OA_060415_AD3017c_3nM6min_pos_no_';
xvals = zeros(1,1);
yvals = zeros(1,1);

fnorm = 0.35; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');
chosencellpos = zeros(1,1);
chosencellno = zeros(1,1);
ind = 1;
m = 0;
d = 0;

dalignedatfiring = cell(1);
malignedatfiring = cell(1);
tdalignedatfiring = cell(1);
tmalignedatfiring = cell(1);

for pos = [3 6 7 8 9 10 11]
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
    
    
    GFPtot = all_obj.tot_nucl_areaR.*all_obj.nuc_Far1 + all_obj.tot_cyt_areaR.*all_obj.cyt_Far1;
    GFPconc = GFPtot./volume;
    GFPconc(find( (isinf(GFPconc)+isnan(GFPconc))>0) )=0;
    
    for i = 1:size(data,1)
        startp = round(data(i,5));
        reentry = round(data(i,6));
        endp = round(data(i,7))+3;
        cellno = data(i,1);
        filteredsig = filtfilt(b,a,GFPconc(cellno,:));
        filteredGFP = filteredsig(startp:endp);
        hold on
        plot(startp:endp,filteredGFP)
        plot([reentry reentry],[min(filteredGFP) max(filteredGFP)],'--','color',[1 0.0 1])
        reentryind = reentry-startp+1;
        Clb5afterarrest = filteredGFP(reentryind:end);
        Clb5atpeak = max(Clb5afterarrest);
        
        sortedclb5 = sort(Clb5afterarrest,'ascend');
        backgr = median(sortedclb5(1:3));

        Clb5increase = Clb5atpeak - backgr;
        Clb5halfmaxindex = min(find(Clb5afterarrest>(Clb5increase*.50) + backgr));
        Clb5incindex = min(find([0 diff(Clb5afterarrest)>0.01].*(Clb5afterarrest>(Clb5increase*.10) + backgr)==1));
        plot([reentry+Clb5incindex-1 reentry+Clb5incindex-1],[min(filteredGFP) max(filteredGFP)],'--','color',[1 0.0 0.0])
       
        checkbutton=1;
        while checkbutton ~= 2;
            if data(i,3) == 0
                [x1,y1,button] = ginput(1);
                if button ==1
                    Clb5incindex = x1 - reentry + 1;
                    Clb5halfmax_duration = (Clb5halfmaxindex - Clb5incindex)*6;
                    Clb5_Whi5_diff = Clb5incindex*6;
                    chosencellpos(ind) = pos;
                    chosencellno(ind) = cellno;
                    plot([reentry+Clb5incindex-1 reentry+Clb5incindex-1],[min(filteredGFP) max(filteredGFP)],'--','color',[1 0.0 0.0])
                    plot([reentry+Clb5halfmaxindex-1 reentry+Clb5halfmaxindex-1],[min(filteredGFP) max(filteredGFP)],'--','color',[1 0.0 0.0])
                    pause(1);
                    
                    checkbutton = 2;
                    ind = ind+1;
                    
                    m = m + 1;
                    peakx = find(filteredsig == Clb5atpeak);
                    alignend = peakx
                    malignedatfiring{m} = filteredsig(startp:alignend);
                    tmalignedatfiring{m} = ((startp:alignend) - round(x1)).*6;

                    
                    
                elseif button == 3
                    checkbutton = 2;
                end
                
            else
                [x1,y1,button] = ginput(1);
                
                if button ==1
                    chosencellpos(ind) = pos;
                    chosencellno(ind) = cellno;
                    Clb5incindex = x1 - reentry + 1;
                    Clb5halfmax_duration = (Clb5halfmaxindex - Clb5incindex)*6;
                    Clb5_Whi5_diff = Clb5incindex*6;
                    plot([reentry+Clb5incindex-1 reentry+Clb5incindex-1],[min(filteredGFP) max(filteredGFP)],'--','color',[1 0.0 0.0])
                    plot([reentry+Clb5halfmaxindex-1 reentry+Clb5halfmaxindex-1],[min(filteredGFP) max(filteredGFP)],'--','color',[1 0.0 0.0])
                    pause(1);
                    
                    checkbutton = 2;
                    ind = ind+1;
                    
                    d = d + 1;
                    peakx = find(filteredsig == Clb5atpeak);
                    alignend = peakx;
                    dalignedatfiring{d} = filteredsig(startp:alignend);
                    tdalignedatfiring{d} = ((startp:alignend) - round(x1)).*6;

                elseif button == 3
                    checkbutton = 2;
                end
            end
            
            
        end
        close all
    end
end

if d > 10
    save(['arrestClb53nm' '_alignment'], 'dalignedatfiring','tdalignedatfiring',...
        'malignedatfiring','tmalignedatfiring','chosencellpos','chosencellno')
end

%%
clear all;clc;close all
load 'arrestClb53nm_alignment'
hold all
td = zeros(1,1);
yd = zeros(1,1);
for i = [1:3 5:12]
    %dalignedatfiring{i} = dalignedatfiring{i}./max(dalignedatfiring{i});
    dalignedatfiring{i} = dalignedatfiring{i} - min(dalignedatfiring{i});
    plot(tdalignedatfiring{i},dalignedatfiring{i})
    curx = tdalignedatfiring{i}(1:end);
    cury = dalignedatfiring{i}(1:end);
    td = [td curx];
    yd = [yd cury];
end
td = td(2:end);
yd = yd(2:end);


[bin binmeds binstds] = makebins(td,yd,-150,100,20);
figure(2)
hold all
ciplot((binmeds-binstds), (binmeds+binstds), bin,'r')
plot(bin,binmeds,'LineWidth',3)
ylim([0 1])
xlim([-150 75])
hold off

    
    