clc; close all; clear all;
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
dClb5duringarrest = zeros(1,1);
mClb5duringarrest = zeros(1,1);
dClb5afterarrest = zeros(1,1);
mClb5afterarrest = zeros(1,1);
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
        startp = round(data(i,5))
        reentry = round(data(i,6))
        endp = round(data(i,7))+3
        cellno = data(i,1);
        filteredGFP = filtfilt(b,a,GFPconc(cellno,startp:endp));

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
        data(i,3)
        pos
        cellno
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
                    mClb5duringarrest = [mClb5duringarrest GFPconc(cellno,startp:reentry)];
                    marrestduration(m) = (reentry-startp)*6;
                    mClb5halfmaxduration(m) = Clb5halfmax_duration;
                    mClb5Whi5diff(m) = Clb5_Whi5_diff;
                    
                    mClb5afterarrest = [mClb5afterarrest Clb5afterarrest];
                    mClb5increase(m) = Clb5increase;
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
                    dClb5duringarrest = [dClb5duringarrest GFPconc(cellno,startp:reentry)];
                    dClb5afterarrest = [dClb5afterarrest Clb5afterarrest];
                    dClb5increase(d) = Clb5increase;
                    darrestduration(d) = (reentry-startp)*6;
                    dClb5halfmaxduration(d) = Clb5halfmax_duration;
                    dClb5Whi5diff(d) = Clb5_Whi5_diff;
                    if cellno == 23 
                        pause(50) 
                    end
                elseif button == 3
                    checkbutton = 2;
                end
            end
            
            
        end
        close all
    end
end
dClb5duringarrest = dClb5duringarrest(2:end);
dClb5afterarrest = dClb5afterarrest(2:end);

if ind > 15
    save(['arrestClb5pr3nM' '_' 'results'],'darrestduration','dClb5duringarrest','dClb5halfmaxduration',...
        'dClb5Whi5diff','dClb5increase','dClb5afterarrest','marrestduration','mClb5duringarrest','mClb5halfmaxduration',...
        'mClb5Whi5diff','mClb5increase','mClb5afterarrest','chosencellpos','chosencellno');
end
