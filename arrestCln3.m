clear all; clc; close all;
backgroundfile = 'OA_041315_OA064_11Abar1del_6nM_pos_no_';
xvals = zeros(1,1);
yvals = zeros(1,1);

fnorm = 0.10; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');
chosencellpos = zeros(1,1);
chosencellno = zeros(1,1);
ind = 1;
hills = zeros(1,1);
for pos = [4:12]
    load([backgroundfile int2str(pos) '_re_exp_vol_analysis']);
    load([backgroundfile int2str(pos) '_re_exp_vol']);
    
    % background subtraction
    area = all_obj.tot_nucl_areaG + all_obj.tot_cyt_areaG;
    background1 = median(all_obj.med_backgr_gfp)/73.*(80.7082+area(area < 1135).*0.0030);%%%background subt.
    j2 = logical((area >= 1135).*(area <1820));
    background2 = median(all_obj.med_backgr_gfp)/73.*(88.2416-area(j2).*0.0036);
    background3 = median(all_obj.med_backgr_gfp)/73.*(77.6698+area(area >= 1820).*0.0020);
    all_obj.nuc_whi5G(area < 1135) = all_obj.nuc_whi5G(area < 1135) - background1;
    all_obj.nuc_whi5G(j2) = all_obj.nuc_whi5G(j2)-background2; %%%%%%%%%%%%
    all_obj.nuc_whi5G(area >= 1820) = all_obj.nuc_whi5G(area >= 1820)-background3; %%%%%%%%%%%%
    all_obj.cyt_whi5G(area < 1135) = all_obj.cyt_whi5G(area < 1135) - background1;
    all_obj.cyt_whi5G(j2) = all_obj.cyt_whi5G(j2)-background2; %%%%%%%%%%%%
    all_obj.cyt_whi5G(area >= 1820) = all_obj.cyt_whi5G(area >= 1820)-background3;
    
    
    
    % Cln3signal
    volume = all_obj.volume_fl;
    volume(find( (isinf(volume)+isnan(volume))>0) )=0;
    
    
    Cln3tot = all_obj.tot_nucl_areaG.*all_obj.nuc_whi5G + all_obj.tot_cyt_areaG.*all_obj.cyt_whi5G;
    Cln3conc = Cln3tot./volume;
    Cln3conc(find( (isinf(Cln3conc)+isnan(Cln3conc))>0) )=0;
    
    for i = 1:size(data,1)
        close all
        hold all
        
        startp = round(data(i,7));
        endp = round(data(i,8));
        cellno = data(i,1)
        xvals = [xvals volume(cellno, startp:endp)];
        yvals = [yvals Cln3conc(cellno,startp:endp)];
        %         maximum = max(volume(cellno,startp:endp));
        %         minimum = min(volume(cellno,startp:endp));
        %         binNo = 50;
        %
        %         bins= linspace(minimum, maximum, binNo);
        %
        %         [h,specBins] = histc(xvals, bins);
        %
        %         for j = 1:binNo
        %             binvalues    = yvals(specBins == j);
        %             binmeds(j)     = median(binvalues);
        %         end
        %
        
        volume_sm = volume(cellno,(startp-5):endp);%filtfilt(b,a,volume(cellno, startp:endp));
        volume_sm_norm = volume_sm./median(volume_sm(1:10));
        cln3conc_sm = Cln3conc(cellno,(startp-5):endp);%filtfilt(b,a,Cln3conc(cellno,startp:endp));
        sortedcln3 = sort(cln3conc_sm,'ascend');
        cln3conc_sm = cln3conc_sm - median(sortedcln3(1:10));
        sortedcln3 = sort(cln3conc_sm,'descend');
        cln3conc_sm_norm = cln3conc_sm./median(sortedcln3(1:10));
        plot(filtfilt(b,a,volume_sm_norm),filtfilt(b,a,cln3conc_sm_norm))
        volume_sm_norm(volume_sm_norm <= 1) = NaN;
        cln3conc_sm_norm(cln3conc_sm_norm <= 0) = NaN;
        plot(volume_sm_norm,cln3conc_sm_norm,'x')
        
        fo_ = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',[0 0  0],'Upper',[15   2 200]);
        ok_ = isfinite(volume_sm_norm) & isfinite(cln3conc_sm_norm);
        st_ = [0.037051126859541816 0.13223694988779544 0.2171532878989586];
        set(fo_,'Startpoint',st_);
        ft_ = fittype('a.*((vol-1).^n)./(Km.^n + (vol-1).^n)',...
            'dependent',{'y'},'independent',{'vol'},...
            'coefficients',{'Km', 'a', 'n'});
        try
            cf_ = fit(volume_sm_norm(ok_)',cln3conc_sm_norm(ok_)',ft_,fo_);
        catch err
            display(err)
        end
        
        h_ = plot(cf_,'fit',0.95);
        
        set(h_(1),'Color',[1 0 0], 'LineStyle','-', 'LineWidth',2, 'Marker',...
            'none', 'MarkerSize',6);
        coeffs = coeffvalues(cf_)
        checkbutton=1;
        ylim([0,1])
        while checkbutton ~= 2;
            [x1,y1,button] = ginput(1);
            if button ==1
                chosencellpos(ind) = pos;
                chosencellno(ind) = cellno;
                hills(ind) = coeffs(3);
                checkbutton = 2;
                ind = ind+1;
                ind
            elseif button == 3
                checkbutton = 2;
            end
            
        end
        hold off
    end
end
if ind > 50
    save(['arrestCln3' '_' 'hills'],'hills','chosencellpos','chosencellno')
end
% xvals = xvals(yvals>0);
% yvals = yvals(yvals>0);
% plot(xvals, yvals,'x')
% maximum = 1e5;
% minimum = 250;
% binNo = 50;
%
% bins= linspace(minimum, maximum, binNo);
%
% [h,specBins] = histc(xvals, bins);
%
% for i = 1:binNo
%     binvalues    = yvals(specBins == i);
%     binmeds(i)     = median(binvalues);
% end
% hold on
% plot(bins,binmeds,'r')
%%
clear all; clc; close all;
backgroundfile = 'OA_041315_OA064_11Abar1del_6nM_pos_no_';
load arrestCln3_hills
indm = 1;
prevpos = chosencellpos(1);
fnorm = 0.20; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');
cln3oscillationdata = zeros(2,2);
pos = 1;
load([backgroundfile int2str(pos) '_re_exp_vol_analysis']);
load([backgroundfile int2str(pos) '_re_exp_vol'])
% background subtraction
area = all_obj.tot_nucl_areaG + all_obj.tot_cyt_areaG;
background1 = median(all_obj.med_backgr_gfp)/73.*(80.7082+area(area < 1135).*0.0030);%%%background subt.
j2 = logical((area >= 1135).*(area <1820));
background2 = median(all_obj.med_backgr_gfp)/73.*(88.2416-area(j2).*0.0036);
background3 = median(all_obj.med_backgr_gfp)/73.*(77.6698+area(area >= 1820).*0.0020);
all_obj.nuc_whi5G(area < 1135) = all_obj.nuc_whi5G(area < 1135) - background1;
all_obj.nuc_whi5G(j2) = all_obj.nuc_whi5G(j2)-background2; %%%%%%%%%%%%
all_obj.nuc_whi5G(area >= 1820) = all_obj.nuc_whi5G(area >= 1820)-background3; %%%%%%%%%%%%
all_obj.cyt_whi5G(area < 1135) = all_obj.cyt_whi5G(area < 1135) - background1;
all_obj.cyt_whi5G(j2) = all_obj.cyt_whi5G(j2)-background2; %%%%%%%%%%%%
all_obj.cyt_whi5G(area >= 1820) = all_obj.cyt_whi5G(area >= 1820)-background3;
% Cln3signal
volume = all_obj.volume_fl;
volume(find( (isinf(volume)+isnan(volume))>0) )=0;


Cln3tot = all_obj.tot_nucl_areaG.*all_obj.nuc_whi5G + all_obj.tot_cyt_areaG.*all_obj.cyt_whi5G;
Cln3conc = Cln3tot./volume;
Cln3conc(find( (isinf(Cln3conc)+isnan(Cln3conc))>0) )=0;
for i = 1:size(chosencellno,2)
    pos = chosencellpos(i);
    if pos ~= prevpos
        load([backgroundfile int2str(pos) '_re_exp_vol_analysis']);
        load([backgroundfile int2str(pos) '_re_exp_vol'])
        % background subtraction
        area = all_obj.tot_nucl_areaG + all_obj.tot_cyt_areaG;
        background1 = median(all_obj.med_backgr_gfp)/73.*(80.7082+area(area < 1135).*0.0030);%%%background subt.
        j2 = logical((area >= 1135).*(area <1820));
        background2 = median(all_obj.med_backgr_gfp)/73.*(88.2416-area(j2).*0.0036);
        background3 = median(all_obj.med_backgr_gfp)/73.*(77.6698+area(area >= 1820).*0.0020);
        all_obj.nuc_whi5G(area < 1135) = all_obj.nuc_whi5G(area < 1135) - background1;
        all_obj.nuc_whi5G(j2) = all_obj.nuc_whi5G(j2)-background2; %%%%%%%%%%%%
        all_obj.nuc_whi5G(area >= 1820) = all_obj.nuc_whi5G(area >= 1820)-background3; %%%%%%%%%%%%
        all_obj.cyt_whi5G(area < 1135) = all_obj.cyt_whi5G(area < 1135) - background1;
        all_obj.cyt_whi5G(j2) = all_obj.cyt_whi5G(j2)-background2; %%%%%%%%%%%%
        all_obj.cyt_whi5G(area >= 1820) = all_obj.cyt_whi5G(area >= 1820)-background3;
        % Cln3signal
        volume = all_obj.volume_fl;
        volume(find( (isinf(volume)+isnan(volume))>0) )=0;
        
        
        Cln3tot = all_obj.tot_nucl_areaG.*all_obj.nuc_whi5G + all_obj.tot_cyt_areaG.*all_obj.cyt_whi5G;
        Cln3conc = Cln3tot./volume;
        Cln3conc(find( (isinf(Cln3conc)+isnan(Cln3conc))>0) )=0;
        
        prevpos=pos;
    end
    
    whichrow= double(find((data(:,1) == chosencellno(i)) > 0));
    dataforchosencell = data(whichrow,:);
    if dataforchosencell(3)  == 0 % if mother cell
        close all
        hold all
        cellno = data(whichrow,1);
        filteredCln3 = filtfilt(b,a,Cln3conc(cellno,:));
        plot(filteredCln3);
        max_v = max(Cln3conc(cellno,:));
        plot([dataforchosencell(5) dataforchosencell(5)],[0 max_v],'color',[1 0 1],'linewidth',2)
        plot([dataforchosencell(6) dataforchosencell(6)],[0 max_v],'color',[1 0 1],'linewidth',2)
        plot([dataforchosencell(7) dataforchosencell(7)],[0 max_v],'color',[1 0 1],'linewidth',2)
        plot([dataforchosencell(8) dataforchosencell(8)],[0 max_v],'color',[1 0 1],'linewidth',2)
        display(dataforchosencell)
        checkbutton=1;
        while checkbutton ~= 2;
            
            [x1,y1,button] = ginput(1);
            if button ==1
                cln3oscillationdata(indm,1)= pos;
                cln3oscillationdata(indm,2)= cellno;
                display('choose maximum point of pre-arrest')
                [x2,y2,button2] = ginput(1);
                plot([x2 x2],[0 max_v],'color',[1 1 0],'linewidth',2)
                x2 = round(x2);
                cln3oscillationdata(indm,3)= mean(filteredCln3(x2-1:x2+1));
                
                display('choose pre-arrest low begin')
                [x3,y3,button3] = ginput(1);
                x3 = round(x3);
                plot([x3 x3],[0 max_v],'color',[1 1 0],'linewidth',2)
                
                display('choose pre-arrest low end')
                [x4,y4,button3] = ginput(1);
                x4 = round(x4);
                plot([x4 x4],[0 max_v],'color',[1 1 0],'linewidth',2)
                cln3oscillationdata(indm,4) = mean(filteredCln3(x3:x4));
                
                display('choose arrest high begin')
                [x5,y5,button3] = ginput(1);
                x5 = round(x5);
                plot([x5 x5],[0 max_v],'color',[1 1 0],'linewidth',2)
                
                display('choose arrest high end')
                [x6,y6,button3] = ginput(1);
                x6 = round(x6);
                plot([x6 x6],[0 max_v],'color',[1 1 0],'linewidth',2)
                cln3oscillationdata(indm,5) = mean(filteredCln3(x5:x6));
                
                
                indm = indm + 1;
                checkbutton = 2;
            elseif button == 3
                checkbutton = 2;
            end
            
        end
    end
end

if size(cln3oscillationdata,1) > 20
    save(['arrestCln3' '_' 'oscillations'],'cln3oscillationdata')
end

