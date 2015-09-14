clear all; clc; close all;


backgroundfile = 'AD_wtc5p_30_17c_oct31_2013_pos_no_';
xvals = zeros(1,1);
yvals = zeros(1,1);

fnorm = 0.10; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');
chosencellpos = zeros(1,1);
chosencellno = zeros(1,1);
ind = 1;
m = 1;
d = 1;
hills = zeros(1,1);
for pos = [1:3 5 7:12 14 15]
    load([backgroundfile int2str(pos) '_vol_analysis']);
    load([backgroundfile int2str(pos) '_vol']);
    
    % background subtraction
    volume = all_obj.volume;
    area= all_obj.area;
    background1 = median(all_obj.med_backgr_gfp)/45.*(74.3126+area(area < 3215).*0.0072);%%%background subt.
    background3 = median(all_obj.med_backgr_gfp)/45.*(100.0599-area(area >= 3215).*0.0011);

    all_obj.nuc_whi5R(area < 3215) = all_obj.nuc_whi5R(area < 3215) - background1;
    all_obj.nuc_whi5R(area >= 3215) = all_obj.nuc_whi5R(area >= 3215)-background3; %%%%%%%%%%%%
    all_obj.cyt_whi5R(area < 3215) = all_obj.cyt_whi5R(area < 3215) - background1;
    all_obj.cyt_whi5R(area >= 3215) = all_obj.cyt_whi5R(area >= 3215)-background3;
    
   
    volume(find( (isinf(volume)+isnan(volume))>0) )=0;
    
    
    GFPtot = all_obj.tot_nucl_areaR.*all_obj.nuc_whi5R + all_obj.tot_cyt_areaR.*all_obj.cyt_whi5R;
    GFPconc = GFPtot./volume;
    GFPconc(find( (isinf(GFPconc)+isnan(GFPconc))>0) )=0;
    
    for i = 1:size(data,1)
        close all
        hold all
        
        startp = round(data(i,5));
        reentry = round(data(i,6));
        endp = round(data(i,7));
        cellno = data(i,1);
        xvals = [xvals volume(cellno, startp:endp)];
        yvals = [yvals GFPconc(cellno,startp:endp)];

        stopp = min(endp+5,170);
        volume_sm = volume(cellno,reentry-17:stopp-7);%filtfilt(b,a,volume(cellno, startp:endp));
        volstart = median(volume(cellno,startp-10:startp));
        volume_sm_norm = volume_sm./volstart;
        clb5conc_sm = GFPconc(cellno,reentry-10:stopp);%filtfilt(b,a,Cln3conc(cellno,startp:endp));
        sortedclb5 = sort(clb5conc_sm,'ascend');
        backgr = median(sortedclb5(1:10));
        clb5conc_sm = clb5conc_sm - backgr;
        sortedclb5 = sort(clb5conc_sm,'descend');
        clb5conc_sm_norm = clb5conc_sm./median(sortedclb5(1:10));
        
        plot(filtfilt(b,a,volume_sm_norm),filtfilt(b,a,clb5conc_sm_norm))
        volume_sm_norm(volume_sm_norm <= 1) = NaN;
        clb5conc_sm_norm(clb5conc_sm_norm <= 0) = NaN;
        
        plot(volume_sm_norm,clb5conc_sm_norm,'x')
        
        fo_ = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',[0 0  0],'Upper',[15   2 200]);
        ok_ = isfinite(volume_sm_norm) & isfinite(clb5conc_sm_norm);
        st_ = [0.037051126859541816 0.13223694988779544 0.2171532878989586];
        set(fo_,'Startpoint',st_);
        ft_ = fittype('a.*((vol-1).^n)./(Km.^n + (vol-1).^n)',...
            'dependent',{'y'},'independent',{'vol'},...
            'coefficients',{'Km', 'a', 'n'});
        try
            cf_ = fit(volume_sm_norm(ok_)',clb5conc_sm_norm(ok_)',ft_,fo_);
        catch err
            display(err)
        end
        
        h_ = plot(cf_,'fit',0.95);
        
        set(h_(1),'Color',[1 0 0], 'LineStyle','-', 'LineWidth',2, 'Marker',...
            'none', 'MarkerSize',6);
        coeffs = coeffvalues(cf_)
        checkbutton=1;
        ylim([0,1.2])
        while checkbutton ~= 2;
            [x1,y1,button] = ginput(1);
            if button ==1
                chosencellpos(ind) = pos;
                chosencellno(ind) = cellno;
                hills(ind) = coeffs(3);
                checkbutton = 2;
                ind = ind+1;
                if data(i,3) == 0
                    hillsm(m) = coeffs(3);
                    m = m + 1;
                else
                    hillsd(d) = coeffs(3);
                    d = d + 1;
                end
                ind
            elseif button == 3
                checkbutton = 2;
            end
            
        end
        hold off
    end
end
if ind > 20
    save(['arrestClb5' '_' 'hills'],'hills','hillsm','hillsd','chosencellpos','chosencellno')
end