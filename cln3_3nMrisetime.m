clear all; clc; close all;
backgroundfile = 'OA_070815_11Abar1d_refurb_3nM6min_pos_no_';
fnorm = 0.10; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');

indd = 1;
indm = 1;
ind = 1;
for pos = [2 3 4 5 7 8 10 11]
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
        close all
        hold all
        startp = data(i,7);
        endp = data(i,8);
        cellno = data(i,1);
        
        from = floor(startp);
        to = ceil(endp);
        
        cln3data = cln3conc(cellno,from:to);
        smoothcln3 = filtfilt(b,a,cln3data);
        
        midcln3= (max(smoothcln3) + min(smoothcln3))/2;
        plot((from:to).*6, cln3data)
        plot((from:to).*6, smoothcln3)

        
        mincln3values(ind) = min(smoothcln3);
        maxcln3values(ind) = max(smoothcln3);
        ind = ind + 1;
        midtimepoint = interp1(smoothcln3,from:to,midcln3);

        pos
        cellno
        %pause(10)
        if data(i,3) == 1
            timetohalfmaxd(indd) = (midtimepoint - startp)*6;
            indd = indd + 1;
        else
            timetohalfmaxm(indm) = (midtimepoint - startp)*6;
            indm = indm + 1;
        end
    end
end
median(timetohalfmaxd)
std(bootstrp(1000,@median,timetohalfmaxd))

