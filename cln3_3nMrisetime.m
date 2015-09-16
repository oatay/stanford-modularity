clear all; clc; close all;
backgroundfile = 'OA_070815_11Abar1d_refurb_3nM6min_pos_no_';
fnorm = 0.10; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');

indd = 1;
indm = 1;
ind = 1;
allcln3alignedd = cell(1);
allxdatad = cell(1);
allcln3alignedm = cell(1);
allxdatam = cell(1);

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
        %close all
        hold all
        startp = data(i,7);
        endp = data(i,8);
        cellno = data(i,1);
        
        from = floor(startp);
        to = ceil(endp);
        
        cln3data = cln3conc(cellno,from:to);
        smoothcln3 = filtfilt(b,a,cln3data);
        
        midcln3= (max(smoothcln3) + min(smoothcln3))/2;
        %plot(cln3data)
        

        
        mincln3values(ind) = min(smoothcln3);
        maxcln3values(ind) = max(smoothcln3);
        ind = ind + 1;
        midtimepoint = interp1(smoothcln3,from:to,midcln3);

        pos
        cellno
        %pause(10)
        if data(i,3) == 1
            timetohalfmaxd(indd) = (midtimepoint - startp)*6;
            allcln3alignedd{indd} = smoothcln3;
            allxdatad{indd} = ((from:to)-from).*6;
            indd = indd + 1;
        else
            timetohalfmaxm(indm) = (midtimepoint - startp)*6;
            allcln3alignedm{indm} = smoothcln3;
            allxdatam{indm} = ((from:to)-from).*6;
            indm = indm + 1;
        end
    end
end
median(timetohalfmaxd)
std(bootstrp(1000,@median,timetohalfmaxd))

if size(allcln3alignedd,2) > 15
    save([backgroundfile,'alignedCln3'],'timetohalfmaxd','timetohalfmaxm',...
        'allcln3alignedd','allxdatad','allcln3alignedm','allxdatam')
end

%%
for i = 1:size(allcln3alignedd,2)
    allcln3alignedd{i} = allcln3alignedd{i}./median(allcln3alignedd{i});
    plot(allxdata{i}, allcln3alignedd{i})
    curx = allxdatad{i}(1:end);
    cury = allcln3alignedd{i}(1:end);
    td = [td curx];
    yd = [yd cury];
end

td = td(2:end);
yd = yd(2:end);

[bin binmeds binstds] = makebins(td,yd,0,300,30);
figure(2)
hold on
ciplot((binmeds-binstds), (binmeds+binstds), bin,'r')
plot(bin,binmeds,'LineWidth',3)
ylim([0 0.5])
xlim([0 200])
hold off