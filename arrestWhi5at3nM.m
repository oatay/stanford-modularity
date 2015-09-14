clc;clear all;close all;
backgroundfile = 'AD_wtc5p_30_17c_oct31_2013_pos_no_';
xd = zeros(1,1); % daughter volume
yd = zeros(1,1); % daughter concenration
td = zeros(1,1); % daughter time

xm = zeros(1,1);
ym = zeros(1,1);
tm = zeros(1,1); % mother time

fnorm = 0.10; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');
totm = 0;
totd = 0;
totre = 0;
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
    
    
    Whi5tot = all_obj.tot_nucl_areaR.*all_obj.nuc_whi5R + all_obj.tot_cyt_areaR.*all_obj.cyt_whi5R;
    Whi5conc = Whi5tot./volume;
    Whi5conc(find( (isinf(Whi5conc)+isnan(Whi5conc))>0) )=0;
    
    for i = 1:size(data,1)
        
        startp = round(data(i,5));
        endp = round(data(i,6));
        cellno = data(i,1);
        
        
        if data(i,3) == 1 %daughter cell
            totd = totd + 1;
            xd = [xd volume(cellno, startp:endp)];
            yd = [yd Whi5conc(cellno,startp:endp)];
            td = [td 3.*(0:endp-startp)];
            coeff = polyfit(3.*(0:endp-startp),volume(cellno, startp:endp),1);
            ppd(totd) = coeff(1);
            ppd2(totd)= coeff(2);
            
            expft=fittype('exp1');
            cf=fit(3.*(0:endp-startp)',volume(cellno, startp:endp)',expft);
            coeffab = coeffvalues(cf);
            ppdexpa(totd) = coeffab(1);
            ppdexpb(totd) = coeffab(2);
            

        else % mother cells
            totm = totm + 1;
            xm = [xm volume(cellno, startp:endp)];
            ym = [ym Whi5conc(cellno, startp:endp)];
            tm = [tm 3.*(0:endp-startp)];
            coeff = polyfit(3.*(0:endp-startp),volume(cellno, startp:endp),1);
            ppm(totm) = coeff(1);
            ppm2(totm) = coeff(2);
            
            expft=fittype('exp1');
            cf=fit(3.*(0:endp-startp)',volume(cellno, startp:endp)',expft);
            coeffab = coeffvalues(cf);
            ppmexpa(totm) = coeffab(1);
            ppmexpb(totm) = coeffab(2);

        end
        % re-entering cells
    end
end


xm = xm(2:end);
ym = ym(2:end);
tm = tm(2:end);
xd = xd(2:end);
yd = yd(2:end);
td = td(2:end);


% xm = xm(xm<4.5e4);
% ym = ym(xm<4.5e4);
% 
% xm = xm(ym<1.6);
% ym = ym(ym<1.6);


xd = xd(yd<4);
yd = yd(yd<4);
td = td(yd<4);


figure(1)
hold on
%plot(xm,ym,'x')
[bin1 binmeds1 binstds1] = makebins(xm,ym,0.8e4,3.5e4,40);
ciplot(binmeds1-binstds1, binmeds1+binstds1, bin1,'r')
ylim([0 4])
xlim([0.8e4 3.5e4])

plot(bin1,binmeds1,'LineWidth',3)
hold off

figure(2)
hold on
%plot(xd,yd,'x')
[bin2 binmeds2 binstds2] = makebins(xd,yd,0.5e4,4e4,40);
ciplot(binmeds2-binstds2, binmeds2+binstds2, bin2,'r')
plot(bin2,binmeds2,'LineWidth',3)
ylim([0 4])
xlim([0.4e4 4e4])
hold off

figure(3)
hold on
plot(td,xd,'x')
[bin4 binmeds4 binstds4] = makebins(td,xd,0,400,40);
%ciplot(binmeds4-binstds4, binmeds4+binstds4, bin4,'r')
%plot(bin4,binmeds4,'LineWidth',3)
plot(bin4,bin4.*median(ppd)+median(ppd2));
hold off

figure(4)
hold on
plot(tm,xm,'x')
[bin5 binmeds5 binstds5] = makebins(tm,xm,0,400,40);
%ciplot(binmeds5-binstds5, binmeds5+binstds5, bin5,'r')
%plot(bin5,binmeds5,'LineWidth',3)
plot(bin4,bin4.*median(ppm)+median(ppm2));
hold off

%109.9257 volume/min for mothers; 113.1496 for daughters. %medians
%median(1/b): 150.6351*ln2 doubling time for daughters; 177.6760*ln2 for mothers
% 104.4123 for daughters; 123.1556 for mothers

save([backgroundfile,'arrestWhi5at3nM'],'xd','yd','td','xm','ym','tm','ppd','ppd2',...
    'ppm','ppm2','ppdexpa','ppdexpb','ppmexpa','ppmexpb');
