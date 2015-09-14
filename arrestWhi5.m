%mother cells
%daughter cells
%re-entry


clear all; clc; close all;
backgroundfile = 'OA_060514_OA040and42_4-5nM_pos_no_';
xd = zeros(1,1);
yd = zeros(1,1);
td = zeros(1,1); % daughter time

xm = zeros(1,1);
ym = zeros(1,1);
tm = zeros(1,1); % mother time


xre = zeros(1,1);
yre = zeros(1,1);

fnorm = 0.10; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');
chosencellpos = zeros(1,1);
chosencellno = zeros(1,1);
ind = 1;
totm = 0;
totd = 0;
totre = 0;

for pos = [2 3 4 5 8 9 11]
    load([backgroundfile int2str(pos) '_re_exp_analysis']);
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
    
    
    Whi5tot = all_obj.tot_nucl_area_w5.*all_obj.nuc_whi5R + all_obj.tot_cyt_area_w5(1:only,:).*all_obj.cyt_whi5R;
    Whi5conc = Whi5tot./volume;
    Whi5conc(find( (isinf(Whi5conc)+isnan(Whi5conc))>0) )=0;
    
    for i = 1:size(data,1)
        
        startp = round(data(i,5));
        endp = round(data(i,6));
        cellno = data(i,1);
        
        
        if data(i,3) == 1 %daughter cell
            totd = totd + 1;
            %plot(6.*(0:endp-startp),volume(cellno, startp:endp))
            %pause(0.5)
            xd = [xd volume(cellno, startp:endp)];
            yd = [yd Whi5conc(cellno,startp:endp)];
            td = [td 6.*(0:endp-startp)];
            coef = polyfit(6.*(0:endp-startp),volume(cellno, startp:endp),1);
            ppd(totd) = coef(1);
            ppd2(totd) = coef(1);
            
            expft=fittype('exp1');
            cf=fit(6.*(0:endp-startp)',volume(cellno, startp:endp)',expft);
            coeffab = coeffvalues(cf);
            ppdexpa(totd) = coeffab(1);
            ppdexpb(totd) = coeffab(2);
            
        else % mother cells
            totm = totm + 1;
            xm = [xm volume(cellno, startp:endp)];
            ym = [ym Whi5conc(cellno, startp:endp)];
            tm = [tm 6.*(0:endp-startp)];
            coef = polyfit(6.*(0:endp-startp),volume(cellno, startp:endp),1);
            ppm(totm) = coef(1);
            ppm2(totm) = coef(2);
            
            expft=fittype('exp1');
            cf=fit(6.*(0:endp-startp)',volume(cellno, startp:endp)',expft);
            coeffab = coeffvalues(cf);
            ppmexpa(totm) = coeffab(1);
            ppmexpb(totm) = coeffab(2);
        end
        % re-entering cells
        if endp < 75
            totre = totre + 1;
            xre = [xre volume(cellno,endp:83)];
            yre = [yre Whi5conc(cellno,endp:83)];
            %plot(Whi5conc(cellno,endp:83))
            %pause(1)
        end
        
        
    end
end

xre = xre(2:end);
yre = yre(2:end);
xm = xm(2:end);
ym = ym(2:end);
xd = xd(2:end);
yd = yd(2:end);
tm = tm(2:end);
td = td(2:end);

xm = xm(xm<4.5e4);
ym = ym(xm<4.5e4);
tm = tm(xm<4.5e4);

xm = xm(ym<1.6);
ym = ym(ym<1.6);
tm = tm(ym<1.6);


xd = xd(yd<1.6);
yd = yd(yd<1.6);
td = td(yd<1.6);

figure(1)
hold on
plot(xm,ym,'x')
[bin1 binmeds1 binstds1] = makebins(xm,ym,0.5e4,4.5e4,40);
ciplot(binmeds1-binstds1, binmeds1+binstds1, bin1,'r')
ylim([0 2])
xlim([0.4e4 4.5e4])

plot(bin1,binmeds1,'LineWidth',3)
hold off

figure(2)
hold on
plot(xd,yd,'x')
[bin2 binmeds2 binstds2] = makebins(xd,yd,0.5e4,4.5e4,40);
ciplot(binmeds2-binstds2, binmeds2+binstds2, bin2,'r')
plot(bin2,binmeds2,'LineWidth',3)
ylim([0 2])
xlim([0.4e4 4.5e4])
hold off

figure(3)
hold on
plot(xre,yre,'x')
[bin3 binmeds3 binstds3] = makebins(xre,yre,0.5e4,4.5e4,40);
ciplot(binmeds3-binstds3, binmeds3+binstds3, bin3,'r')
plot(bin3,binmeds3,'LineWidth',3)
ylim([0 2])
xlim([1.5e4 4.5e4])
hold off

figure(4)
hold on
plot(td,xd,'x')
[bin4 binmeds4 binstds4] = makebins(td,xd,0,400,40);
plot(bin4,median(ppdexpa).*exp(median(ppdexpb).*bin4))
%ciplot(binmeds4-binstds4, binmeds4+binstds4, bin4,'r')
%plot(bin4,binmeds4,'LineWidth',3)
hold off

figure(5)
hold on
plot(tm,xm,'x')
[bin5 binmeds5 binstds5] = makebins(tm,xm,0,400,40);
%ciplot(binmeds5-binstds5, binmeds5+binstds5, bin5,'r')
%plot(bin5,binmeds5,'LineWidth',3)
plot(bin5,median(ppmexpa).*exp(median(ppmexpb).*bin5))
plot(bin5,median(ppm2)+bin5.*median(ppm))
hold off
% 67.0196 for mothers; 79.6352 for daughters; linear fit medians
% median(1/b) doubling times: 231.6890*ln2 for mothers; 204.7307*ln2 for daughters
% 160.59 min for mothers; 141.9085 for daughters
save([backgroundfile,'arrestWhi5_OA040_at4_5nM'],'xd','yd','td','xm','ym','tm','ppd','ppd2',...
    'ppm','ppm2','ppdexpa','ppdexpb','ppmexpa','ppmexpb');
