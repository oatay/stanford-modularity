%mother cells
%daughter cells
%re-entry


clear all; clc; close all;
backgroundfile = 'OA_061015_OA040_3nM6min_pos_no_';
% for reentering cells
xd = zeros(1,1);
yd = zeros(1,1);
td = zeros(1,1); % daughter time

xm = zeros(1,1);
ym = zeros(1,1);
tm = zeros(1,1); % mother time

% not reentering cells
xdnre = zeros(1,1);
ydnre = zeros(1,1);
tdnre = zeros(1,1);

xmnre = zeros(1,1);
ymnre = zeros(1,1);
tmnre = zeros(1,1);

zerovol = zeros(1,1);

fnorm = 0.15; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');
chosencellpos = zeros(1,1);
chosencellno = zeros(1,1);
ind = 1;
totm = 0;
totd = 0;
totmnre = 0;
totdnre = 0;


for pos = [2 3 4 5 6 7 10 12 13]
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
    
    
    Whi5tot = all_obj.tot_nucl_area_w5.*all_obj.nuc_whi5R + all_obj.tot_cyt_area_w5(1:only,:).*all_obj.cyt_whi5R;
    Whi5conc = Whi5tot./volume;
    Whi5conc(find( (isinf(Whi5conc)+isnan(Whi5conc))>0) )=0;
    
    for i = 1:size(data,1)
        
        startp = round(data(i,5));
        endp = round(data(i,6));
        cellno = data(i,1);
        % reentering cells
        if endp < 135
            if data(i,3) == 1 %daughter cell
                totd = totd + 1;
                %plot(6.*(0:endp-startp),volume(cellno, startp:endp))
                %pause(0.5)
                xd = [xd volume(cellno, startp:endp)];
                yd = [yd Whi5conc(cellno,startp:endp)];
                td = [td 6.*(0:endp-startp)];
                darresttime(totd) = (endp-startp)*6;
                
                coef = polyfit(6.*((startp:endp)-startp),volume(cellno, (startp):endp),1);
                smoothvol = filtfilt(b,a,volume(cellno, (startp):endp));
                % plot(((startp:endp)-startp),smoothvol)
                %pause(2)
                doubletime(totd) = interp1(smoothvol,((startp:endp)-startp),2*mean(smoothvol(1:3)));
                zerovol = [zerovol smoothvol];
                
                hold all
                smoothened = filtfilt(b,a,volume(cellno, 1:endp));
                normfv = smoothvol(1);
                plot(6.*(0:endp-1),volume(cellno, 1:endp)./normfv)
                plot(6.*(0:endp-1),smoothened./normfv)
                ylim([0 5])
                xlim([0 400])
                pause(10)
                hold off
                close all
                
                ppd(totd) = coef(1);
                ppd2(totd) = coef(2);
                halfmaxtimesd(totd) = coef(2)/coef(1);
                
                % hold on
                %plot((1:140).*6,Whi5conc(cellno,:))
                %plot(6.*((startp:endp)-startp),volume(cellno,startp:endp))
                expft=fittype('exp1');
                cf=fit(6.*(0:endp-startp)',volume(cellno, startp:endp)',expft);
                coeffab = coeffvalues(cf);
                ppdexpa(totd) = coeffab(1);
                ppdexpb(totd) = coeffab(2);
                %plot((1:140).*6,volume(cellno,:))
                
%                 linfit = fittype('poly1');
%                 cf=fit(6.*(0:endp-startp)',volume(cellno, startp:endp)',linfit);
%                 coeffab = coeffvalues(cf);
%                 ppdlina(totd) = coeffab(1);
%                 ppdlinb(totd) = coeffab(2);
        
            else % mother cells
                totm = totm + 1;
                xm = [xm volume(cellno, startp:endp)];
                ym = [ym Whi5conc(cellno, startp:endp)];
                tm = [tm 6.*(0:endp-startp)];
                marresttime(totm) = (endp-startp)*6;
                coef = polyfit(6.*(0:endp-startp),volume(cellno, startp:endp),1);
                ppm(totm) = coef(1);
                ppm2(totm) = coef(2);

                expft=fittype('exp1');
                cf=fit(6.*(0:endp-startp)',volume(cellno, startp:endp)',expft);
                coeffab = coeffvalues(cf);
                ppmexpa(totm) = coeffab(1);
                ppmexpb(totm) = coeffab(2);
                
%                 linfit = fittype('poly1');
%                 cf=fit(6.*(0:endp-startp)',volume(cellno, startp:endp)',linfit);
%                 coeffab = coeffvalues(cf);
%                 ppmlina(totm) = coeffab(1);
%                 ppmlinb(totm) = coeffab(2);
            end
        % not reentering cells (for more than five hours)
        elseif ((endp > 135) && (startp < 90))
            if data(i,3) == 1 %daughter cell
                totdnre = totdnre + 1;
                %plot(6.*(0:endp-startp),volume(cellno, startp:endp))
                %pause(0.5)
                xdnre = [xdnre volume(cellno, startp:endp)];
                ydnre = [ydnre Whi5conc(cellno,startp:endp)];
                tdnre = [tdnre 6.*(0:endp-startp)];

                coef = polyfit(6.*(0:endp-startp),volume(cellno, startp:endp),1);
                ppdnre(totdnre) = coef(1);
                ppd2nre(totdnre) = coef(1);
                
                expft=fittype('exp1');
                cf=fit(6.*(0:endp-startp)',volume(cellno, startp:endp)',expft);
                coeffab = coeffvalues(cf);
                ppdexpanre(totdnre) = coeffab(1);
                ppdexpbnre(totdnre) = coeffab(2);
                
            else % mother cells
                totmnre = totmnre + 1;
                xmnre = [xmnre volume(cellno, startp:endp)];
                ymnre = [ymnre Whi5conc(cellno, startp:endp)];
                tmnre = [tmnre 6.*(0:endp-startp)];
                coef = polyfit(6.*(0:endp-startp),volume(cellno, startp:endp),1);
                ppmnre(totmnre) = coef(1);
                ppm2nre(totmnre) = coef(2);
                
                expft=fittype('exp1');
                cf=fit(6.*(0:endp-startp)',volume(cellno, startp:endp)',expft);
                coeffab = coeffvalues(cf);
                ppmexpanre(totmnre) = coeffab(1);
                ppmexpbnre(totmnre) = coeffab(2);
            end
        end
        
        
    end
end

close all;
xm = xm(2:end); xmnre = xmnre(2:end);
ym = ym(2:end); ymnre = ymnre(2:end);
xd = xd(2:end); xdnre = xdnre(2:end);
yd = yd(2:end); ydnre = ydnre(2:end);
tm = tm(2:end); tmnre = tmnre(2:end);
td = td(2:end); tdnre = tdnre(2:end);
zerovol = zerovol(2:end);

xm = xm(xm<4.5e4); xmnre = xmnre(xmnre<4.5e4);
ym = ym(xm<4.5e4); ymnre = ymnre(xmnre<4.5e4);
tm = tm(xm<4.5e4); tmnre = tmnre(tmnre<4.5e4);

xm = xm(ym<1.5); xmnre = xmnre(ymnre<1.5);
ym = ym(ym<1.5); ymnre = ymnre(ymnre<1.5);
tm = tm(ym<1.5); tmnre = tmnre(ymnre<1.5);

xd = xd(xd<4.5e4); xdnre = xdnre(xdnre<4.5e4);
yd = yd(xd<4.5e4); ydnre = ydnre(xdnre<4.5e4);
td = td(xd<4.5e4); tdnre = tdnre(xdnre<4.5e4);
zerovol = zerovol(xd<4.5e4);

xd = xd(yd<1.5); xdnre = xdnre(ydnre<1.5);
yd = yd(yd<1.5); ydnre = ydnre(ydnre<1.5);
td = td(yd<1.5); tdnre = tdnre(ydnre<1.5);
zerovol = zerovol(yd<1.5);

figure(1)
hold on
plot(xm,ym,'x')
[bin1 binmeds1 binstds1] = makebins(xm,ym,0.5e4,4.5e4,30);
ciplot(binmeds1-binstds1, binmeds1+binstds1, bin1,'r')
ylim([0 2])
xlim([0.4e4 4.5e4])

plot(bin1,binmeds1,'LineWidth',3)
hold off

figure(2)
hold on
%plot(xd,yd,'x')
[bin2 binmeds2 binstds2] = makebins(xd,yd,0.5e4,3.5e4,30);
bin2 = bin2./min(bin2);
normd = binmeds2(1);
binmeds2 = binmeds2./normd;
bindstds2 = binstds2./normd;

ciplot(binmeds2-binstds2, binmeds2+binstds2, bin2,'r')
plot(bin2,binmeds2,'LineWidth',3)
ylim([0 1.5])
xlim([1 7])
hold off

figure(3)
hold on
plot(xdnre,ydnre,'x')
[bin3 binmeds3 binstds3] = makebins(xdnre,ydnre,0.5e4,4.5e4,30);
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
plot(bin4,mean(ppd2)+bin4.*mean(ppd),'r')
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


% save([backgroundfile,'arrestWhi5_OA040_at3nM'],'xd','yd','td','xm','ym','tm','ppd','ppd2',...
%     'ppm','ppm2','ppdexpa','ppdexpb','ppmexpa','ppmexpb');

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

xvals = 0:1:300;
yvals = a.*exp(b.*xvals);
highyvals = (a+stda).*exp((b+stdb).*xvals);
lowyvals = (a-stda).*exp((b-stdb).*xvals);
normfac = min(yvals);
yvals = yvals./normfac;
highyvals = highyvals./normfac;
lowyvals = lowyvals./normfac;


lower_er = yvals - lowyvals;
upper_er = highyvals - yvals;
figure(6)
errorbar(xvals,yvals,lower_er,upper_er)
xlim([0 300])
  
%linear volume fit


med = @(x)median(x);    

a = median(ppd);
b = median(ppd2);
bootd3nMlina = bootstrp(1000,med,ppd);
bootd3nMlinb = bootstrp(1000,med,ppd2);
stda = std(bootd3nMlina);
stdb = std(bootd3nMlinb);

xvals = 150:1.5:600;
yvals = a+xvals.*b;
highyvals = (a+stda)+xvals.*(b+stdb);
lowyvals = (a-stda)+xvals.*(b-stdb);
normfac = min(yvals);
yvals = yvals./normfac;
highyvals = highyvals./normfac;
lowyvals = lowyvals./normfac;


lower_er = yvals - lowyvals;
upper_er = highyvals - yvals;
figure(7)
errorbar(xvals-150,yvals,lower_er,upper_er)
xlim([0 400]) 
  

figure(8)
hold on
td = td(zerovol < 5e4);
yd = yd(zerovol < 5e4);
zerovol = zerovol(zerovol < 5e4);
% plot(td,zerovol,'x')
[bin8 binmeds8 binstds8] = makebins(td,zerovol,0,300,100);
normf = binmeds8(1); 
ciplot((binmeds8-binstds8)./normf, (binmeds8+binstds8)./normf, bin8,'r')
plot(bin8,binmeds8./normf,'LineWidth',3)

figure(9)
hold on
%plot(td,yd,'x')
[bin9 binmeds9 binstds9] = makebins(td,yd,0,200,20);
bin9 = bin9;
normd = 1; % binmeds9(1);
binmeds9 = binmeds9./normd;
bindstds9 = binstds9./normd;

ciplot(binmeds9-binstds9, binmeds9+binstds9, bin9,'r')
plot(bin9,binmeds9,'LineWidth',3)
ylim([0 1])
xlim([0 200])
hold off