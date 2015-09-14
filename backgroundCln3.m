clear all; clc; close all;
backgroundfile = 'OA_051315_JS163_8d-4_Venus3min6nM_pos_no_';
medbckgrGFP = zeros(1,1);
j = 1;
xvals = zeros(1,1);
yvals = zeros(1,1);



for pos = [2 3 4 5 6 9 11 12]
    load([backgroundfile int2str(pos) '_re_exp_analysis']);
    load([backgroundfile int2str(pos) '_re_exp']);

    totautofluorescence = all_obj.cyt_whi5G.*all_obj.tot_cyt_areaG + ...
    all_obj.nuc_whi5G.*all_obj.nuc_whi5G;
    area = all_obj.tot_cyt_areaG + all_obj.tot_nucl_areaG;
    medbckgrGFP(j) = median(all_obj.med_backgr_gfp);
    j = j +1 ;
    for i = 1:size(data,1)
        startp = round(data(i,7));
        endp = round(data(i,8));
        cellno = data(i,1);
        xvals = [xvals area(cellno, startp:endp)];
        yvals = [yvals totautofluorescence(cellno,startp:endp)./area(cellno, startp:endp)];        
    end
end
xvals = xvals(yvals>0);
yvals = yvals(yvals>0);
plot(xvals, yvals,'x')
maximum = 2500; 
minimum = 250;
binNo = 50;

bins= linspace(minimum, maximum, binNo);

[h,specBins] = histc(xvals, bins);

for i = 1:binNo
    binvalues    = yvals(specBins == i);
    binmeds(i)     = median(binvalues);
end
hold all
plot(bins,binmeds,'r')

p1 = polyfit(bins(5:27), binmeds(5:27),1);
p2 = polyfit(bins(27:34), binmeds(27:34),1);
p3 = polyfit(bins(34:47), binmeds(34:47),0);
p1y = p1(1).*bins(5:27)+p1(2);
p2y = p2(1).*bins(27:34)+p2(2);
p3y = 0.*bins(34:47)+p3(1);
plot(bins(5:27),p1y)
plot(bins(27:34),p2y)
plot(bins(34:47),p3y)

% breaks = [1000 2000];
% pp = splinefit(bins(1:99), binmeds(1:99),breaks);
% xx = linspace(0,2500,50);
% yy = ppval(pp,xx);
% plot(xx,yy)
% 
% save([backgroundfile '_' 'background'],'pp','p1','p2','bins','binmeds')
