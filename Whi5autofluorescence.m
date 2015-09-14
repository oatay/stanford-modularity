clear all; clc; close all;
% Whi5-mkok autofluorescence - 750 ms
backgroundfile = 'OA_011514_JS163-8dand29_25_5-120nMto3nM_pos_no_';
medbckgrGFP = zeros(1,1);
j = 1;
xvals = zeros(1,1);
yvals = zeros(1,1);



for pos = [4 6 7 8]
    load([backgroundfile int2str(pos) '_1_analysis']);
    load([backgroundfile int2str(pos) '_1']);
    only = size(all_obj.volume,1);
    totautofluorescence = all_obj.tot_nucl_area_w5.*all_obj.nuc_whi5R + all_obj.tot_cyt_area_w5(1:only,:).*all_obj.cyt_whi5R;
    area = all_obj.tot_nucl_area_w5+all_obj.tot_cyt_area_w5(1:only,:);
    medbckgrGFP(j) = median(all_obj.med_backgr_w5);
    j = j +1 ;
    for i = 1:size(data,1)
        startp = round(data(i,7));
        endp = round(data(i,8));
        cellno = data(i,1);
        if pos == 8
            cons = 29/27;
        else
            cons = 1;
        end
        xvals = [xvals all_obj.area(cellno, startp:endp)];
        yvals = [yvals cons.*(totautofluorescence(cellno,startp:endp)./area(cellno, startp:endp))];

    end
end
xvals = xvals(yvals>0);
yvals = yvals(yvals>0);
plot(xvals, yvals,'x')
maximum = 6000; 
minimum = 750;
binNo = 50;

bins= linspace(minimum, maximum, binNo);

[h,specBins] = histc(xvals, bins);

for i = 1:binNo
    binvalues    = yvals(specBins == i);
    binmeds(i)     = median(binvalues);
end
hold on
plot(bins,binmeds,'r')

ppt = 24;
ppt2 = 33;
ppt3 = 49;
p2 = polyfit(bins((ppt):ppt2), binmeds((ppt):ppt2),1);
p1 = polyfit(bins(1:ppt), binmeds(1:ppt),1);
p3 = polyfit(bins((ppt2):ppt3), binmeds((ppt2):ppt3),1);
p1y = p1(1).*bins(1:ppt)+p1(2);
p2y = p2(1).*bins((ppt):ppt2)+p2(2);
p3y = p3(1).*bins((ppt2):ppt3)+p3(2);
plot(bins(1:ppt),p1y)
plot(bins((ppt):ppt2),p2y)
plot(bins((ppt2):ppt3),p3y)
