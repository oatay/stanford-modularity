clear all; clc; close all;
backgroundfile = 'OA_012414_OA037and38_12nMaFcontrol_pos_no_';
data_end_point = 70;
medbckgrGFP = zeros(1,1);
j = 1;
xvals = zeros(1,1);
yvals = zeros(1,1);


for pos = [2 3 6 10 11 12]
    load([backgroundfile int2str(pos) '_re_exp_analysis']);
    load([backgroundfile int2str(pos) '_re_exp']);
    totautofluorescence = all_obj.tot_nucl_areaR.*all_obj.nuc_Far1 + all_obj.tot_cyt_areaR.*all_obj.cyt_Far1;
    area = all_obj.tot_nucl_areaR+all_obj.tot_cyt_areaR;
    medbckgrGFP(j) = median(all_obj.med_backgr_Far1);
    j = j +1 ;
    for i = 1:size(data,1)
        startp = round(data(i,5));
        endp = round(data(i,6));
        cellno = data(i,1);
        if pos == 3
            cons = 45/42;
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
ppt2 = 49;
p2 = polyfit(bins((ppt):ppt2), binmeds((ppt):ppt2),1);
p1 = polyfit(bins(3:ppt), binmeds(3:ppt),1);
p1y = p1(1).*bins(1:ppt)+p1(2);
p2y = p2(1).*bins((ppt):ppt2)+p2(2);
plot(bins(1:ppt),p1y)
plot(bins((ppt):ppt2),p2y)

