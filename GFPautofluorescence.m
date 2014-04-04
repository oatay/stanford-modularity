clear all; clc; close all;
backgroundfile = 'OA_012414_OA037and38_12nMaFcontrol_pos_no_';
data_end_point = 70;
analysisfile = 'OA_030714_OA039L381-L392-Longpulse_pos_no_';
m = 0;
p1s = zeros(1,1);
p2s = zeros(1,1);
for pos = [2 3 6 10 11 12]
    m = m + 1;
    load([backgroundfile int2str(pos) '_re_exp_analysis']);
    load([backgroundfile int2str(pos) '_re_exp']);
    t = 0;
    length = zeros(1,1);
    nucfar1s = zeros(1,1);
    far1conc = zeros(1,1);
    nucfar1W5 = zeros(1,1);
    for i = 1:size(data,1)
        whi5_entry = round(data(i,5));
        whi5_exit = round(data(i,6));
        cellno = interesting_cells(i);
        for j = whi5_entry:whi5_exit
            t = t + 1;
            length(t) = all_obj.volume(cellno,j)/all_obj.area(cellno,j);
            nucfar1s(t) = all_obj.nuc_Far1(cellno,j);
            nucfar1W5(t) = all_obj.total_nuclear_Far1_W5(cellno,j)/all_obj.nuclear_area_W5(cellno,j);
            far1conc(t) = all_obj.Far1_appr_conc(cellno,j);
            
        end
    end
    noNaNindex = ~(isnan(length) | isnan(nucfar1W5)) ;
    p = polyfit(length(noNaNindex), nucfar1W5(noNaNindex), 1);
    p1s(m) = p(1);
    p2s(m) = p(2);
end





k = 0;
perdiff = zeros(1,1,1);
for posanalysis = [1 2 3]
    k = k + 1;
    
    load([analysisfile int2str(posanalysis) '_re_exp_analysis']);
    load([analysisfile int2str(posanalysis) '_re_exp']);

    j = 0;
    for i = 1:size(data,1)
        far1_st = round(data(i,7));
        far1_end = round(data(i,8));
        cellno = interesting_cells(i);
        if far1_end < data_end_point
            j = j +1;
            for l = 1:m
                p1 = p1s(l); p2 = p2s(l);
                zdim = all_obj.volume(cellno,:)./all_obj.area(cellno,:);
                autofl_nucfar1 = p1.*zdim + p2;
                nucfar1 = all_obj.nuc_Far1(cellno,:) - autofl_nucfar1;
                diff = nucfar1(far1_st) - nucfar1(far1_end);
                perdiff(j,l,k) = diff / nucfar1(far1_st);
            end
        end
    end
end


