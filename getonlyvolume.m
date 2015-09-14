clear all; clc; close all;
filename = 'AD_wtc5p_30_17c_oct31_2013_pos_no_';

for pos = [1:5 7:19]
    load([filename int2str(pos) '_2']);
    no_obj = size(all_obj.nuc_whi5R,1);
    pos
    all_obj.volume              = zeros(no_obj,numbM);
    all_obj.area                = zeros(no_obj,numbM);
    all_obj.proj_area           = zeros(no_obj,numbM);%projected area
    all_obj.majAp_ax            = zeros(no_obj,numbM);
    all_obj.minAp_ax_med        = zeros(no_obj,numbM);
    all_obj.minAp_ax_mean       = zeros(no_obj,numbM);
    all_obj.minAp_ax_sd         = zeros(no_obj,numbM);
    all_obj.minAp_ax_len        = zeros(no_obj,numbM);

    for i = 1:no_obj
        for c_time = 1:170
            new_c_Image2 = (all_obj.cells(:,:,c_time) == i);
            [vol_c,area_c,maj_appr_axis,min_appr_axis]=get_area2(new_c_Image2,c_time,i); %i=cell number
            all_obj.volume(i,c_time)    = vol_c;
            all_obj.area(i,c_time)      =  area_c;
            all_obj.proj_area(i,c_time) = sum(new_c_Image2(:)>0);
            all_obj.majAp_ax(i,c_time)= maj_appr_axis;
            all_obj.minAp_ax_med(i,c_time)= median(min_appr_axis(min_appr_axis>0));
            all_obj.minAp_ax_mean(i,c_time)= mean(min_appr_axis(min_appr_axis>0));
            all_obj.minAp_ax_sd(i,c_time)  =sqrt(var(min_appr_axis(min_appr_axis>0)));
            all_obj.minAp_ax_len(i,c_time)  =length(min_appr_axis(min_appr_axis>0));
        end
    end
    name1=[filename num2str(pos) '_vol'];
    save(name1, 'all_obj', 'prefix', 'suffix2', 'pos_num', 'suffix', 'numbM', 'type', 'LcellsO' ...
        ,'cell_exists','max_size_vs_largest_cell','max_area_increase_per_tp','higher_threshold','lower_threshold' ...
        ,'threshold_increase_factor','phase_subtraction_factor','GFP_subtraction_factor');
end
