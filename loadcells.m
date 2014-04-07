clear all
positions = [8 11:19 21];
filename = 'OA_030714_OA039L381-L392-Longpulse_pos_no_5_re_exp';
load(filename)
lim_h=[11 41 46 56];   % [31 61 66 76]         % <- change for each new exp
% [21 26]
% [31 61 66 76]
% [21 41 42 46 47 51 52 56]
% [11 41 46 56]
% [21 41 42 46 47 51 52 56 57 61]
% [21 41 44 51 54 61 64 71]
no_cells=length(all_obj.max_nucl_int(:,1));
interesting_cells = zeros(0,0);
add_data=1;
if add_data==1
    disp(['no of cells here = ' num2str(no_cells)])
    i4=1; 
    for no_cell_of_interest= 1:no_cells
        initial_tp=lim_h(2);
        intr_cell= monitor_time_series_OA_1('c3',initial_tp,no_cell_of_interest,prefix,pos_num,suffix,type,suffix2,numbM,all_obj,lim_h);
        if sum(intr_cell) > 0
            interesting_cells(i4)=no_cell_of_interest;
            data(i4,:) = intr_cell;
            i4=i4+1;
        end
    end
    if ~isempty(interesting_cells)
        save([filename '_' 'analysis'],'all_obj', 'data', 'interesting_cells','lim_h')
    end
    
    
    %----------------------------------------------------------
    %[cell_values] = get_f1w5_data_1(interesting_cells,all_obj,lim_h,numbM);%<-this is for 240 steps
    %   [cell_values] = get_f1w5_data_2(interesting_cells,all_obj,lim_h,numbM);%<- this is for cycling cells
    %-   save([prefix '_' num2str(pos_num) '_est_half_times_ENDV_1'], 'cell_values','lim_h');
    %   save([prefix '_' num2str(pos_num) '_est_half_times_clb56n_1'], 'cell_values','lim_h');
else
end