clear all; clc; close all;
backgroundfile = 'OA_070915_RV200_refurb_3nM6min_pos_no_';
fnorm = 0.25; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');

ind = 1;
indd = 1;
indm = 1;
dalignedatentry = cell(1);
talignedatentry = cell(1);
for pos = [1 2 3 4 5 6 10]
    load([backgroundfile int2str(pos) '_re_exp_Sic1drop']);
    load([backgroundfile int2str(pos) '_re_exp']);
    
    totsic1 = all_obj.nuc_Far1.*all_obj.tot_nucl_areaR+all_obj.cyt_Far1.*all_obj.tot_cyt_areaR;
    
    for i = 1:size(data,1)
        close all
        hold all
        
        startp = round(data(i,6));
        endp = ceil(data(i,7));
        cellno = data(i,1);
        %         totsic1cell = totsic1(cellno,:);
        %         sic1conc = totsic1cell./all_obj.volume(cellno,:);
        %         sic1conc = sic1conc - min(sic1conc(endp:end));
        curr_plot_Sic1nuc = all_obj.nuc_Far1(cellno,:);
        p1 = 7.493; p2 = 47.27;
        zdim = all_obj.volume(cellno,:)./all_obj.area(cellno,:);
        autofl_nucsic1 = p1.*zdim + p2;
        sic1nuc = curr_plot_Sic1nuc-autofl_nucsic1;
        sic1nuc = 0.1.*sic1nuc./zdim;
        filteredsig = filtfilt(b,a,sic1nuc);
        plot(sic1nuc)
        plot(filteredsig)
        plot([startp startp],[0 max(sic1nuc)])
        plot([endp endp],[0 max(sic1nuc)])

        if data(i,3) == 0
            %             plot(sic1conc)
            %             plot([startp startp],[0 max(sic1conc)])
            %             plot([endp endp],[0 max(sic1conc)])
            indm = indm+1;
        else
            %             plot(sic1conc)
            %             plot([startp startp],[0 max(sic1conc)])
            %             plot([endp endp],[0 max(sic1conc)])
            [x1,y1,button] = ginput(1);
            x1 = round(x1);
            dalignedatentry{indd} = sic1nuc(x1:endp);
            tdalignedatentry{indd} = ((x1:endp) - endp).*6 + 6;

            indd = indd+1;
        end
        
        ind = ind + 1;
        
    end
end

if indd > 10
    save(['sic1alignment' '_daughters'], 'dalignedatentry','tdalignedatentry')
end
%%
clear all;clc;close all
load 'sic1alignment_daughters'
hold all
td = zeros(1,1);
yd = zeros(1,1);
for i = 1:16
    dalignedatentry{i} = dalignedatentry{i}- nanmin(dalignedatentry{i});
    plot(tdalignedatentry{i},dalignedatentry{i})
    curx = tdalignedatentry{i}(1:end);
    cury = dalignedatentry{i}(1:end);
    td = [td curx];
    yd = [yd cury];
end
td = td(2:end);
yd = yd(2:end);


[bin binmeds binstds] = makebins(td,yd,-200,20,30);
figure(2)
hold all
ciplot((binmeds-binstds), (binmeds+binstds), bin,'r')
plot(bin,binmeds,'LineWidth',3)
ylim([0 1])
xlim([-200 10])
hold off
