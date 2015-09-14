clear all; clc; close all;
backgroundfile = 'OA_070915_RV200_refurb_3nM6min_pos_no_';
fnorm = 0.25; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');

ind = 1;
indd = 1;
indm = 1;

for pos = 3% [1 2 3 4 5 6 10]
    load([backgroundfile int2str(pos) '_re_exp_Sic1drop']);
    load([backgroundfile int2str(pos) '_re_exp']);
    
    for i = 1:size(data,1)
        close all
        hold all
        startp = data(i,6);
        endp = data(i,7);
        
        
        cellno = data(i,1);
        
        from = floor(startp);
        to = ceil(endp);       
        
        curr_plot_Sic1nuc = all_obj.nuc_Far1(cellno,:);
        p1 = 7.493; p2 = 47.27;
        zdim = all_obj.volume(cellno,:)./all_obj.area(cellno,:);
        autofl_nucsic1 = p1.*zdim + p2;
        sic1nuc = curr_plot_Sic1nuc-autofl_nucsic1;
        if cellno == 13
            sic1nuc = 0.1.*sic1nuc./zdim;
            plot(sic1nuc)
            plot([startp startp],[0 max(sic1nuc)])
            plot([endp endp],[0 max(sic1nuc)])
            ylim([-0.2 1.2*max(sic1nuc)])
            pos
            cellno
            pause(15)
        end
        
        sic1nuc = sic1nuc(from:to);

       
        
        midsic1= (max(sic1nuc) + min(sic1nuc))/2;
        if sic1nuc(2) > sic1nuc(1)
            sic1nuc = sic1nuc(2:end);
            from = from + 1;
        end
        minsic1nucvalues(ind) = min(sic1nuc);
        maxsic1nucvalues(ind) = max(sic1nuc);
        midtimepoint = interp1(sic1nuc,from:to,midsic1);
        (midtimepoint - startp)*6
        
        if data(i,3) == 0
            mtimetohalfmax(indm) = (midtimepoint - startp)*6;
            mbeforedecaytimes(indm) = (data(i,6)-data(i,5))*6;
            indm = indm+1;
        else
            dtimetohalfmax(indd) = (midtimepoint - startp)*6;
            dbeforedecaytimes(indd) = (data(i,6)-data(i,5))*6;
            indd = indd+1;
        end
        ind = ind + 1;

    end
end

median(mtimetohalfmax)
std(bootstrp(100,@median,mtimetohalfmax))

median(dtimetohalfmax)
std(bootstrp(100,@median,dtimetohalfmax))

timetohalfmax = [mtimetohalfmax dtimetohalfmax];
median(timetohalfmax)
std(bootstrp(100,@median,timetohalfmax))
