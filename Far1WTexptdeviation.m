clear all; clc; close all;
backgroundfile = 'OA_061015_OA040_3nM6min_pos_no_';

fnorm = 0.25; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');

indd = 1;
indm = 1;
perdegrd = zeros(1,1);
perdegrm = zeros(1,1);


for pos = [2 3 4 6 7 10 12 13]
    load([backgroundfile int2str(pos) '_re_exp_Far1analysis']);
    load([backgroundfile int2str(pos) '_re_exp']);
    
    for i = 1:size(data,1)
        close all
        hold all
        
        cellno = data(i,1);
        

        curr_plot_Far1nuc = all_obj.nuc_Far1(cellno,:);
        p1 = 7.493; p2 = 47.27;
        zdim = all_obj.volume(cellno,:)./all_obj.area(cellno,:);
        autofl_nucfar1 = p1.*zdim + p2;
        far1nuc = curr_plot_Far1nuc-autofl_nucfar1;
        far1nuc(isnan(far1nuc)) = 0;
        whi5nuc = (all_obj.nuc_whi5R(cellno,:)-all_obj.cyt_whi5R(cellno,:))./zdim;

        totnucFar1 = all_obj.total_nuclear_Far1_F1(cellno,:);%-autofl_nucfar1.*all_obj.nuclear_area_F1(cellno,:);
        far1nucnorm = far1nuc./zdim;
        plot(far1nucnorm)
        plot(whi5nuc,'r')
        plot([data(i,5) data(i,5)], [0 max(far1nucnorm)])
        
        checkbutton=1;
        while checkbutton ~= 2;
            [x1,y1,button] = ginput(1);
            startp = round(x1);
            plot([startp startp], [0 max(far1nucnorm)])
            
            [x2,y2,button] = ginput(1);
            endp = round(x2);
            plot([endp endp], [0 max(far1nucnorm)])
            
            [x3,y3,button] = ginput(1);
            backgr = round(x3);
            plot([endp endp], [0 max(far1nucnorm)])
            
            far1Catarrest = far1nucnorm(startp)-far1nucnorm(backgr);
            far1Cattheend = far1nucnorm(endp)-far1nucnorm(backgr);

            
            perdegr = (far1Catarrest - far1Cattheend)/far1Catarrest;
            if perdegr < 0
                perdegr = 0;
            end
            perdegr
            [x4,y4,button] = ginput(1);
            if button == 1
                if data(i,3) == 1 % daughter
                    chosencellposd(indd) = pos;
                    chosencellnod(indd) = cellno;

                    perdegrd(indd) = perdegr;
                    indd = indd + 1;
                    checkbutton = 2;
                else % mother
                    chosencellposm(indm) = pos;
                    chosencellnom(indm) = cellno;
                    
                    perdegrm(indm) = perdegr;

                    indm = indm + 1;
                    checkbutton = 2;
                end
            elseif button == 3
                checkbutton = 2;
            end
            
        end
        hold off
    end
end


if size(perdegrd,2) > 10
    save([backgroundfile,'wtexptdeviation'],'perdegrd','perdegrm',...
        'chosencellposd','chosencellnod');
end