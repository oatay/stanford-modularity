clear all; clc; close all;
backgroundfile = 'OA_071015_OA042_refurb_3nM6min_pos_no_';

fnorm = 0.25; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');

indd = 1;
indm = 1;

allydata =cell(1);
allxdata =cell(1);
allydatam = cell(1);
allxdatam = cell(1);

for pos = [2 3 5 7 8 10 11 12]
    load([backgroundfile int2str(pos) '_re_exp_Far1L92Pnucentry']);
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
%             [x1,y1,button] = ginput(1);
%             startp = round(x1);
            startp = round(data(i,5));
            plot([startp startp], [0 max(far1nucnorm)])
            
%             [x2,y2,button] = ginput(1);
%             endp = round(x2);
            endp = round(data(i,6));
            plot([endp endp], [0 max(far1nucnorm)])
            
%             [x3,y3,button] = ginput(1);
%             backg = round(x3);
            %backg = round(data(i,7)); % not a very good backgr for far1nucnorm
            
          
            [x4,y4,button] = ginput(1);
            backg = round(x4);
            plot([backg backg], [0 max(far1nucnorm)])
            
            if button == 1
                if data(i,3) == 1 % daughter
                    chosencellposd(indd) = pos;
                    chosencellnod(indd) = cellno;
                    
                    ydata = far1nucnorm((startp):endp+4)-far1nucnorm(backg);
                    xdata = startp:(endp+4);
                    allydata{indd} = ydata;
                    allxdata{indd} = xdata;
                    
                    indd = indd + 1;
                    checkbutton = 2;
                else 
                    chosencellposm(indm) = pos;
                    chosencellnom(indm) = cellno;
                    
                    % for mothers take 4 data points before arrest as well
                    % align at just before reentry but take 4 more points
                    ydata = far1nucnorm((startp-4):(endp+4))-far1nucnorm(backg);
                    xdata = (startp-4):(endp+4);
                    allydatam{indm} = ydata;
                    allxdatam{indm} = xdata;
                    
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

if size(allydata,2) > 15
    save([backgroundfile,'alignedFar1L92P'],'allxdata','allydata','allxdatam','allydatam',...
        'chosencellposd','chosencellnod','chosencellposm','chosencellnom');
end

%%
clear all; clc; close all;
load OA_071015_OA042_refurb_3nM6min_pos_no_alignedFar1L92P
hold all
td = zeros(1,1);
yd = zeros(1,1);

for i = 1:size(allydata,2)
    curx = allxdata{i};
    curendp = curx(end);
    allxdata{i} = (allxdata{i} - curendp).*6;
    plot(allxdata{i}, allydata{i})
    curx = allxdata{i}(1:end);
    cury = allydata{i}(1:end);
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