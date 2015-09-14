clear all; clc; close all;
backgroundfile = 'OA_071615_OA049_50GalFar13nM_pos_no_';

fnorm = 0.25; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');

indd = 1;
indm = 1;
hills = zeros(1,1);
for pos = [4:7 10]
    load([backgroundfile int2str(pos) '_re_exp_Far1drop']);
    load([backgroundfile int2str(pos) '_re_exp']);
    
    for i = 1:size(data,1)
        close all
        hold all
        
        startp = round(data(i,5));
        endp = round(data(i,6));
        backgp = round(data(i,7));
        cellno = data(i,1);
        
        autofluorescence = all_obj.Far1_appr_conc(cellno,backgp);

        far1sig = all_obj.Far1_appr_conc(cellno,startp:endp) - autofluorescence;
        far1conc =  all_obj.Far1_appr_conc(cellno,:)-autofluorescence;
        time = (startp:endp).*6;
        
        expft=fittype('exp1');
        cf=fit(time',far1sig',expft);
        coeffab = coeffvalues(cf);
        decaya = coeffab(1);
        decayb = coeffab(2);
        figure(1)
        plot(time,far1sig,'b')
        plot(time,decaya.*exp(decayb.*time),'r')
        figure(2)

        1/decayb
        checkbutton=1;
        while checkbutton ~= 2;
            [x1,y1,button] = ginput(1);
            if button == 1
                if data(i,3) == 1 % daughter
                    chosencellposd(indd) = pos;
                    chosencellnod(indd) = cellno;
                    decaydaughters(indd) = decayb;
                    indd = indd + 1;
                    checkbutton = 2;
                else % mother
                    chosencellposm(indm) = pos;
                    chosencellnom(indm) = cellno;
                    decaymothers(indm) = decayb;
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

decaytimemothers = 1./decaymothers;
decaytimedaughters = 1./decaydaughters;
if size(decaytimedaughters,2) > 10
    save([backgroundfile,'decaytimes'],'decaytimemothers','decaytimedaughters',...
        'chosencellposd','chosencellnod');
end