clear all; clc; close all;
backgroundfile = 'OA_061015_OA040_3nM6min_pos_no_';
fnorm = 0.25; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');
indd = 1;
indm = 1;
ind = 1;
inda = 1;
for pos = [2 3 4 6 7 10 12 13]
    load([backgroundfile int2str(pos) '_re_exp_Far1analysis']);
    load([backgroundfile int2str(pos) '_re_exp']);
    
    for i = 1:size(data,1)
        close all
        hold all
        startp = data(i,5);
        endp = data(i,6);
        totaldecaytimes(ind) = (endp-startp)*6;
        cellno = data(i,1);

        from = floor(startp);
        to = ceil(endp);       
        
        curr_plot_Far1nuc = all_obj.nuc_Far1(cellno,:);
        p1 = 7.493; p2 = 47.27;
        zdim = all_obj.volume(cellno,:)./all_obj.area(cellno,:);
        autofl_nucfar1 = p1.*zdim + p2;
        far1nuc = curr_plot_Far1nuc-autofl_nucfar1;
        far1nucnorm = far1nuc./zdim;
        
        plot(far1nucnorm)
        plot([startp startp],[0 max(far1nucnorm)])
        plot([endp endp],[0 max(far1nucnorm)])
        
        whi5nuc = all_obj.nuc_whi5R(cellno,:)-all_obj.cyt_whi5R(cellno,:);
        whi5nuc = 1.2.*max(whi5nuc)/max(far1nuc).*far1nuc;
        %plot(whi5nuc)
        pos
        cellno
        pause(5)
%         checkbutton=1;
%         while checkbutton ~= 2;
%             
%             [x1,y1,button] = ginput(1);
%             plot([x1 x1],[0 max(far1nuc)],'color',[1 1 0],'linewidth',2)
% 
%             if button == 1 && data(i,3) == 1 %daughterarresttimes
%                 display('choose whi5in')
%                 x1 = round(x1);
%                 totalarrestd(inda) = (x1-data(i,5))*6
%                 inda = inda + 1;
%                 checkbutton = 2;
%             else
%                 checkbutton=2;
%             end
%         
%         end
        
        far1nuc = far1nuc(from:to);
%         pos
%         cellno
%         pause(2)
        
        
        midfar1= (max(far1nuc) + min(far1nuc))/2;
        if far1nuc(2) > far1nuc(1)
            far1nuc = far1nuc(2:end);
            from = from + 1;
        end
        if data(i,3) == 1
            minfar1nucvaluesd(indd) = min(far1nuc);
            maxfar1nucvaluesd(indd) = max(far1nuc);
            midtimepoint = interp1(far1nuc,from:to,midfar1);
            (midtimepoint - startp)*6
            timetohalfmaxd(indd) = (midtimepoint - startp)*6;
            indd = indd + 1;
        else
            minfar1nucvaluesm(indm) = min(far1nuc);
            maxfar1nucvaluesm(indm) = max(far1nuc);
            midtimepoint = interp1(far1nuc,from:to,midfar1);
            (midtimepoint - startp)*6
            timetohalfmaxm(indm) = (midtimepoint - startp)*6;
            indm = indm + 1;
        end

    end
end
median(timetohalfmaxm)
std(bootstrp(100,@median,timetohalfmaxm))

dt = median(timetohalfmaxd)
dstd = std(bootstrp(100,@median,timetohalfmaxd))

timetohalfmax = [timetohalfmaxd timetohalfmaxm];
median(timetohalfmax)
std(bootstrp(100,@median,timetohalfmax))


barwitherr([ 34.36 dstd], [201 dt])

%%
close all; clear all; clc;

pos = 13;%6;%13;%3;% 6;
cellno = 8;%1;%8;%42;%1;
i = 1;
backgroundfile = 'OA_061015_OA040_3nM6min_pos_no_';
load([backgroundfile int2str(pos) '_re_exp_Far1analysis']);
load([backgroundfile int2str(pos) '_re_exp']);
area = all_obj.area;

background1 = median(all_obj.med_backgr_w5)/29.*(30.7743+area(area < 3215).*0.0010);%%%background subt.
j2 = logical((area >= 3215).*(area <4180));
background2 = median(all_obj.med_backgr_w5)/29.*(35.0318-area(j2).*0.0004);
background3 = median(all_obj.med_backgr_w5)/29.*(30.8450+area(area >= 4180).*0.0006);

all_obj.nuc_whi5R(area < 3215) = all_obj.nuc_whi5R(area < 3215) - background1;
all_obj.nuc_whi5R(j2) = all_obj.nuc_whi5R(j2)-background2; %%%%%%%%%%%%
all_obj.nuc_whi5R(area >= 4180) = all_obj.nuc_whi5R(area >= 4180)-background3; %%%%%%%%%%%%
all_obj.cyt_whi5R(area < 3215) = all_obj.cyt_whi5R(area < 3215) - background1;
all_obj.cyt_whi5R(j2) = all_obj.cyt_whi5R(j2)-background2; %%%%%%%%%%%%
all_obj.cyt_whi5R(area >= 4180) = all_obj.cyt_whi5R(area >= 4180)-background3;


curr_plot_Far1nuc = all_obj.nuc_Far1(cellno,:);
p1 = 7.493; p2 = 47.27;
zdim = all_obj.volume(cellno,:)./all_obj.area(cellno,:);
autofl_nucfar1 = p1.*zdim + p2;
far1nuc = curr_plot_Far1nuc-autofl_nucfar1;
far1nuc(isnan(far1nuc)) = 0;
whi5nuc = 0.5.*(all_obj.nuc_whi5R(cellno,:)-all_obj.cyt_whi5R(cellno,:))./zdim;


fnorm = 0.80; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');

totnucFar1 = all_obj.total_nuclear_Far1_F1(cellno,:);%-autofl_nucfar1.*all_obj.nuclear_area_F1(cellno,:);
far1nucnorm = far1nuc./zdim;

figure(1)
hold on
xlim([0 139*6])
plot(0:6:139*6, far1nucnorm(1:end),'r','LineWidth',2)
plot(0:6:139*6, whi5nuc(1:end),'b','LineWidth',2)
plot([120 120],[0 max(far1nucnorm)/2])
hold off
figure(2)
plot(228:6:139*6, filtfilt(b,a,far1nuc(39:end)))
% plot([startp startp],[0 max(far1nuc)])
% plot([endp endp],[0 max(far1nuc)])
%%
filename = 'OA_061015_OA040_3nM6min_pos_no_13_re_exp';
load(filename)
cellno = 8; 
%pos=2; cellno =1; t 25;50;75;
initial_tp=56; % t=113-116
c_num   ='c1';%phase
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);
I=imread(im_name);
c_num   ='c3';
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);
IN=imread(im_name);

c_num   ='c2';
[im_name] = get_image_name(prefix,pos_num,suffix,initial_tp, c_num,type,suffix2,numbM);
IG=imread(im_name);

outline=(255.*uint8(bwmorph(all_obj.cells(:,:,initial_tp)==cellno,'remove')));
figure(1)
imshow(I+outline);title('phase')
figure(2)
imshow(outline+IN);title('mKok')
figure(3)
imshow(outline+IG);title('GFP')
%%

close all;
pos = 6;%3;%6;%13;%3;% 6;
cellno = 1;%42;%1;%8;%42;%1;
i = 1;
backgroundfile = 'OA_061015_OA040_3nM6min_pos_no_';
load([backgroundfile int2str(pos) '_re_exp_Far1analysis']);
load([backgroundfile int2str(pos) '_re_exp']);
area = all_obj.area;

background1 = median(all_obj.med_backgr_w5)/29.*(30.7743+area(area < 3215).*0.0010);%%%background subt.
j2 = logical((area >= 3215).*(area <4180));
background2 = median(all_obj.med_backgr_w5)/29.*(35.0318-area(j2).*0.0004);
background3 = median(all_obj.med_backgr_w5)/29.*(30.8450+area(area >= 4180).*0.0006);

all_obj.nuc_whi5R(area < 3215) = all_obj.nuc_whi5R(area < 3215) - background1;
all_obj.nuc_whi5R(j2) = all_obj.nuc_whi5R(j2)-background2; %%%%%%%%%%%%
all_obj.nuc_whi5R(area >= 4180) = all_obj.nuc_whi5R(area >= 4180)-background3; %%%%%%%%%%%%
all_obj.cyt_whi5R(area < 3215) = all_obj.cyt_whi5R(area < 3215) - background1;
all_obj.cyt_whi5R(j2) = all_obj.cyt_whi5R(j2)-background2; %%%%%%%%%%%%
all_obj.cyt_whi5R(area >= 4180) = all_obj.cyt_whi5R(area >= 4180)-background3;


curr_plot_Far1nuc = all_obj.nuc_Far1(cellno,:);
p1 = 7.493; p2 = 47.27;
zdim = all_obj.volume(cellno,:)./all_obj.area(cellno,:);
autofl_nucfar1 = p1.*zdim + p2;
far1nuc = curr_plot_Far1nuc-autofl_nucfar1;
far1nuc(isnan(far1nuc)) = 0;
whi5nuc = 0.5.*(all_obj.nuc_whi5R(cellno,:)-all_obj.cyt_whi5R(cellno,:))./zdim;


fnorm = 0.80; % (arbitrary) decrease to increase smoothing
[b,a] = butter(3, fnorm, 'low');

totnucFar1 = all_obj.total_nuclear_Far1_F1(cellno,:);%-autofl_nucfar1.*all_obj.nuclear_area_F1(cellno,:);
far1nucnorm = far1nuc./zdim;

figure(1)
hold on
xlim([65*6 139*6])
plot(70*6:6:139*6, far1nucnorm(71:140),'r','LineWidth',2)
%plot(70*6:6:139*6,whi5nuc)
hold off
%figure(2)
%plot(228:6:139*6, filtfilt(b,a,far1nuc(39:end)))
% plot([startp startp],[0 max(far1nuc)])
% plot([endp endp],[0 max(far1nuc)])