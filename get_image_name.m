function [im_name] = get_image_name(prefix,pos_num,suffix,numb,...
    c_num,type,suffix2,numbM)

% prefix  ='AD_22apr09';
% pos_num =12;
% suffix  ='_gfp';
% numb    =10;
% c_num   =c4;
% type    ='.jpg'

if numbM<70 + 30
    if numb<10
        rnumb=['00' num2str(numb)];%00
    elseif numb<100
        rnumb=['0' num2str(numb)];%0
    else
        rnumb=num2str(numb);
    end
elseif numbM>69 + 30
    if numb<10
        rnumb=['00' num2str(numb)];%00
    elseif numb<100
        rnumb=['0' num2str(numb)];%0
    else
        rnumb=num2str(numb);
    end
elseif numbM
    disp('strange number given')
else
end

if length(suffix2)==1
im_name=[prefix '-000' suffix2 '_Position(' num2str(pos_num)...
    ')' suffix 't' rnumb num2str(c_num) type ];
else
    im_name=[prefix '-00' suffix2 '_Position(' num2str(pos_num)...
    ')' suffix 't' rnumb num2str(c_num) type ];
end
%AD_22apr09            -0003_Position(12)_gfp_t001c4.JPG
%AD_growth_test_290409 -0003_Position(10)_t001(c1)

%imread(im_name);