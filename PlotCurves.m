for i = 18
    volume_ax = all_obj.volume_ax(i,:);
    volume_fl = all_obj.volume_fl(i,:);
    volume = all_obj.area(i,:).^1.5;
    close all
    hold on
    plot(volume_ax,'r')
    plot(volume_fl,'g')
    plot(volume,'b')
    pause(1)
    hold off
end