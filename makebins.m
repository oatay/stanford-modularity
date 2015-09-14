function [bin binmeds binstds] = makebins(x,y,minimum,maximum,binNo)

bins= linspace(minimum, maximum, binNo);

[h,specBins] = histc(x, bins);
binmeds = zeros(1,1);
binstds = zeros(1,1);
bin = zeros(1,1);

ind = 1;
for i = 1:binNo
    binvalues    = y(specBins == i);
    if (length(binvalues) > 3 && mean(binvalues) > 0)
        binmeds(ind)     = mean(binvalues);        
        binstds(ind)     = std(binvalues)/sqrt(length(binvalues));
        bin(ind) = bins(i);
        ind = ind+1;
    end
end
