function [rN,rdur] = ripple_stats(ripples)

rN = size(ripples,1);                                                       % Number of ripples in this session
rdur = ripples(:,2)-ripples(:,1);                                           % and ripple durations
rmax = max(rdur);                                                           % Maximum duration
disp([num2str(rN),' ripples of average duration ', num2str(mean(rdur)*1000),' msec and maximum duration ',num2str(rmax*1000),' msec',' (SD = ',num2str(std(rdur)*1000),')'])
end