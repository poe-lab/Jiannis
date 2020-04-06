function detect_ripples(session,refs)

%--------------------------------------------------------------------------
% Detects ripples based on the RMS of the "Fripple" bandpassed LFP from a
% reference shank/channel, when it goes above a threshold.
% Ripple edges are detected based on a second threshold. Only segments
% corresponding to rat immobility are considered.
% The ripple beginning/ending/peak times are stored as columns in a matrix.
% The ripples are aligned over their peak times (for plotting only).
%
% Created by Jiannis Taxidis, Caltech, USA, February 2013
%--------------------------------------------------------------------------

refshank = refs(1);                                                         % The first element of "refs" is the reference shank
refchannel = refs(2);                                                       % the second one is the reference channel

Fs = 1017.35;
Fripple = [120,200];                                                        % Ripple frequency range to bandpass LFP

bindur = 5;                                                              % Bin duration for RMS (4.8 msec) --> Now 4.9 ms
disp(['The bin duration for the rms calculation is: ',num2str(bindur/Fs*1000),' msec.']);
crit1 = 2;     
crit2 = 1;
critV = 10;
periripdur = 0.15;                                                          % Duration around ripple (300 msec total window length (like Kamran's));                                                        
mindur = 0.02;                                                              % Minimum accepted ripple duration (20 msec)
mindist = 0.05;                                                             % Minimum accepted distance between ripples (50 msec)

peririp = ceil(2*periripdur*Fs);                                    

figure;

%% LOAD LFPS AND BANDPASS THEM
%-------
% load your LFP here 
% LFP = ......;                                                    % Keep the LFP from the reference channel
%----------

LFPbp = bandpass_c(LFP, Fs, Fripple(1), Fripple(2), 'pass');                                % Filter the LFP

N = length(LFP);
time = (0:N-1)/Fs;
bins = 1 : bindur/2 : (N-bindur);                                           % Array of locations of bin starting times
time_b = time(bins + bindur/2);

%% PLOT POWER SPECTRUM
seg = 10;
[f,~,~] = sp2a2_m1(0,LFP,LFP,Fs,seg,'r2 t2 h');

subplot(4,5,[11,16])
hold on
plot(f(:,1), (f(:,2)+log10(2*pi)),'Linewidth',1.5,'Color','k')
xlabel('Frequency (Hz)','Fontsize',20)
xlim([0 500])
set(gca,'Fontsize',16)

%% CALCULATE ROOT-MEAN-SQUARE OVER 50% OVERLAPPING BINS
RMS = zeros(length(bins),1);
for i = 1:length(bins)
    RMS(i) = sum(LFPbp(bins(i):(bins(i)+bindur)).^2);                       % Calculate rms for each bin (IT IS THE SUM OF SQUARES ONLY!!! units = mV^2)
end


%% DETECT ALL SEGMENTS ABOVE THE LFP THRESHOLD AND THEIR LIMITS
stdr = std(RMS);                                                            % Standard deviation of rms
CRIT1 = crit1*stdr;                                                         % Set the detection criterion
CRIT2 = crit2*stdr;                                                         % Set the limits detection criterion
RMS(1) = 0;                                                                 % Set the first and last points to zero
RMS(end) = 0;                                                               % to avoid detecting on the edges.

cross = (RMS >= CRIT1);                                                     % Find the array locations over the threshold
crossbd = find((cross - circshift(cross,1)) == 1);                          % Find locations whose preceding location is below the threshold (ripple beginning)
crossbd = [crossbd find(cross - circshift(cross,1) == -1)];                 % Add locations whose following location is below the threshold (ripple end)

ripples = zeros(length(crossbd(:,1)),3);
for i = 1:length(crossbd(:,1))                                              % For each ripple beginnning
    rtmp = RMS(crossbd(i,1):crossbd(i,2));                                  % Keep the rms
    rippeak = find(rtmp == max(rtmp)) + crossbd(i,1);                       % Find the overall location of the rms peak
    ripples(i,1) = find(RMS(1:rippeak) < CRIT2,1,'last');                   % Keep the last location before ripple beginnning below the limit threshold
    ripples(i,2) = find(RMS(rippeak:end) < CRIT2,1,'first') + rippeak - 1;  % and the first location after ripple end below the limit threshold
    ripples(i,3) = rippeak;                                                 % and the location of the ripple peak
end
ripples = time_b(ripples);

[~,m] = unique(ripples(:,1), 'rows');                                       % Keep unique events by checking repeating ripple beginning times
ripples = ripples(m,:);
[~,m] = unique(ripples(:,2), 'rows');                                       % Keep unique events by checking repeating ripple end times
ripples = ripples(m,:);

disp('Initially :')
[~,rdur] = ripple_stats(ripples);

%% REMOVE SHORT AND OVERLAPPING RIPPLES
rshort = (rdur < mindur);                                                   % If any ripple lasts less than 20 msec discard it
ripples(rshort,:) = [];

k = ripples(2:end,1) - ripples(1:(end-1),2);                                % Find the distances between all ripples
k = (k < mindist);                                                          % Find when the distance is not at least 50 msec (covers overlaps too)
k1 = [k; 0];
k2 = [0; k];                                                                % Add a zero for the first ripple that was not included
ripples(k1==1,2) = ripples(k2==1,2);                                        % Unite them (use the ending of the second ripple as end of the first one...
ripples(k2==1,:) = [];                                                      % ...and delete the times of the second one)

disp('After close-ripples removal: ')
[rN,~] = ripple_stats(ripples);

%% REMOVE THETA-TIME RIPPLES
% load detect_movement data here
load(['../Sessions/',session,'/Movement_data.mat']);                                 % Load the movement analysis of this session

movein = [];
moveout = [];
for r = 1:rN
    k = find(move_in(:,1) < ripples(r,3),1,'last');                         % Find the last incoming motion segment beginning before the ripple
    if move_in(k,2) > ripples(r,3)                                          % If the segment's end is after the ripple then the ripple is inside the segment
        movein = [movein, r];                                               % Store it
    end
    k = find(move_out(:,1) < ripples(r,3),1,'last');                        % Same for outgoing segments
    if move_out(k,2) > ripples(r,3)
        moveout = [moveout, r];
    end
end
move = unique([movein,moveout]);                                            % Unite all motion-ripples
ripples(move,:) = [];                                                       % and remove them

immob = (V <= critV);                                                       % Find where the velocity is below the threshold
move = [];
for r = 1:size(ripples,1)
    z = find(timeV <= ripples(r,3),1,'last');                               % Find the last velocity-timepoint before the ripple
    t1 = immob(z);                                                          % Check if it is immobile
    z = find(timeV >= ripples(r,3),1,'first');                              % Same for the last timepoint after the ripple
    t2 = immob(z);
    if t1 == 0 && t2 == 0                                                   % If both are not immobile timepoints
        move = [move, r];                                                   % Keep the ripple's number
    end
end
ripples(move,:) = [];                                                       % Remove non-immobile ripples

disp('After running-ripples removal: ')
[rN,~] = ripple_stats(ripples);


%% PLOT LFPS AND RIPPLE TIMES
subplot(4,5,1:10)
hold on;
plot(timeV, V/100+abs(min(LFP))+max(LFPbp)+max(LFP),'k')
line([time(1),time(end)], [critV,critV]/100+abs(min(LFP))+max(LFPbp)+max(LFP),'Color','k')
plot(time, LFP+abs(min(LFP))+max(LFPbp),'b');
plot(time, LFPbp,'b');
plot(time_b, RMS,'k','LineWidth',1.5);


line([time(1),time(end)], [CRIT1,CRIT1],'Color','k')
line([time(1),time(end)], [CRIT2,CRIT2],'Color','k')
z = max(V/100)+abs(min(LFP))+max(LFPbp)+max(LFP);
for i = 1:rN
    line([ripples(i,1),ripples(i,1)],[min(LFPbp) z],'Color','b');
    line([ripples(i,2),ripples(i,2)],[min(LFPbp) z],'Color','r');
end

%% ALIGN RIPPLES
rips = round(ripples*Fs);
LFPrip = zeros(rN,peririp+1);
LFPbprip = zeros(rN,peririp+1);

rmax = zeros(1,rN);
rout = [];
for r = 1:rN                                                                % For each event
    lfptmp = LFPbp(rips(r,1):rips(r,2));                                    % Keep the corresponding bandpassed signal
    [~,rm] = max(lfptmp);                                                   % Find the position of its maximum value
    rmax(r) = rips(r,1) - 1 + rm;                                           % Store the overall array position of the peak
    if rmax(r)-peririp/2 <= 0 || rmax(r)+peririp/2 > N                      % If the perievent duration hits the LFP boundaries
        rout = [rout r];                                                    % Note the number of the event
    end
end

ripples(rout,:) = [];
rmax(rout) = [];
rN = size(ripples,1);

for r = 1:rN
    LFPrip(r,:) = LFP(rmax(r) + (-peririp/2 : peririp/2));
    LFPbprip(r,:) = LFPbp(rmax(r) + (-peririp/2 : peririp/2));
end

%% PLOT ALIGNED RIPPLES
for r = 1:rN
    subplot(4,5,12:15)
    hold on
    plot((-peririp/2:peririp/2)/Fs,LFPrip(r,:),'Color',[0.5 0.5 0.5])
    subplot(4,5,17:20)
    hold on
    plot((-peririp/2:peririp/2)/Fs,LFPbprip(r,:),'Color',[0.5 0.5 0.5])
end
subplot(4,5,12:15)
plot((-peririp/2:peririp/2)/Fs,mean(LFPrip),'LineWidth',3,'Color','k')
axis tight
subplot(4,5,17:20)
plot((-peririp/2:peririp/2)/Fs,mean(LFPbprip),'LineWidth',3,'Color','k')
axis tight
ylabel({'Average';'Ripple'},'Fontsize',20)
xlabel('Lags (sec)','Fontsize',20)

%% PLOT SPECTRUM OF RIPPLE SECTIONS
LFPtemp = LFPrip';
LFPtemp = LFPtemp(:);
[f,~,~] = sp2a2_m1(0,LFPtemp,LFPtemp,Fs,seg,'r2 t2 h');
subplot(4,5,[11,16])
plot(f(:,1), (f(:,2)+log10(2*pi)),'Linewidth',1.5,'Color','b')

%% SAVE DATA AND FIGURE
save(['Ripple_data',num2str(refshank),'_',num2str(refchannel),'.mat'],'ripples','LFPbp','time','time_b','LFPrip','LFPbprip','rmax','periripdur');

%--------------------------------------------------------------------------
%% ------------------------------------------------------------------------

function [rN,rdur] = ripple_stats(ripples)

rN = size(ripples,1);                                                       % Number of ripples in this session
rdur = ripples(:,2)-ripples(:,1);                                           % and ripple durations
rmax = max(rdur);                                                           % Maximum duration
disp([num2str(rN),' ripples of average duration ', num2str(mean(rdur)*1000),' msec and maximum duration ',num2str(rmax*1000),' msec',' (SD = ',num2str(std(rdur)*1000),')'])
