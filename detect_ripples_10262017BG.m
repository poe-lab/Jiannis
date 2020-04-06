function detect_ripples_10262017BG(session,refs)
%--------------------------------------------------------------------------
% Detects ripples based on the RMS of the "Fripple" bandpassed LFP from a
% reference shank/channel, when it goes above a threshold.
% Ripple edges are detected based on a second threshold. Only segments
% corresponding to rat immobility are considered.
% The ripple beginning/ending/peak times are stored as columns in a matrix.
% The ripples are aligned over their peak times (for plotting only).
%
% Created by Jiannis Taxidis, Caltech, USA, February 2013
%
% THIS VERSION MODIFIED FOR SLEEP ANALYSES (10/2017)
%--------------------------------------------------------------------------
%% Select Stage Scored File:
working_dir=pwd;
current_dir='C:\SleepData';
cd(current_dir);
scoredCheck = 0;
while isequal(scoredCheck, 0)
    [scoredFile, scoredPath] = uigetfile({'*.xls','Excel 1997-2003 File (*.xls)'},...
        'Select the Sleep Scored File');
    if isequal(scoredFile,0) || isequal(scoredPath,0)
        uiwait(errordlg('You need to select a file. Please try again',...
            'ERROR','modal'));
    else
        cd(working_dir);
        stageScoredFile= fullfile(scoredPath, scoredFile);
        %Load sleep scored file:
        try
            [numData, stringData] = xlsread(stageScoredFile);
            scoredCheck = 1;
        catch %#ok<*CTCH>
            % If file fails to load, it will notify user and prompt to
            % choose another file.
            uiwait(errordlg('Check if the scored file is saved in Microsoft Excel format.',...
             'ERROR','modal'));
         scoredCheck = 0;
        end

    end
end

%% Detect if states are in number or 2-letter format:
if isequal(size(numData,2),3)
    scoredStates = numData(:,2:3);
    clear numData stringData
else
    scoredStates = numData(:,2);
    clear numData
    stringData = stringData(3:end,3);
    [stateNumber] = stateLetter2NumberConverter(stringData);
    scoredStates = [scoredStates stateNumber];
    clear stateNumber stringData
end
epochInSeconds = scoredStates(2,1) - scoredStates(1,1);
startTime = scoredStates(1,1) * 10^6;
endTime = (scoredStates(end,1) + epochInSeconds) * 10^6;

%% Load Neuralynx continuous file
[CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
        'Pick CSC files.'},'Select Continuously Sampled Channel File');
cscFile = fullfile(CSCFilePath, CSCFilename);
[Timestamps, Samples] = Nlx2MatCSC(cscFile, [1 0 0 0 1], 0, 4, [startTime endTime]);

% Reshape LFP array:
[m1,n1]=size(Samples);  
newM = m1*n1;
nsamp = m1;
LFP = reshape(Samples, newM, 1);
clear Samples

% Fill in the time stamps:
eelen=length(Timestamps);
precise_timestamps = zeros(eelen*nsamp, 1);
idx = 1;
for i = 1:eelen
  if i < eelen
    t1 = Timestamps(i);
    t2 = Timestamps(i+1);
    interval = (t2-t1)/nsamp;
    trange =([t1 : interval : t2-interval]);
    precise_timestamps(idx:idx+nsamp-1,1) = trange;
  else
    t1 = Timestamps(i);
    t2 = t1+interval*nsamp;
    trange =([t1 :interval : t2-interval]);
    precise_timestamps(idx:idx+nsamp-1,1) = trange;
  end
  idx = idx + nsamp;
end
clear Timestamps eelen t1 t2 interval trange
% Convert from usec to seconds:
time = precise_timestamps/1000000;
clear precise_timestamps
Fs = 1/(time(2)-time(1));
% refshank = refs(1);                                                         % The first element of "refs" is the reference shank
% refchannel = refs(2);                                                       % the second one is the reference channel

%% BANDPASS LFP in ripple band
Fripple = [120,200];                                                        % Ripple frequency range to bandpass LFP
LFPbp = bandpass_c(LFP, Fs, Fripple(1), Fripple(2), 'pass');                % Filter the LFP

%% PLOT POWER SPECTRUM
figure;
seg = 10;
[f,~,~] = sp2a2_m1(0,LFP,LFP,Fs,seg,'r2 t2 h');

subplot(4,5,[11,16])
hold on
plot(f(:,1), (f(:,2)+log10(2*pi)),'Linewidth',1.5,'Color','k')
xlabel('Frequency (Hz)','Fontsize',20)
xlim([0 500])
set(gca,'Fontsize',16)

%% CALCULATE ROOT-MEAN-SQUARE OVER 50% OVERLAPPING BINS
bindur = 6;                                                                 % Bin duration for RMS (4.8 msec) --> Now 5.8977 ms
disp(['The bin duration for the rms calculation is: ',num2str(bindur/Fs*1000),' msec.']);
N = length(LFP);
bins = 1 : bindur/2 : (N-bindur);                                           % Array of locations of bin starting times
RMS = zeros(length(bins),1);
for i = 1:length(bins)
    RMS(i) = sum(LFPbp(bins(i):(bins(i)+bindur)).^2);                       % Calculate rms for each bin (IT IS THE SUM OF SQUARES ONLY!!! units = mV^2)
end

%% DETECT ALL SEGMENTS ABOVE THE LFP THRESHOLD AND THEIR LIMITS
crit1 = 2;     
crit2 = 1;
stdr = std(RMS);                                                            % Standard deviation of rms
CRIT1 = crit1*stdr;                                                         % Set the detection criterion
CRIT2 = crit2*stdr;                                                         % Set the limits detection criterion
clear crit1 crit2

RMS(1) = 0;                                                                 % Set the first and last points to zero
RMS(end) = 0;                                                               % to avoid detecting on the edges.

cross = (RMS >= CRIT1);                                                     % Find the array locations over the threshold
crossbd = find((cross - circshift(cross,1)) == 1);                          % Find locations whose preceding location is below the threshold (ripple beginning)
crossbd = [crossbd find(cross - circshift(cross,1) == -1)];                 % Add locations whose following location is below the threshold (ripple end)

time_b = time(bins + bindur/2);
ripples = zeros(length(crossbd(:,1)),3);

for i = 1:length(crossbd(:,1))                                              % For each ripple beginnning
    rtmp = RMS(crossbd(i,1):crossbd(i,2));                                  % Keep the rms
    rippeak = find(rtmp == max(rtmp)) + crossbd(i,1);                       % Find the overall location of the rms peak
    ripples(i,1) = find(RMS(1:rippeak) < CRIT2,1,'last');                   % Keep the last location before ripple beginnning below the limit threshold
    ripples(i,2) = find(RMS(rippeak:end) < CRIT2,1,'first') + rippeak - 1;  % and the first location after ripple end below the limit threshold
    ripples(i,3) = rippeak;                                                 % and the location of the ripple peak
end

[~,m] = unique(ripples(:,1), 'rows');                                       % Keep unique events by checking repeating ripple beginning times
ripples = ripples(m,:);
[~,m] = unique(ripples(:,2), 'rows');                                       % Keep unique events by checking repeating ripple end times
ripples = ripples(m,:);

ripple.startIdx = ripples(:,1);
ripple.stopIdx = ripples(:,2);
ripples = time_b(ripples);
disp('Initially :')
[~,rdur] = ripple_stats(ripples);

%% ASSIGN scored states to each ripple
ripple.states = zeros(rN,1);
for i = 1:size(scoredStates, 1)
    if isequal(i, size(scoredStates, 1))
        subIntervalIndx = find(ripples(:,3) >= scoredStates(i,1) & ripples(:,3) < (scoredStates(i,1) + 10));
    else
        subIntervalIndx = find(ripples(:,3) >= scoredStates(i,1) & ripples(:,3) < scoredStates(i+1,1));
    end
    if ~isempty(subIntervalIndx)
        ripple.states(subIntervalIndx) = scoredStates(i,2);
        clear subIntervalIndx
    end
end

%% KEEP only NREM ripples
targetRipples = ripple.states == 2| ripple.states == 6;                     % NREM and Intermediate Stage (IS)
ripples = ripples(targetRipples, :);
ripple.startIdx = ripple.startIdx(targetRipples);
ripple.stopIdx = ripple.stopIdx(targetRipples);
ripple.states = ripple.states(targetRipples);                               % Preserving distinction between NREM and IS sleep
disp('After removal of waking and REM: ')
[rN,~] = ripple_stats(ripples);

%% REMOVE SHORT AND OVERLAPPING RIPPLES
mindur = 0.02;                                                              % Minimum accepted ripple duration (20 msec)
mindist = 0.05;                                                             % Minimum accepted distance between ripples (50 msec)

rshort = (rdur < mindur);                                                   % If any ripple lasts less than 20 msec discard it
ripples(rshort,:) = [];
ripple.startIdx(rshort) = [];
ripple.stopIdx(rshort) = [];
ripple.states(rshort) = [];

k = ripples(2:end,1) - ripples(1:(end-1),2);                                % Find the distances between all ripples
k = (k < mindist);                                                          % Find when the distance is not at least 50 msec (covers overlaps too)
k1 = [k; 0];
k2 = [0; k];                                                                % Add a zero for the first ripple that was not included
ripples(k1==1,2) = ripples(k2==1,2);                                        % Unite them (use the ending of the second ripple as end of the first one...
ripples(k2==1,:) = [];                                                      % ...and delete the times of the second one)
ripple.startIdx(k2==1) = [];
ripple.stopIdx(k2==1) = [];
ripple.states(k2==1) = [];

disp('After close-ripples removal: ')
[rN,~] = ripple_stats(ripples);

%% REMOVE THETA-TIME RIPPLES --Not needed for sleep data
% load detect_movement data here
% load(['../Sessions/',session,'/Movement_data.mat']);                        % Load the movement analysis of this session
% critV = 10; % Not needed for sleep recording.
% movein = [];
% moveout = [];
% for r = 1:rN
%     k = find(move_in(:,1) < ripples(r,3),1,'last');                         % Find the last incoming motion segment beginning before the ripple
%     if move_in(k,2) > ripples(r,3)                                          % If the segment's end is after the ripple then the ripple is inside the segment
%         movein = [movein, r];                                               % Store it
%     end
%     k = find(move_out(:,1) < ripples(r,3),1,'last');                        % Same for outgoing segments
%     if move_out(k,2) > ripples(r,3)
%         moveout = [moveout, r];
%     end
% end
% move = unique([movein,moveout]);                                            % Unite all motion-ripples
% ripples(move,:) = [];                                                       % and remove them
% 
% immob = (V <= critV);                                                       % Find where the velocity is below the threshold
% move = [];
% for r = 1:size(ripples,1)
%     z = find(timeV <= ripples(r,3),1,'last');                               % Find the last velocity-timepoint before the ripple
%     t1 = immob(z);                                                          % Check if it is immobile
%     z = find(timeV >= ripples(r,3),1,'first');                              % Same for the last timepoint after the ripple
%     t2 = immob(z);
%     if t1 == 0 && t2 == 0                                                   % If both are not immobile timepoints
%         move = [move, r];                                                   % Keep the ripple's number
%     end
% end
% ripples(move,:) = [];                                                       % Remove non-immobile ripples
% 
% disp('After running-ripples removal: ')
% [rN,~] = ripple_stats(ripples);

%% PLOT LFPS AND RIPPLE TIMES
% subplot(4,5,1:10)
% hold on;
% plot(timeV, V/100+abs(min(LFP))+max(LFPbp)+max(LFP),'k')
% line([time(1),time(end)], [critV,critV]/100+abs(min(LFP))+max(LFPbp)+max(LFP),'Color','k')
% plot(time, LFP+abs(min(LFP))+max(LFPbp),'b');
% plot(time, LFPbp,'b');
% plot(time_b, RMS,'k','LineWidth',1.5);
% 
% 
% line([time(1),time(end)], [CRIT1,CRIT1],'Color','k')
% line([time(1),time(end)], [CRIT2,CRIT2],'Color','k')
% z = max(V/100)+abs(min(LFP))+max(LFPbp)+max(LFP);
% for i = 1:rN
%     line([ripples(i,1),ripples(i,1)],[min(LFPbp) z],'Color','b');
%     line([ripples(i,2),ripples(i,2)],[min(LFPbp) z],'Color','r');
% end

%% ALIGN RIPPLES
periripdur = 0.15;                                                          % Duration around ripple (300 msec total window length (like Kamran's));                                                        
peririp = ceil(2*periripdur*Fs);                                    
LFPrip = zeros(rN,peririp+1);
LFPbprip = zeros(rN,peririp+1);

rmax = zeros(1,rN);
rout = [];
for r = 1:rN                                                                % For each event
    lfptmp = LFPbp(ripple.startIdx(r):ripple.stopIdx(r));                   % Keep the corresponding bandpassed signal
    [~,rm] = max(lfptmp);                                                   % Find the position of its maximum value
    rmax(r) = ripple.startIdx(r) - 1 + rm;                                  % Store the overall array position of the peak
    if rmax(r)-peririp/2 <= 0 || rmax(r)+peririp/2 > N                      % If the perievent duration hits the LFP boundaries
        rout = [rout r];                                                    % Note the number of the event
    end
end

ripples(rout,:) = [];
rmax(rout) = [];
ripple.startIdx(rout) = [];
ripple.stopIdx(rout) = [];
ripple.states(rout) = [];
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
%Request user to name output file:
prompt = {'Enter the filename you want to save it as: (just the name)'};
def = {'Rat#_Day'};
dlgTitle = 'Save .MAT file';
lineNo = 1;
answer = inputdlg(prompt,dlgTitle,lineNo,def);
filename = char(answer(1,:));
save(['Ripples', filename, '.mat'], 'ripples','LFPbp','time','time_b',...
    'LFPrip','LFPbprip','rmax','periripdur', 'ripple', 'scoredFile', 'CSCFilename');