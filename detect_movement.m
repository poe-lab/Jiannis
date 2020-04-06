function detect_movement(session,refs)

%--------------------------------------------------------------------------
% Detects the time segments when the rat is moving along the linear track,
% based on its x-position. Also calculates its velocity and the direction
% of movement. Plots the x-position, velocity, LFP and time
% segments corresponding to each transversing the linear track. Also plots
% power spectrum of total LFP and of movement-sections.
%
% Created by Jiannis Taxidis, Caltech, USA, February 2013
%--------------------------------------------------------------------------

refshank = refs(1);                                                         % The first element of "refs" is the reference shank
refchannel = refs(2);                                                       % the second one is the reference channel (for LFP plotting)
Fs = 1252;

critV1 = 3;                                                                 % Criterion for movement detection (ratio of velocity mean)
critV2 = 10;                                                                % Criterion of immobility
mindur = 1;                                                                 % Minimum accepted duration of a transversal
mindist = 1;                                                                % Minimum accepted time distance between two transversals


figure;

%% LOAD LFPS (ONLY FOR PLOTTING)
%-------
% load your LFP here
% LFP = ....;                                                    % Keep the LFP from the reference channel
%------

N = length(LFP);
time = (0:N-1)/Fs;                                                          % Timescale for LFP plotting

%% CALCULATE RAT'S VELOCITY AND POSITION

%-----
% Load your position file here
% x = ....;
% y = ....;
%-----

F = 60;                                                                     % Recoringd frequency (Hz)
Dt = 1/F;                                                                
V = sqrt(diff(x).^2+diff(y).^2);                                            % Calculate velocity;
V = [0; V];                                                                 % Add an initial zero since the first point was lost by diff(x)
V(end) = 0;                                                                 % Set last velocity poitn to zero to avoid detecting on the edges.
V = V/Dt;                                                                   % Velocity (cm/sec)

timeV = 0 : Dt : (length(V)-1)*Dt;                                              % Create a time array

V = bandpass_c(V,60,0,1,'pass');

critV1 = critV1*mean(V);                                                    % Movement detection criterion

%% DETECT ALL SEGMENTS ABOVE THE THRESHOLD AND THEIR LIMITS
cross = (V >= critV1);                                                      % Find the array locations over the threshold
crossbd = find((cross - circshift(cross,1)) == 1);                          % Find locations whose preceding location is below the threshold (movement beginning)
crossbd = [crossbd find(cross - circshift(cross,1) == -1)];                 % Add locations whose following location is below the threshold (movement end)
rN = length(crossbd(:,1));                                                  % Number of detected segements

move = zeros(rN,3);                                                         % Matrix with the beginning/end/velocity peak array-locations of movement
for i = 1:rN                                                                % For each segment
    Vtmp = V(crossbd(i,1):crossbd(i,2));                                    % Keep the velocity
    Vpeak = find(Vtmp == max(Vtmp),1,'first') + crossbd(i,1)-1;             % Find the overall location of the velocity peak
    move(i,1) = find(V(1:Vpeak) < critV2,1,'last');                         % Keep the last location before segment-beginnning below the limit threshold
    move(i,2) = find(V(Vpeak:end) < critV2,1,'first') + Vpeak - 1;          % and the first location after segment-end below the limit threshold
    move(i,3) = Vpeak;                                                      % and the location of the velocity peak
end

[~,m] = unique(move(:,1), 'rows');                                          % Keep unique events by checking repeating segment beginning entries
move = move(m,:);
[~,m] = unique(move(:,2), 'rows');                                          % Keep unique events by checking repeating segment end entries
move = move(m,:);

move_t = timeV(move);                                                       % Get time points from array locations

disp('Initially :')
[~,rdur] = theta_stats(move_t);                                             % Report

%% REMOVE SHORT AND OVERLAPPING SEGMENTS
rshort = (rdur < mindur);                                                   % If any segment lasts less than 1 sec discard it
move_t(rshort,:) = [];
move(rshort,:) = [];                                                        % Remove corresponding entries as well

k = move_t(2:end,1) - move_t(1:(end-1),2);                                  % Find the time distances between all segments
k = (k < mindist);                                                          % Find when the distance is not at least 1 sec (covers overlaps too)
k1 = [k; 0];                                                                % Keep locations of the preceding segment
k2 = [0; k];                                                                % and of the following segment
move_t(k1==1,2) = move_t(k2==1,2);                                          % Unite them (use the ending of the following segment as end of the preceding one...
move_t(k2==1,:) = [];                                                       % and delete the times of the following one)
move(k1==1,2) = move(k2==1,2);                                              % Same for entry-list
move(k2==1,:) = [];                                                         %

disp('After close-ripples removal: ')
rN = theta_stats(move_t);                                                   % Report

%% GET DIRECTION OF MOTION AND REMOVE AMBIGUOUS SEGMENTS
dir = zeros(rN,2);                                                          % Logical matrix with direction of motion
for r = 1:rN                                                                % For each segment
    t1 = move(r,1);                                                         % Keep segment beginning entry
    t2 = move(r,2);                                                         % Ending entry
    t3 = move(r,3);                                                         % Velocity peak entry
    if x(t1) < x(t3) && x(t3) < x(t2)                                       % If the x-position at the peak is larger than the beginning and smaller than the ending one
        dir(r,:) = [1 0];                                                   % Mark the segment as incoming
    elseif x(t1) > x(t3) && x(t3) > x(t2)                                   % If it is smaller than the beginning and larger than the ending one
        dir(r,:) = [0 1];                                                   % Mark the segment as outgoing
    end
end
amb = (sum(dir,2) == 0);                                                    % Find the ambiguous segments that were not classified
move_t(amb, :) = [];                                                        % Remove them from time-list
move(amb, :) = [];                                                          % and from entry list
dir(amb, :) = [];                                                           % and from direction list


%% GET FINAL TIMES AND POSITIONS DURING EVENTS
ins = find(dir(:,1));                                                       % Get segment numbers of incoming motion
outs = find(dir(:,2));                                                      % Get segment numbers of outgoing motion
move_in = move(ins,:);                                                      % Keep corresponding entries from the entry list                       
move_out = move(outs,:);

pos_in = cell(length(ins),1);                                               % Cell with the x-positions of each incoming segment
for r = 1:length(ins)                                                       % For each incoming segment
    t1 = move_in(r,1);                                                      % Keep the entry of the beginning
    t2 = move_in(r,2);                                                      % and of the end
    pos_in{r} = x(t1:t2);                                                   % Store the corresponding x-postions
end
pos_out = cell(length(outs),1);                                             % Same for outgoing segments
for r = 1:length(outs)
    t1 = move_out(r,1);
    t2 = move_out(r,2);
    pos_out{r} = x(t1:t2);
end

move_in = timeV(move_in);                                                   % Finally get the time points of the incoming
move_out = timeV(move_out);                                                 % and outgoing segments

disp('After direction detecting: ')
rN = theta_stats(move_t);                                                   % Report
disp([num2str(length(ins)),' sessions are incoming and ', num2str(length(outs)),' sessions are outgoing'])

%% PLOT LFPS AND MOTION TIMES
subplot(1,5,2:5)
hold on;
plot(timeV, x/400 + max(LFP) + max(V/50) + min(x/400), 'r')
plot(timeV, V/50 + max(LFP),'k')
line([time(1),time(end)], [critV1,critV1]/50 + max(LFP),'Color','k')
line([time(1),time(end)], [critV2,critV2]/50 + max(LFP),'Color','k')
plot(time, LFP,'b');
z = max(LFP) + max(V/50) + max(x/400);
for i = 1:rN
    line([move_t(i,1),move_t(i,1)],[min(LFP) z+0.2],'Color','b');
    line([move_t(i,2),move_t(i,2)],[min(LFP) z+0.2],'Color','r');
end
line([move_in(:,1), move_in(:,2)]', repmat([z, z]+0.2,length(ins),1)','Color','r','LineWidth',1.5)
line([move_out(:,1), move_out(:,2)]', repmat([z, z]+0.2,length(outs),1)','Color','k','LineWidth',1.5)
xlim([time(1),time(end)])

%% PLOT LFP POWER SPECTRUM
th = round(move_t*Fs);                                                      % Get the corresponding time points of the LFP
if th(1) == 0
    th(1) = 1;
end
if th(end) > N
    th(end) = N;
end

LFPth = [];
for r = 1:rN                                                                % For each event (irrelevantly of motion direction)
    LFPth = [LFPth; LFP(th(r,1):th(r,2))];                                  % Keep the corresponding LFP segment
end

seg = 11;
[f,~,~] = sp2a2_m1(0,LFP,LFP,Fs,seg);                             % Total LFP spectrum
subplot(1,5,1)
hold on
plot(f(:,1), (f(:,2)+log10(2*pi)),'Linewidth',1.5,'Color','k')

[f,~,~] = sp2a2_m1(0,LFPth,LFPth,Fs,seg);                         % Motion-LFP spectrum
plot(f(:,1), (f(:,2)+log10(2*pi)),'Linewidth',1.5,'Color','b')
xlabel('Frequency (Hz)','Fontsize',20)
xlim([0 500])
set(gca,'Fontsize',16)

%% SAVE
save(['Movement_data.mat'],'move_in','move_out','pos_in','pos_out','LFPth','time','timeV','V','x');
%--------------------------------------------------------------------------
%% ------------------------------------------------------------------------


function [rN,rdur] = theta_stats(ripples)

rN = size(ripples,1);                                                       % Number of ripples in this session
rdur = ripples(:,2)-ripples(:,1);                                           % and ripple durations
rmax = max(rdur);                                                           % Maximum duration
disp([num2str(rN),' running sessions of average duration ', num2str(mean(rdur)*1000),' msec and maximum duration ',num2str(rmax*1000),' msec',' (SD = ',num2str(std(rdur)*1000),')'])

