function plot_data(session,LFP,event_flag,time_lims,ripple_refs)

%--------------------------------------------------------------------------
% Plots the LFP from shank "refshank" over the time window given as a
% 2-entries vector "time" in sec. Also draws lines over the time limits of
% the user-defined events ("event_flag") and ripples. The spike times of
% all incoming-encoging place cells are shown above the LFP (red) and those
% of the outgoing-encoding cells are shown below (black).
%
%
% Created by Jiannis Taxidis, Caltech, USA, November 2012
%--------------------------------------------------------------------------

%% INITIALIZE
if nargin < 5
    ripple_refs = LFP;
end
refshank = ripple_refs(1);                                                         % The first element of "refs" is the reference shank
refchannel = ripple_refs(2);                   

shank = LFP(1);
channel = LFP(2);

%---------
% Load LFP, Place cell spike data, and replay data here
% LFP = .....;                                % Load the LFP of the reference shank
% load(['PlSpikes.mat']);                    % Load spikes of incoming template
% load(['Replays',num2str(refshank),'_',num2str(refchannel),'.mat']);
%---------

cells = length(PlSpikeData)                                              % Number of incoming-cells


Fs = 1252;
timepoint = round(time_lims*Fs);                                                 % Turn time limits from sec to arraypoints
if timepoint(2) > size(LFP,1)
    error(['Maximum time limit = ',num2str(size(LFP,1)/Fs)])
end
LFP = LFP(timepoint(1):timepoint(2),channel);                               % Get the LFP within the given time points

ripples = Event_times{1};
events = Event_times{event_flag};
events = [events(:,3)-periripdur , events(:,3)+periripdur];                 % Peri-rip times 
k = (events(:,1) > time_lims(1) & events(:,2) < time_lims(2));                        % Find event times within the time window
find(k)
if isempty(k)
    return
end

events = events(k,:);
k = (ripples(:,1) > time_lims(1) & ripples(:,2) < time_lims(2));                        % Find event times within the time window
ripples = ripples(k,:);

%% PLOT
figure;
hold on
plot(time_lims(1):1/Fs:time_lims(2), LFP);                                            % Plot the LFP over time

for c = 1:cells                                                           % For each incoming cell
    k = (PlSpikeData{c} > time_lims(1) & PlSpikeData{c} < time_lims(2));              % Find its spikes within the time window
    spikes{c} = PlSpikeData{c}(k);
end

[minl,maxl] = plot_spikes(spikes,1:cells,LFP);

for ev = 1:size(events,1)                                                   % For each event within the time window
    line([events(ev,1) events(ev,1)],[minl-0.1*cells maxl],'Color','k'); % Draw lines around its borders
    line([events(ev,2) events(ev,2)],[minl-0.1*cells maxl],'Color','k');
end
for ev = 1:size(ripples,1)                                                  % For each ripple within the time window
    line([ripples(ev,1) ripples(ev,2)],[max(LFP) max(LFP)],'Color','k','LineWidth',2);% Draw a line over the LFP along its duration
end



