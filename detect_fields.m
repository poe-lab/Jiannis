function detect_fields(session,inshanks)

%--------------------------------------------------------------------------
% Detects place fields from all cells in "inshanks".
% Finds the outermost x-position limits of the motion and splits each
% directed motion into position bins. Counts the spikes falling in each bin
% and normalized by the time spent in the bin. Place fields are calculated
% based on the smoothed average firing rate of each cell in each direction.
% Only cells with significantly above zero firing peaks are considered (and
% with peaks below "rmax"). Also cells whose rate is not modulated enough
% around the peak are discarded.
%
% Created by Jiannis Taxidis, Caltech, USA, February 2013
%--------------------------------------------------------------------------

%% INITIALIZE

%------------
% load your spike files here, and the detect_movement output file
% SpikeData = .....; 
% Here, SpikeData is a cell array, where each cell entry corresonds to 
% each silone probe and contains another cell array where each entry contains the spikes
% of each unit in that probe.
%
% load([Movement_data.mat']);                        % Load file with the directed motion segments
%-------------

Dt = diff(timeV);
Dt = Dt(1);                                                                 % Velocity time step

bin_x = 5;                                                                  % Spatial bin length (cm) for spike counting
pval = 0.05;

disp(' ')
disp(['The bin length for the position of firing rate peaks is: ',num2str(bin_x),' cm.']);
disp(['The shanks used for place cells search are: ',num2str(inshanks),'.']);
disp(' ')

shanks = length(inshanks);                                                  % Keep the number of spikes and cells (some cells are empty entries in the SpikeData cell array)

Rates = cell(shanks,1);                                                  % Cell arrays with the firing rates during incoming motion

cells = zeros(1,shanks);
for sh = 1:shanks
    cells(sh) = length(SpikeData{inshanks(sh)}) - 2;                        % number of cells in the shank (EXCLUDE MUA AND BAD UNITS! - first two entries)
    Rates{sh} = cell(cells(sh),1);
end
allcells = sum(cells);
disp(['Initial number of cells = ',num2str(allcells)])

mRates = cell(allcells,1);                                               % and corresponding average rates
sigin = zeros(allcells,1);                                                  % Matrix with above-zero significance of incoming rate peak

cellvector = [];
for sh = 1:shanks
    cellvector = [cellvector; [inshanks(sh)*ones(cells(sh),1) , [1:cells(sh)]'] ];   % Create a 2-column cell-location vector with the shank number and cell number of each cell
end                                                                                        % for locating each incoming place cell that will be detected


%% GET OVERALL X-POSITION BINS
evs = size(move_in,1);                                                   % Number of incoming segments 
disp([num2str(evs),' incoming and ',num2str(evs_out),' outgoing sessions'])

k = zeros(evs,2);
for i = 1:evs                                                            % For each incoming segment
    k(i,:) = [pos_in{i}(1) pos_in{i}(end)];                                 % Keep the outermost x-positions (limits)
end
xlims = [min(k(:,1)), max(k(:,2))];                                      % Find the total x-position limits of the incoming motion

bins = xlims(1) : bin_x : xlims(2);                                         % Create bin vector with bin starting times for incming segments
bins = bins(2:end-1);

disp(['Position limits for all sessions are: ',num2str(bins(1)),' to ', num2str(bins(end)),' cm.'])


%% GET FIRING RATES AND RESHAPE IN COLUMNS
count = 1;
for sh = 1:shanks                                                           % For each shank
    Spikes = SpikeData{inshanks(sh)}(3:end);                                % Keep only the spikes from the included shanks and EXCLUDE MUA AND BAD UNITS!
    
    for c = 1:cells(sh)                                                     % For each cell
        Rates{sh}{c} = firing_rates(Spikes{c},pos_in,move_in,bins,Dt);   % Compute its x-positions/times/firing rates for position bins over every incoming segment
        
        mRates{count} = mean(Rates{sh}{c}(:,5,:),3);                  % Get its mean rate over all incoming segments
        
        k = find(mRates{count} == max(mRates{count}));                % Find the location of the average incoming rate peak
        sigin(count) = ttest(squeeze(Rates{sh}{c}(k(1),5,:)),0,pval,'right'); % Check if it is significantly above zero
         
        count = count + 1;
    end
end

%% SMOOTH AVERAGE RATES
for c = 1:allcells                                                             % For each cell
    mRates{c} = smooth(mRates{c},3,'moving');                                    % Smooth its average firing rate over incoming segments
end

%% REMOVE INTERNEURONS
%------
% INs = IN_clusters(session);   % READ THE INTERNEURON INDEXES
% ------

[~,kout] = intersect([cellvector(:,1), cellvector(:,2) + 2],INs,'rows');
mRates(kout) = [];                                                      
cellvector(kout,:) = [];
sigin(kout) = [];


allcells = allcells - length(kout);

disp('After interneurons removal:')
disp(['Number of incoming place cells = ',num2str(length(mRates))])

%% REMOVE SPARSE SPIKING
kin = [];
for c = 1:allcells                                                             % For each cell
    if  sigin(c) ~= 1 %|| mean(mRates_in{c}) > rmax %any(mRates_in{c} > rmax)                           % If its incoming firing rate peak isnt significantly nonzero, or if it exceed the maximum rate
        kin = [kin, c];                                                     % Keep the cell's number
    end
end
mRates(kin) = [];                                                        % Delete these cells from the incoming rates array
cellvector(kin,:) = [];                                                  % And from the cell-location vector

disp('After sparse firing removal:')
disp(['Number of incoming place cells = ',num2str(length(mRates))])

f1 = figure;
ribbon_plot(mRates,bins,subplot(2,3,1))                                  % Plot all average firing rates

%% REMOVE CONTINUOUS SPIKING (NON-FIELD-MODULATED)
k = [];
for c = 1:length(mRates)                                                 % For each (potential) incoming place cell
    m = mean(mRates{c});                                                 % Find the mean value of its mean firing rate
    if sum(mRates{c} > m) >= sum(mRates{c} <= m) / 2                  % If the values above mean are more than half those below mean (not modulated)
        k = [k, c];                                                         % Keep its number (SAME AS SAYING THAT POINTS ABOVE MEAN MUST BE LESS THAN 1/3 OF TOTAL POINTS!)
    end
end
mRates(k) = [];                                                          % Delete these cells
cellvector(k,:) = [];


disp('After continuous firing removal:')
disp(['Number of incoming place cells = ',num2str(length(mRates))])

ribbon_plot(mRates,bins,subplot(2,3,2))                                  % Plot all average firing rates


%% FIND PLACE FIELDS OF EACH CELL AND REMOVE EDGE-PLACE CELLS
[mPlace_rates, plcells, fields] = firing_peaks(mRates, cellvector, bins);                   % Find the peak location of the average firing rate of each incoming cell

lc = size(plcells,1);                                                % Final number of incoming place cells

disp('After edge fields removal:')
disp(['Number of incoming place cells = ',num2str(lc)])

ribbon_plot(mPlace_rates , bins , subplot(2,3,3));                        % Plot all average firing rates


%% SET UP FINAL VECTORS
plcells = [plcells, fields];                                           % Combine the cell locations with the field locations


PlSpikeData = cell(lc,1);                                             % Allocate memory for the spikes of incoming place cells
for c = 1:lc                                                          % For each incoming cell
    PlSpikeData{c} = SpikeData{plcells(c,1)}{plcells(c,2) + 2};           % Keep its spikes
end


Place_rates = cell(lc,1);
for c = 1:lc
    Place_rates{c} = Rates{find(inshanks == plcells(c,1))}{plcells(c,2)};
end

plcells(:,2) = plcells(:,2) + 2;                                          % Add +2 to each cell number to account for the ignored clusters 0 and 1

plcells


%% PLOT FINAL FIELDS
f2 = figure;
trajectories_plot(pos_in , move_in , Dt , bins , subplot(3,2,1));                   % Plot trajectories and bins
ribbon_plot(mPlace_rates , bins , subplot(3,2,[3,5]));                        % Plot all average firing rates

%% SAVE DATA AND FIGURE
save(['PlaceFields_data.mat'],'Place_rates','mPlace_rates','Rates','plcells','bins');
save(['PlSpikes.mat'],'PlSpikeData','pcells');

%--------------------------------------------------------------------------
%% ------------------------------------------------------------------------








function binned_rates = firing_rates(spikes,pos,times,bins,dt)

%--------------------------------------------------------------------------
% Splits each motion segment into x-position bins and calculates the firing
% rate in each bin (spike number normalized by the time spent on the bin).
% Returns a matrix with the x-limits, the time-limits and the firing rate
% at each x-position bin, stacked in 5 columns. Each event is in the 3rd
% dimension.
%
% Output: binned_rates = (bins x 5 x events)
%--------------------------------------------------------------------------

evs = size(times,1);                                                        % Number of events
binsl = length(bins);                                                       % Number of bins
dx = diff(bins);
dx = dx(1);                                                                 % Bin length (cm)
binned_rates = zeros(binsl,5,evs);

for r = 1:evs                                                               % For each event
    t = times(r,1) : dt : times(r,2);                                       % Keep the time array of the event
    for b = 1:binsl                                                         % For each bin
        lims = find(pos{r} >= bins(b) & pos{r} <= bins(b)+dx);              % Find the entries where the position is inside the bin
        if length(lims) > 1                                                 % If more than one entry (to avoid having one limit and thus zero time-duration (NaN rates))
            xlims = [pos{r}(lims(1)), pos{r}(lims(end))];                   % Keep x-positions at the edges of the bin
            tlims = [t(lims(1)), t(lims(end))];                             % Keep the corresponding time points
            rate = (spikes >= tlims(1) & spikes <= tlims(2));               % Find number of spikes within the bin
            rate = sum(rate) / (tlims(2) - tlims(1));                       % Normalize by time spent in bin
            
            binned_rates(b,:,r) = [xlims, tlims, rate];                     % Store
        end
    end
end
%--------------------------------------------------------------------------


function [rates,vector,peaks] = firing_peaks(rates,vector,bins)

%--------------------------------------------------------------------------
% Locates the position where each cell has its firing rate peak. Rearrages
% cells according to the peak locations.
%--------------------------------------------------------------------------

cells = length(rates);                                                      % Numbe of cells
peaks = zeros(cells,1);                                                     % Vector with the peak positions

dx = diff(bins);
dx = dx(1);                                                                 % Bin length (cm)

for c = 1:cells
    peak = find(rates{c} == max(rates{c}),1,'first');                       % Find the first entry of the maximum rate
    peaks(c) = bins(peak) + dx/2;                                           % Store the corresponding bin middle point
%     centroid = (bins+dx/2)*rates{c} / sum(rates{c});
%     dc = abs(rates{c} - centroid);
%     dc = find(dc == min(dc));
%     dc = dc(1);
%     peaks(c) = bins(dc) + dx/2;
end
[peaks,order] = sort(peaks);                                                % Sort the peak locations
vector = vector(order,:);                                                   % Rearrange cell positions accordingly
rates = rates(order);                                                       % And their rates
%--------------------------------------------------------------------------


function ribbon_plot(data,bins,sub)

Y = [];
for i = 1:numel(data)
    if ~isempty(data{i})
        Y = [Y data{i}];
    else
        Y = [Y zeros(size(Y,1),1)];
    end
end
y = repmat(bins',1,size(Y,2));
subplot(sub)
ribbon(y,Y)
view(-62,30)
colormap hsv
ylabel('Position on track')
xlabel('Cell number')
zlabel('Firing rates (Hz)')
%--------------------------------------------------------------------------


function trajectories_plot(pos,time,Dt,bins,sub)
subplot(sub)
evs = length(pos);
c = zeros(1,evs);
hold on
for r = 1:evs
    t = time(r,1):Dt:time(r,2);                                             % Time array of the segment
    dm = abs(pos{r} - mean(pos{r}));                                        % Get the distance of each x-position from the mean value of the x-position
    tm = t(dm == min(dm));                                                  % Find the timepoint where the x-position is closest to its mean
    plot(t-tm(1) , pos{r},'k');                                             % Plot time vs x-position with the ~mean(x-position) at time-0
    c(r) = abs(t(1) - tm(1));                                               % Keep the time-shift
end
c = 2*mean(c);                                                      
line(repmat([-c;c],1,length(bins)) , [bins; bins], 'Color',[0.8,0.8,0.8])   % Plot all bins
%--------------------------------------------------------------------------


