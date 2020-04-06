function [rep_cand,r_high,r_low, numCellsInSeq] = candidateSpindleReplays_11102017(spikeData,timeBins,crit)

%--------------------------------------------------------------------------
% Find all spike sequences around each ripple and store the sequence of
% first spikes, mean spike times, and median spikes. Also return the cases
% where the cell spiking is too low (according to "crit") or high enough.
%--------------------------------------------------------------------------

numEvents = size(timeBins,1);                                               % Number of time bins of target events (e.g. spindles, ripples)
cells = size(spikeData,1);                                                  % Number of place cells (incoming or outgoing)
if isequal(cells, 1)                                                        % Make sure cell array has cells in rows.
    spikeData = spikeData';
    cells = size(spikeData,1);
end

maxsp = crit(1); %Orginal =max(crit(1),crit(2)*cells);                                         % Maximum between minimum amount of cells or minimum percentage of cells for accepting a replay
rep_cand = cell(numEvents,1);                                               % Cell array with the candidate replays
r_high = [];                                                                
r_low = [];

numCellsInSeq = zeros(numEvents,1);
last_spike1 = 0;                                                            % Time of the last spike of the previous sequence

for n = 1:numEvents                                                         % For each time bin
    spikes1 = zeros(1,cells);                                               % Matrix with the first-spike sequences
    spikesm = zeros(1,cells);                                               % Matrix with the mean-spike-time sequences
    spikesmd = zeros(1,cells);                                              % Matrix with the median-spike sequences
    
    for c = 1:cells                                                         % For each place cell
        spikes = spikeData{c}/1000000;                                              % Keep all its spikes
        spikes = spikes(spikes >= timeBins(n,1) & spikes <= timeBins(n,2)); % Keep all its spikes that fall within the target interval
        if ~isempty(spikes)                                                 % If there are any spikes
            spikes1(c) = spikes(1);                                         % Keep the first spike
            spikesm(c) = mean(spikes);                                      % Keep the mean spike time
            spikesmd(c) = median(spikes);                                   % Keep the median spike
        end
    end
    
    [spikes1sort, cells1sort] = sort(spikes1);                              % Sort the first spikes in ascending order and get the cell sequence corresponding to this order
    [spikesmsort, cellsmsort] = sort(spikesm);                              % Same for mean spike times
    [spikesmdsort, cellsmdsort] = sort(spikesmd);                           % and median spikes.
    
    cells1sort(spikes1sort == 0) = [];                                      % Remove the cells that did not spike
    spikes1sort(spikes1sort == 0) = [];                                     % Remove their first spikes too
    cellsmsort(spikesmsort == 0) = [];                                      % Same for mean spike time
    spikesmsort(spikesmsort == 0) = [];                                     %
    cellsmdsort(spikesmdsort == 0) = [];                                    % And for median spike
    spikesmdsort(spikesmdsort == 0) = [];                                   %
    
    numCellsInSeq(n) = size(spikes1sort,2);
    
    if length(cells1sort) >= floor(maxsp) & spikes1sort(1) > last_spike1    % If the number of cells that spiked is above threshold and the first spike of the sequence is after the last spike of the previous sequence
        rep_cand{n} = [cells1sort' spikes1sort' cellsmsort'...
            spikesmsort' cellsmdsort' spikesmdsort'];                       % Keep the first-spike sequence and the first spikes etc...
        r_high = [r_high, n];                                               % Keep the number of the ripple as a high-participation ripple (candidate replay)
        last_spike1 = spikes1sort(end);                                     % And update the "last spike of the previous replay"
    else                                                                    % Otherwise 
        rep_cand{n} = [];                                                   % Leave the replay cell empty
        r_low = [r_low, n];                                                 % and keep it as a low-participation ripple
    end
end