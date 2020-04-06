function detect_replays(session,refs)

%--------------------------------------------------------------------------
% Detects place cell replays during detected SWRs from a
% reference shank/channel.
% 
%
% Created by Jiannis Taxidis, Caltech, USA, February 2013


% condition for ripple to be considered for replay:(more than 5 cells or 30% of cells spiked, 
% whichever is larger)

% Created by Jiannis Taxidis, Caltech, USA, February 2013
%--------------------------------------------------------------------------

refshank = refs(1);                                                  % The first element of "refs" is the reference shank
refchannel = refs(2);                                                % the second one is the reference channel

%-----------
% load ripple and place cell spike data here
% session = strcat([num2str(session(1)),' ',num2str(session(2)),' ',num2str(session(3))]);
% load(['Ripple_data',num2str(refshank),'_',num2str(refchannel),'.mat']); % Load ripple times
% load(['PlSpikes.mat']);                                                 % Load incomind place cell spikes
%---------

mincells = [5, 0.3];                                                        % Criterion for detecting too little spiking in a replay (at least 5 cells or 30% of cells - whichever is larger- have to spike)
seq_crit = 1;                                                               % Criterion for getting the replay sequence (1 = first spike, 3 = mean spike time, 5 = median spike)
repet = 500;                                                                % Number of random permutations for each replay

hbins = -1 : 0.05 : 1;                                                      % Bins for correlation histograms
rN = size(ripples,1);
disp(' ')
disp([num2str(rN),' ripples in total'])
disp(' ')

%% DETECT SPIKE SEQUENCES
disp('Detecting spike sequences from all ripples')
[Sequences,r_high,r_low] = cand_replays(SpikeData,ripples,periripdur,mincells);     % Get all candidate incoming replays from ripples
disp(' ')

Seq_cand = Sequences;
Seq_cand(r_low) = [];                                                 % Remove replays with low participation

disp([num2str(length(r_high)),' ripples contain candidate replays'])
disp([num2str(length(r_low)/rN*100),'% of ripples contain too low spiking and were removed'])
disp(' ')

%% CREATE RANDOM PERMUTATIONS
disp('Creating random permutations of spike sequences')
Seq_perm = rand_perm(Seq_cand,seq_crit,repet);                        % Create random permutations of each incoming replay
disp(' ')

%% COMPUTE CORRELATIONS
disp('Computing correlations of sequences with place cell sequence')
[Rs,~] = rep_corr(Seq_cand);                                          % Compute correlation of each incoming replay with the place cell sequence
disp(' ')

figure
plot_hist(Rs(:,1),hbins,[0.5,0.5,0.5],'w',[4,4,1])                       % Histogram plots
plot_hist(Rs(:,3),hbins,[0.5,0.5,0.5],'w',[4,4,2])
plot_hist(Rs(:,5),hbins,[0.5,0.5,0.5],'w',[4,4,3])

%% COMPUTE CORRELATIONS OF PERMUTATIONS
disp('Computing correlations of permuted sequences with place cell sequence')
Rs_perm = perm_corr(Seq_perm);                                    % Compute correlation of each permutation with place cell sequence
disp(' ')  

plot_hist(Rs_perm(:),hbins,'k','k',[4,4,4])                              % Histogram plots

%% COMPUTE SIGNIFICANCE AND GET DIRECTION OF REPLAY
disp('Computing significant correlations')
[forw, rev] = compare_corrs(Rs,Rs_perm,seq_crit);               % Find the incoming replays that are significantly correlated with the place cell sequence (get direction)
disp(' ')

Forw = r_high(forw);                                               % Keep the incoming replays that are forward
Rev = r_high(rev);                                                 % and those that are reverse

if ~isempty(intersect(Forw, Rev)) 
    error('Some replays were identified as both forward AND reverse!')
end

disp([num2str(length(unique([Forw,Rev]))),' ripples contain significant replays'])
disp(' ')

Rs_f = Rs(forw,seq_crit);                                          % Keep the correlations of the incoming forward replays
Rs_r = Rs(rev,seq_crit);                                           % and the reverse ones

plot_hist(Rs_f,hbins,'k','k',[4,4,9])                                    % Histogram plots
plot_hist(Rs_r,hbins,'k','k',[4,4,10])

%% FIND "DOUBLE-IDENTITY" REPLAYS
% You dont need that since you only have unidirectional motion

%% REMOVE "DOUBLE IDENTITY" REPLAYS 
% You dont need that since you only have unidirectional motion


%% SET UP FINAL MATRICES
disp([num2str(length(unique([Forw,Rev]))),' ripples contain significant replays'])
disp([num2str(length(Forw)),' are forward replays of incoming cells'])
disp([num2str(length(Rev)),' are reverse replays of incoming cells'])
disp(' ')  

plot_hist(Rs_f,hbins,'w','k',[4,4,11])                                   % Histogram plots
plot_hist(Rs_r,hbins,'w','k',[4,4,12])

Replays_f = ripples(Forw,:);                                          % Keep the final ripple times that correspond to forward incoming replays
Replays_r = ripples(Rev,:);                                           %

Reps = sort([Forw,Rev]);
Replays_all = ripples(Reps,:);           % And for all replays together

Noreps = 1:rN;
Noreps(Reps) = [];
No_replays = ripples(Noreps,:);

Sequences_f = Sequences(Forw);                                     % Keep the final sequences that correspond to forward incoming replays
Sequences_r = Sequences(Rev);                                      %
 
Sequences_noreps = Sequences(Noreps); 

Event_indexes = {1:rN; Noreps ; Reps; Forw; Rev};
Event_times = {ripples; No_replays; Replays_all; Replays_f; Replays_r};

%% SAVE
save(['Replays',num2str(refshank),'_',num2str(refchannel),'.mat'],'Event_times','Event_indexes',...
    'Sequences_f','Sequences_r','Sequences_noreps','periripdur');
save(['Replays_data',num2str(refshank),'_',num2str(refchannel),'.mat'],'r_high',...
    'Seq_perm','Rs','Rs_perm','Rs_f','Rs_r');

% -------------------------------------------------------------------------
%% ------------------------------------------------------------------------
        






function [rep_cand,r_high,r_low] = cand_replays(SpikeData,ripples,periripdur,crit)

%--------------------------------------------------------------------------
% Find all spike sequences around each ripple and store the sequence of
% first spikes, mean spike times, and median spikes. Also return the cases
% where the cell spiking is too low (according to "crit") or high enough.
%--------------------------------------------------------------------------

rN = size(ripples,1);                                                       % Number of ripples
cells = size(SpikeData,1);                                                  % Number of place cells (incoming or outgoing)
maxsp = max(crit(1),crit(2)*cells);                                         % Maximum between minimum amount of cells or minimum percentage of cells for accepting a replay
rep_cand = cell(rN,1);                                                      % Cell array with the candidate replays
r_high = [];                                                                
r_low = [];

last_spike1 = 0;                                                            % Time of the last spike of the previous sequence

for r = 1:rN                                                                % For each ripple
    spikes1 = zeros(1,cells);                                               % Matrix with the first-spike sequences
    spikesm = zeros(1,cells);                                               % Matrix with the mean-spike-time sequences
    spikesmd = zeros(1,cells);                                              % Matrix with the median-spike sequences
    
    for c = 1:cells                                                         % For each place cell
        spikes = SpikeData{c};                                              % Keep all its spikes
        spikes = spikes(spikes >= ripples(r,3)-periripdur & spikes <= ripples(r,3)+periripdur); % Keep all its spikes that fall within the peri-ripple interval
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
    
    if length(cells1sort) >= floor(maxsp) & spikes1sort(1) > last_spike1           % If the number of cells that spiked is above threshold and the first spike of the sequence is after the last spike of the previous sequence
        rep_cand{r} = [cells1sort' spikes1sort' cellsmsort' spikesmsort' cellsmdsort' spikesmdsort']; % Keep the first-spike sequence and the first spikes etc...
        r_high = [r_high, r];                                               % Keep the number of the ripple as a high-participation ripple (candidate replay)
        last_spike1 = spikes1sort(end);                                     % And update the "last spike of the previous replay"
    else                                                                    % Otherwise 
        rep_cand{r} = [];                                                   % Leave the replay cell empty
        r_low = [r_low, r];                                                 % and keep it as a low-participation ripple
    end
end
% -------------------------------------------------------------------------


function rep_perm = rand_perm(rep_cand,crit,repetitions)

%--------------------------------------------------------------------------
% Creates random permutations of the candidate replays based on the "crit"
% criterion (1 = first spike, 3 = mean spike time, 5 = median spike).
%--------------------------------------------------------------------------

if ~any(crit == [1,3,5])
    error('Choose one of the three criteria for the permutation!')
end

rh = length(rep_cand);                                                      % Number of candidate replays
rep_perm = cell(rh,1);                                                      % Cell array with the permutations of each replay (in a cell)

for r = 1:rh                                                                % For each replay
    cells1sort = rep_cand{r}(:,crit);                                       % Get the sequence of cells (based on the selected criterion)
    rep_perm{r} = zeros(length(cells1sort) , repetitions);                  % Allocate memory for the permutations (stacked in columns)
    for i = 1:repetitions                                                   % For each repetition
        z = randperm(length(cells1sort));                                   % Randomly permute the ENTRIES of the cell sequence
        rep_perm{r}(:,i) = cells1sort(z);                                   % Store the shuffled cells
    end
end
%% ------------------------------------------------------------------------


function [Rs,Ps] = rep_corr(rep_cand)

%--------------------------------------------------------------------------
% Computes Spearman rank-order correlations between the candidate replays
% and the place cell sequence based on first spike, mean spike time and 
% median spike.
%--------------------------------------------------------------------------

rh = length(rep_cand);                                                      % Number of candidate replays
Rs = zeros(rh,5);                                                           % Matrix with the correlations (for first spike, mean spike time, median spike)
Ps = zeros(rh,5);                                                           % Matrix with the p-values

for r = 1:rh                                                                % For each replay
    cells1sort = rep_cand{r}(:,1);                                          % Get the cells that participated, sorted in ascending order of first spike
    cellsmsort = rep_cand{r}(:,3);                                          % Get the cells that participated, sorted in ascending order of mean spike time
    cellsmdsort = rep_cand{r}(:,5);                                         % Get the cells that participated, sorted in ascending order of median spike
    
    [rs,p] = corr([cells1sort , cellsmsort , cellsmdsort], [sort(cells1sort), sort(cellsmsort), sort(cellsmdsort)], 'type','Spearman'); % Get the correlation of the sequence with the sequence of sorted cell numbers
    
    Rs(r,:) = [rs(1,1), 0, rs(2,2), 0, rs(3,3)];                                  % Store the autocorrelations (diagonal
    Ps(r,:) = [p(1,1), 0, p(2,2), 0, p(3,3)];                                     % And the p-values.
end
%% ------------------------------------------------------------------------


function Rperm = perm_corr(rep_perm)

%--------------------------------------------------------------------------
% Computes Spearman rank-order correlations between the shuffled replays
% and the place cell sequence.
%--------------------------------------------------------------------------

rh = length(rep_perm);                                                      % Number of candidate replays
repetitions = size(rep_perm{1},2);                                          % Number of shuffle repetitions

Rperm = zeros(rh,repetitions);                                              % Matrix with the correlations
% Pperm = zeros(rh,repetitions);                                              % Matrix with the p-values

for r = 1:rh                                                                % For each replay
    cellssorted = sort(rep_perm{r}(:,1));                                   % Get the cells that participated, sorted in ascending order
    R = repmat(cellssorted,1,repetitions);
    R = corr(rep_perm{r},R, 'type','Spearman'); % Get the correlation of the shuffled and the prototype sequence
    Rperm(r,:) = diag(R)';
end
%% ------------------------------------------------------------------------


function [forw, rev] = compare_corrs(Rs,Rperm,crit)

%--------------------------------------------------------------------------
% Checks if the correlation of each sequence (based on "crit") is
% significant (above or below 95% of the random permutations correlations)
%--------------------------------------------------------------------------

if ~any(crit == [1,3,5])
    error('Choose one of the three criteria for the permutation!')
end

[rh,repet] = size(Rperm);                                                   % Number of candidate replays and repetitions
forw = [];
rev = [];

for r = 1:rh                                                                % For each replay
    Rless = sum(Rperm(r,:) <= Rs(r,crit));                                  % Find how many permutations of the replay have less correlation with the prototype sequence than the original sequence
    Rless = Rless/repet*100;                                                % Turn it to percentage
    Rmore = sum(Rperm(r,:) >= Rs(r,crit));                                  % Same for permutations with higher correlation than the original
    Rmore = Rmore/repet*100;                                                % 
    
    if Rless >= 95                                                          % If there are more than 95% permuations with correlation less than the original sequence
        forw = [forw; r];                                                   % Mark this replay as a forward replay
    elseif Rmore >= 95                                                      % If there are more than 95% permutations with correlation higher than the original sequence
        rev = [rev; r];                                                     % Mark this replay as a reverse replay
    end
end
%% ------------------------------------------------------------------------
