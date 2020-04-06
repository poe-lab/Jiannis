function detect_replays_10302017BG%(session,refs)

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

% refshank = refs(1);                                                  % The first element of "refs" is the reference shank
% refchannel = refs(2);                                                % the second one is the reference channel

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
rN = size(EschenkoSpindle.timestamp,1);
disp(' ')
disp([num2str(rN),' spindles in total'])
disp(' ')

%% DETECT SPIKE SEQUENCES IN SPINDLES
[Sequences,r_high,r_low] = candidateSpindleReplays_10302017(spikeData,spindleHalf2,mincells);     % Get all candidate replays
Seq_cand = Sequences;
Seq_cand(r_low) = [];                                                       % Remove replays with low participation

%% CREATE FORWARD AND REVERSE ORDERED SEQUENCES FOR CIRCULAR TRACKS
numSeq = size(Seq_cand,1);
forwardSeq = cell(numSeq,1);
reverseSeq = cell(numSeq,1);

for i = 1:numSeq
    A = Seq_cand{i}(:,seq_crit);
    sortA = sort(A);   
    forwardSeq{i} = [sortA(sortA>=A(1)); sortA(sortA<A(1))];
    sortA = sort(A, 'descend');
    reverseSeq{i} = [sortA(sortA<=A(1)); sortA(sortA>A(1))];
    clear A sortA
end

%% COMPUTE CORRELATIONS IN FORWARD AND REVERSE ORDERED SEQUENCES
%--------------------------------------------------------------------------
% Computes Spearman rank-order correlations between the candidate replays
% and the place cell sequence based on seq_crit (1 = first spike, 3 = mean
% spike time, 5 = median spike)
%--------------------------------------------------------------------------
disp('Computing correlations of sequences with place cell sequence')
Rs = zeros(numSeq,2);                                                       % Matrix with the correlations (for first spike, mean spike time, median spike)
Ps = zeros(numSeq,2);                                                       % Matrix with the p-values

for r = 1:numSeq                                                            % For each replay, get the correlation of the sequence with...
    [rsFwd,pFwd] = corr(Seq_cand{r}(:,seq_crit), forwardSeq{r}, 'type','Spearman','tail', 'right'); % ...forward sorted cell numbers
    [rsRev,pRev] = corr(Seq_cand{r}(:,seq_crit), reverseSeq{r}, 'type','Spearman','tail', 'right'); % ...reverse sorted cell numbers
    
    Rs(r,:) = [rsFwd rsRev];                                                % Store the autocorrelations (diagonal
    Ps(r,:) = [pFwd pRev];                                                  % And the p-values.
end

disp(' ')

figure
% histogram(Rs(:,1),hbins,[0.5,0.5,0.5],'w',[4,4,1])                          % Histogram plot for foward replay correlations
histogram(Rs(:,1),hbins,'FaceColor',[0.5,0.5,0.5])
% histogram(Rs(:,2),hbins,[0.5,0.5,0.5],'w',[4,4,2])                          % Histogram plot for reverse replay correlations
histogram(Rs(:,2),hbins,'FaceColor',[0.5,0.5,0.5])
%% CREATE RANDOM PERMUTATIONS
%--------------------------------------------------------------------------
% Creates random permutations of the candidate replays based on the "crit"
% criterion (1 = first spike, 3 = mean spike time, 5 = median spike).
%--------------------------------------------------------------------------
disp('Creating random permutations of spike sequences')
Seq_perm = cell(numSeq,1);                                                  % Cell array with the permutations of each replay (in a cell)

for r = 1:numSeq                                                            % For each replay
    cells1sort = Seq_cand{r}(:, seq_crit);                                  % Get the sequence of cells (based on the selected criterion)
    Seq_perm{r} = zeros(length(cells1sort), repet);                         % Allocate memory for the permutations (stacked in columns)
    for i = 1:repet                                                         % For each repetition
        z = randperm(length(cells1sort));                                   % Randomly permute the ENTRIES of the cell sequence
        Seq_perm{r}(:,i) = cells1sort(z);                                   % Store the shuffled cells
    end
end
disp(' ')

%% COMPUTE CORRELATIONS OF PERMUTATIONS
%--------------------------------------------------------------------------
% Computes Spearman rank-order correlations between the shuffled replays
% and the place cell sequence.
%--------------------------------------------------------------------------
disp('Computing correlations of permuted sequences with place cell sequence')
RpermFwd = zeros(numSeq,repet);                                             % Matrix with the forward correlations
% Pperm = zeros(rh,repetitions);                                              % Matrix with the p-values
RpermRev = zeros(numSeq,repet);                                             % Matrix with the reverse correlations
for r = 1:numSeq                                                            % For each replay
    for i = 1:repet                                                         % For each permutation
        RpermFwd(r,i) = corr(Seq_perm{r}(:,i),forwardSeq{r}, 'type','Spearman','tail', 'right'); % Get the correlation of shuffled and the forward sequence
        RpermRev(r,i) = corr(Seq_perm{r}(:,i),reverseSeq{r}, 'type','Spearman','tail', 'right'); % Get the correlation of shuffled and the reverse sequence
    end
end
disp(' ')  
% plot_hist(Rs(:,5),hbins,[0.5,0.5,0.5],'w',[4,4,3])
histogram(RpermFwd,hbins,'FaceColor',[0.5,0.5,0.5])
% plot_hist(Rs_perm(:),hbins,'k','k',[4,4,4])                              % Histogram plots
histogram(RpermRev,hbins,'FaceColor',[0.5,0.5,0.5])

%% COMPUTE SIGNIFICANCE AND GET DIRECTION OF REPLAY
%--------------------------------------------------------------------------
% Checks if the correlation of each sequence (based on "crit") is
% significant (above or below 95% of the random permutations correlations)
%--------------------------------------------------------------------------
disp('Computing significant correlations')
forw = [];
rev = [];

for r = 1:numSeq                                                            % For each replay
    %Forward replay
    RlessFwd = sum(RpermFwd(r,:) <= Rs(r,1));                               % Find how many permutations of the replay have less correlation with the prototype sequence than the original sequence
    RlessFwd = RlessFwd/repet*100;                                          % Turn it to percentage
%     RmoreFwd = sum(RpermFwd(r,:) >= Rs(r,1));                                  % Same for permutations with higher correlation than the original
%     RmoreFwd = RmoreFwd/repet*100;                                                % 
    %Reverse replay
    RlessRev = sum(RpermRev(r,:) <= Rs(r,2));                               % Find how many permutations of the replay have less correlation with the prototype sequence than the original sequence
    RlessRev = RlessRev/repet*100;                                          % Turn it to percentage
    if RlessFwd >= 95                                                          % If there are more than 95% permuations with correlation less than the original sequence
        forw = [forw; r];                                                   % Mark this replay as a forward replay
    end
    if RlessRev >= 95                                                      % If there are more than 95% permutations with correlation higher than the original sequence
        rev = [rev; r];                                                     % Mark this replay as a reverse replay
    end
end

disp(' ')

Forw = r_high(forw);                                               % Keep the incoming replays that are forward
Rev = r_high(rev);                                                 % and those that are reverse

if ~isempty(intersect(Forw, Rev)) 
    error('Some replays were identified as both forward AND reverse!')
end

disp([num2str(length(unique([Forw,Rev]))),' spindles contain significant replays'])
disp(' ')

Rs_f = Rs(forw,1);                                          % Keep the correlations of the incoming forward replays
Rs_r = Rs(rev,2);                                           % and the reverse ones

% plot_hist(Rs_f,hbins,'k','k',[4,4,9])                                    % Histogram plots
% plot_hist(Rs_r,hbins,'k','k',[4,4,10])

%% SET UP FINAL MATRICES
disp([num2str(length(unique([Forw,Rev]))),' spindles contain significant replays'])
disp([num2str(length(Forw)),' are forward replays of incoming cells'])
disp([num2str(length(Rev)),' are reverse replays of incoming cells'])
disp(' ')  

% plot_hist(Rs_f,hbins,'w','k',[4,4,11])                                   % Histogram plots
% plot_hist(Rs_r,hbins,'w','k',[4,4,12])

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