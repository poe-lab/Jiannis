function [minl, maxl] = plot_spikes(spikesin,cellsin,lfp)

lin = length(cellsin);

maxl = max(lfp) + 0.1*lin;
minl = min(lfp);

for c = cellsin                                                       % For each incoming cell
    sin = ones(size(spikesin{c}));
    line([spikesin{c} spikesin{c}]' , (sin * (maxl-([lin+1, lin] - c)*0.1))','Color','r'); % Draw a line for every spike around the event (red, above)
end
