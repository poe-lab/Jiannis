function [sigFilt,x1,x2,y1,y2] = bandpass_c(signal, Fs, freq1, freq2, filt_flag)

%--------------------------------------------------------------------------
% Bandstops a signal between freq1 and freq2 by removing the corresponding
% Fourrier components.
%
% INPUT VARIABLES
% signal: Time-series to be analyzed
% freq1/2: Frequencies for which power is set to zero
% Fs: Sampling rate (Hz)
% OUTPUT VARIABLES
% sigFilt: Lowpass-filtered time-series (same sampling rate)
% x1/2: Locations in discrete frequency spectrum closest to "freq1/2"
% y1/2: Exact frequencies at x1/2
%
% Developed by Costas Anastassiou. Modified by Jiannis Taxidis
% November 2012, Caltech, USA
%--------------------------------------------------------------------------

if nargin == 4
    filt_flag = 'pass';
    disp('The filter was set to bandpass by default')
end
if ~any([strcmp(filt_flag,'pass'),strcmp(filt_flag,'stop')])
    error('The filter was not defined')
end

L = length(signal);                                                         % Length of signal
NFFT = 2^nextpow2(L);                                                       % Next power of 2 from length of y
Y = fft(signal,NFFT)/L;                                                     % Apply Fourier transform
f = Fs/2*linspace(0,1,NFFT/2+1);                                            % Frequency spectrum

% Locations in discrete frequency spectrum closest to "freq1/2"
x1 = find(f>freq1,1,'first');
x2 = find(f<freq2,1,'last');

% Frequencies at x1/2
y1 = f(x1);
y2 = f(x2);

% Set power between frequencies "freq1/2" equal to zero
% If freq2 was set to zero, then set all frequencies above freq1 to zero
% (lowpass)
if strcmp(filt_flag,'pass')
    if freq2 == 0
        Y(1:x1) = 0;
    elseif freq1 == 0
        Y(x2:end) = 0;
    else
        Y(1:x1) = 0;
        Y(x2:end) = 0;
    end
elseif strcmp(filt_flag,'stop')
    if freq2 == 0
        Y(x1:end) = 0;
    elseif freq1 == 0
        Y(1:x2) = 0;
    else
        Y(x1:x2) = 0;
    end
end

% Filtered time-series
dummy = L*ifft(Y,'symmetric');
sigFilt = dummy(1:L);

% Frequency analysis of the filtered signal
clear Y f
% Y = fft(sigFilt,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1);
% ampl = 2*abs(Y(1:NFFT/2+1));


% Figures
% if(nargout == 0)
%     lcolor = 'k';
%     figure; set(gcf,'color','w');
%     subplot(2,1,1);
%     plot(signal(1:L),'-c','linewidth',2); hold on;
%     plot(sigFilt, '-', 'color', lcolor, 'linewidth',2);
%     subplot(2,1,2);
%     semilogy(f, ampl, '-', 'color', lcolor, 'linewidth',2); hold on;
% end

