
function [f,t,cl,sc] = sp2a2_m1(sp_type,dat1,dat2,varargin);
% Type 0: [f,t,cl,sc] =  sp2a2_m1(0,dat1,dat2,samp_rate,seg_pwr,opt_str)
% Type 1: [f,t,cl,sc] =  sp2a2_m1(1,dat1,dat2,trig_times,duration,samp_rate,seg_pwr,opt_str)
% Type 2: [f,t,cl,sc] =  sp2a2_m1(2,dat1,dat2,trig_times,offset,seg_pts,samp_rate,seg_pwr,opt_str)
% Function to calculate spectra, coherence, phase & cumulant for 2 time series.
% Type 0, 1, 2 analysis - using weighted periodogram estimates.
%
% Copyright (C) 2008, David M. Halliday.
% This file is part of NeuroSpec.
%
%    NeuroSpec is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    NeuroSpec is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with NeuroSpec; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%    NeuroSpec is available at:  http://www.neurospec.org/
%    Contact:  contact@neurospec.org
%
% Input arguments
%  sp_type    Analysis type (0, 1, 2).
%  dat1       Channel 1  (input) time series vector.
%  dat2       Channel 2 (output) time series vector.
%
%  Additional arguments
%  Type 0:
%   samp_rate  Sampling rate (samples/sec).
%   seg_pwr    Segment length - specified as power of 2.
%   opt_str    Options string.
%
%  Type 1:
%   trig_times Vector: List of trigger times defining start of each data segment (in samples).
%   duration   Vector: Duration of each data segment (in samples).
%   samp_rate  Sampling rate (samples/sec).
%   seg_pwr    Segment length - specified as power of 2.
%   opt_str    Options string.
%
%  Type 2:
%   trig_times Vector: List of trigger times defining start of each data segment (in samples).
%   offset     Vector: List of offset values from trigger times to start of each data segment
%                      A separate analysis is done for each offset value.
%   seg_pts    Scalar: fixed number of data points in each segment (in samples).
%   samp_rate  Sampling rate (samples/sec).
%   seg_pwr    Segment length - specified as power of 2.
%   opt_str    Options string.
%
%
% Output arguments
%  f    column matrix with frequency domain parameter estimates.
%  t    column matrix with      time domain parameter estimates.
%  cl   single structure with scalar values related to analysis.
%  sc   optional matrix for spectral coefficients.
%
% Input Options
%  c  Include Periodogram COV test (Not valid for type 1 analysis).
%  h  Additional smoothing of spectra using hanning filter.
%  i  Invert channel reference for phase and cumulant density.
%  m  Mains/line frequency suppression.
%  n  Normalise to unit variance within each segment.
%  r  Rectification - this options requires an argument, valid arguments are:[0, 1, 2].
%                      r0:  Rectify Input  channel  (dat1)
%                      r1:  Rectify Output channel  (dat2)
%                      r2:  Rectify Both   channels (dat1 & dat2)
%  s  System Identification: Estimates gain and impulse response.
%  t  linear De-trend - this options requires an argument, valid arguments are:[0, 1, 2].
%                      t0:  De-trend Input  channel  (dat1)
%                      t1:  De-trend Output channel  (dat2)
%                      t2:  De-trend Both   channels (dat1 & dat2)
%  w  Use data window  - this options requires an argument, valid arguments are:[1, 2, 3, 4].
%                      w1: Split cosine taper applied to each segment,  5% tapered at each end.
%                      w2: Split cosine taper applied to each segment, 10% tapered at each end.
%                      w3: Split cosine taper applied to each segment, 25% tapered at each end.
%                      w4: Full  cosine taper applied to each segment, 50% tapered at each end.
%
% Options examples:
%  to set all options on (5% split taper), set opt_str='c h i m n r2 s t2 w1'
%  to rectify and de-trend both channels, set opt_str='r2 t2'
%  to rectify channel 2 and de-trend both channels, set opt_str='r1 t2'
%  to rectify channel 1 and de-trend both channels, set opt_str='r0 t2'
%
% Output parameters
%  f column 1       frequency in Hz.
%  f column 2       Log input/dat1  spectrum.
%  f column 3       Log output/dat2 spectrum.
%  f column 4       Coherence.
%  f column 5       Phase.
%  f column 6       Log of gain magnitude (with s option).
%  f column 6 or 7  Periodogram COV test  input channel (with c option).
%  f column 7 or 8  Periodogram COV test output channel (with c option).
%
%  t column 1       Lag in ms.
%  t column 2       Cumulant density.
%  t column 3       Impulse response (with s option).
%
%  cl.type          Analysis type (0, 1, 2)
%  cl.seg_size      Segment length.
%  cl.seg_tot       Number of segments.
%  cl.seg_tot_var   Effective no of segments, used to calculate confidence limits.
%  cl.samp_tot      Number of samples analysed.
%  cl.samp_rate     Sampling rate of data (samps/sec).
%  cl.dt            Time domain bin width (ms).
%  cl.df            Frequency domain bin width (Hz).
%  cl.f_c95         95% confidence limit for Log spectral estimates.
%  cl.ch_c95        95% confidence limit for coherence.
%  cl.q_c95         95% confidence limit for cumulant density.
%  cl.a_c95         95% confidence limit for impulse response      (with s option).
%  cl.col_g         column in f matrix containing gain estimate    (with s option).
%  cl.col_a         column in t matrix containing impulse response (with s option).
%  col_cova         Column containing Periodogram COV test on  input channel (with c option).
%  col_covb         Column containing Periodogram COV test on output channel (with c option).
%  cl.seg_pts       Data points in segment       (Type 2 only).
%  cl.offset        Offset value(s) for analysis (Type 2 only).
%  cl.opt_str       Copy of options string.
%  cl.what          Text label, used in plotting routines.
%
%  sc column 1      f11.
%  sc column 2      f22.
%  sc column 3      f21 (Complex).
%
% References:
% 1. Halliday D.M., Rosenberg J.R., Amjad A.M., Breeze P., Conway B.A. & Farmer S.F.
%    Progress in Biophysics and molecular Biology, 64, 237-278, 1995.
% 2. Bloomfield, P. Fourier Analysis of Time Series: An Introduction.
%    2nd edition. Wiley, New York, 2000.
% 3. Nielsen J.B., Conway B.A., Halliday D.M., Perreault M-C. and Hultborn H.
%    Journal of Physiology, 569, 291-304, 2005.
%
% Type 0: [f,t,cl,sc] =  sp2a2_m1(0,dat1,dat2,samp_rate,seg_pwr,opt_str)
% Type 1: [f,t,cl,sc] =  sp2a2_m1(1,dat1,dat2,trig_times,duration,samp_rate,seg_pwr,opt_str)
% Type 2: [f,t,cl,sc] =  sp2a2_m1(2,dat1,dat2,trig_times,offset,seg_pts,samp_rate,seg_pwr,opt_str)

% Check number of output arguments
if (nargout<3)
    error(' Not enough output arguments');
end

% Check sp_type is specified as scalar value.
if (max(size(sp_type)) ~= 1)
    error('Non scalar value for analysis type');
end

% Check for single column data
[nrow,ncol]=size(dat1);
if (ncol~=1)
    error(' Input NOT single column: dat1')
end
[nrow,ncol]=size(dat2);
if (ncol~=1)
    error(' Input NOT single column: dat2')
end

pts_tot=length(dat1);           % Determine size of data vector.
if (length(dat2)~=pts_tot)      % Check that input vectors are equal length.
    error (' Unequal length data arrays');
end

% Default options string.
opt_str='';

% Check and assign parameters according to analysis type.
switch (sp_type)
    case 0                                   % Type 0.
        if (nargin<5)                          % Check numbers of arguments.
            error(' Type 0 - Not enough input arguments');
        end
        samp_rate=varargin{1};   % Assign parameters.
        seg_pwr=varargin{2};
        if (nargin>5)
            opt_str=varargin{3};
        end
        if (max(size(samp_rate)) ~= 1)
            error(' Type 0 - Non scalar value for: samp_rate');
        end
        if (max(size(seg_pwr)) ~= 1)
            error(' Type 0 - Non scalar value for: seg_pwr');
        end
        seg_size=2^seg_pwr;                    % DFT segment length (S).
        seg_tot=fix(pts_tot/seg_size);         % Number of complete segments (L).
        samp_tot=seg_tot*seg_size;             % Number of samples to analyse: R=LS.
        seg_samp_start=(1:seg_size:samp_tot)'; % Start times for each segment.
        seg_samp_no(1:seg_tot,1)=seg_size;     % Fixed number of data points per segment, T = S.
        offset=0;                              % No offset values for type 0 analysis.
    case 1              % Type 1.
        if (nargin<7)     % Check numbers of arguments.
            error(' Type 1 - Not enough input arguments');
        end
        trig_times=varargin{1}; % Assign parameters.
        duration=varargin{2};
        samp_rate=varargin{3};
        seg_pwr=varargin{4};
        if (nargin>7)
            opt_str=varargin{5};
        end
        [nrow,ncol]=size(trig_times);
        if (ncol~=1)
            error(' Input NOT single column: trig_times')
        end
        [nrow,ncol]=size(duration);
        if (ncol~=1)
            error(' Input NOT single column: duration')
        end
        if (min(trig_times)<1)
            error(' Type 1 - Negative or zero trig_times')
        end
        if (min(duration)<1)
            error(' Type 1 - Negative or zero duration')
        end
        if (length(trig_times) ~= length(duration))
            error(' Type 1 - Unequal numbers of trig_times, duration')
        end
        if (max(trig_times>pts_tot))
            error(' Type 1 - trig_times exceed data length')
        end
        if (max(trig_times+duration)-1>pts_tot)
            error(' Type 1 - trig_times+duration exceed data length')
        end
        if (max(size(samp_rate)) ~= 1)
            error(' Type 1 - Non scalar value for: samp_rate');
        end
        if (max(size(seg_pwr)) ~= 1)
            error(' Type 1 - Non scalar value for: seg_pwr');
        end
        seg_size=2^seg_pwr; % DFT segment length (S).
        seg_min=0.05;       % Define minimum percentage of data points allowed in each segment.
        seg_samp_min=round(seg_min*seg_size); % Convert to minimum number of data points.
        seg_tot=0;          % Counts number of segments in analysis.
        samp_tot=0;         % Counts number of samples  in analysis.
        offset=0;           % No offset values for type 1 analysis.
        for ind=1:length(trig_times)  % Loop through all trigger times.
            seg_start_offset=0;         % Offset within each block of data.
            while ((duration(ind)-seg_start_offset)>=seg_samp_min)
                seg_tot=seg_tot+1;        % Additional segment in this block.
                seg_samp_start(seg_tot,1)=trig_times(ind)+seg_start_offset;  % Start of segment.
                if ((duration(ind)-seg_start_offset)>seg_size)
                    seg_samp_no(seg_tot,1)=seg_size;                        % Data for complete segment.
                else
                    seg_samp_no(seg_tot,1)=duration(ind)-seg_start_offset;  % Data for part segment only.
                end
                samp_tot=samp_tot+seg_samp_no(seg_tot,1);                 % Update number of samples.
                seg_start_offset=seg_start_offset+seg_samp_no(seg_tot,1); % Update start offset in block.
            end
        end
    case 2                            % Type 2.
        if (nargin<8)                   % Check numbers of arguments.
            error(' Type 2 - Not enough input arguments');
        end
        trig_times=varargin{1}; % Assign parameters.
        offset=varargin{2};
        seg_pts=varargin{3};
        samp_rate=varargin{4};
        seg_pwr=varargin{5};
        if (nargin>8)
            opt_str=varargin{6};
        end
        if (max(size(samp_rate)) ~= 1)
            error(' Type 2 - Non scalar value for: samp_rate');
        end
        if (max(size(seg_pwr)) ~= 1)
            error(' Type 2 - Non scalar value for: seg_pwr');
        end
        seg_tot=length(trig_times);       % Number of segments, L, from no of triggers.
        samp_tot=seg_tot*seg_pts;         % Number of samples to analyse.
        seg_samp_start=trig_times;        % Start times for each segment, from trigger times.
        seg_samp_no(1:seg_tot,1)=seg_pts; % Fixed number of samples per segment, T=seg_pts.
        seg_size=2^seg_pwr;               % DFT segment length (S).
        [nrow,ncol]=size(trig_times);
        if (ncol~=1)
            error(' Type 2 - Input NOT single column: trig_times')
        end
        [nrow,ncol]=size(offset);
        if (ncol~=1)
            error(' Type 2 - Input NOT single column: offset')
        end
        if (min(trig_times)<1)
            error(' Type 2 - Negative or zero trig_times')
        end
        if (max(trig_times)>pts_tot)
            error(' Type 2 - trig_times exceed data length')
        end
        if (seg_pts<1)
            error(' Type 2 - Negative or zero seg_pts')
        end
        if (seg_pts > seg_size)
            error(' Type 2 - seg_pts exceeds DFT segment size')
        end
        if (max(trig_times)+seg_pts-1>pts_tot)
            error(' Type 2 - trig_times+seg_pts exceed data length')
        end
        if (min(trig_times)+min(offset)<1)
            error(' Type 2 - Negative or zero trig_times+offset')
        end
        if (max(trig_times)+max(offset)+seg_pts-1>pts_tot)
            error(' Type 2 - trig_times+offset+seg_pts exceed data length')
        end
    otherwise
        error ([ 'Type error -- ',num2str(sp_type)]);  % Illegal type.
end

% Create data matrices, S rows, L columns.
rd1=zeros(seg_size,seg_tot);
rd2=zeros(seg_size,seg_tot);

% Process options.
flags.sp_type=sp_type;   % Type passed to periodogram subroutine.
if (sp_type == 2)
    flags.seg_pts=seg_pts; % Type 2: seg_pts used to define start frequency in f_out.
end
flags.line=0;            % Set defaults - options off.
flags.inv=0;
flags.han=0;
flags.wind=0;
flags.gain=0;
flags.cov=0;
norm_chan=0;
trend_chan_1=0;
trend_chan_2=0;
rect_chan_1=0;
rect_chan_2=0;
options=deblank(opt_str);
while (any(options))              % Parse individual options from string.
    [opt,options]=strtok(options);
    optarg=opt(2:length(opt));      % Determine option argument.
    switch (opt(1))
        case 'c'             % Include Periodogram COV test.
            flags.cov=1;
            if (sp_type==1)
                disp('Warning - Periodogram COV test not compatible with Type 1 analysis. Option disabled')
                flags.cov=0;
            end
        case 'h'             % Use additional hanning smoothing on spectra.
            flags.han=1;
        case 'i'             % Channel reference invert option.
            flags.inv=1;
        case 'm'             % Mains/line frequency suppression option.
            flags.line=1;
        case 'n'             % Normalisation to unit variance within each segment.
            norm_chan=1;
        case 'r'             % Rectification option.
            i=str2num(optarg);
            if (i<0 | i>2)
                error(['error in option argument -- r',optarg]);
            end
            if (i~=1)
                rect_chan_1=1;     % Rectify ch 1.
            end
            if (i>=1)
                rect_chan_2=1;     % Rectify ch 2.
            end
        case 's'             % System Identification option (gain and impulse response).
            flags.gain=1;
        case 't'             % Linear de-trend option.
            i=str2num(optarg);
            if (i<0 | i>2)
                error(['error in option argument -- t',optarg]);
            end
            if (i~=1)
                trend_chan_1=1;    % De-trend ch 1.
            end
            if (i>=1)
                trend_chan_2=1;    % De-trend ch 2.
            end
        case 'w'             % Data window (taper) option
            i=str2num(optarg);
            if (i<1 | i>4)
                error(['error in option argument -- w',optarg]);
            end
            flags.wind=i;
        otherwise
            error (['Illegal option -- ',opt]);  % Illegal option.
    end
end

if (trend_chan_1 | trend_chan_2)       % Index for fitting data with polynomial.
    trend_x=(1:seg_size)';
end

for ind_off=1:length(offset)            % Loop across all offset values.
    for ind=1:seg_tot                        % Loop across columns/segments.
        seg_pts=seg_samp_no(ind);                       % No of data points in segment.
        seg_start=seg_samp_start(ind)+offset(ind_off);  % Start sample in data vector.
        seg_stop=seg_start+seg_pts-1;                   % Stop  sample in data vector.
        dat_seg1=dat1(seg_start:seg_stop);              % Extract segment from dat1.
        dat_seg2=dat2(seg_start:seg_stop);              % Extract segment from dat2.
        md1=mean(dat_seg1);                             % Mean of segment from dat1.
        md2=mean(dat_seg2);                             % Mean of segment from dat2.
        
        rd1(1:seg_pts,ind)=dat_seg1-md1;        % Subtract mean from ch 1.
        rd2(1:seg_pts,ind)=dat_seg2-md2;        % Subtract mean from ch 2.
        if rect_chan_1
            rd1(1:seg_pts,ind)=abs(rd1(1:seg_pts,ind));            % Rectification of ch 1 (Full wave).
        end
        if rect_chan_2
            rd2(1:seg_pts,ind)=abs(rd2(1:seg_pts,ind));            % Rectification of ch 2 (Full wave).
        end
        if trend_chan_1                               % Linear trend removal.
            p=polyfit(trend_x(1:seg_pts,1),rd1(1:seg_pts,ind),1);                 % Fit 1st order polynomial.
            rd1(1:seg_pts,ind)=rd1(1:seg_pts,ind)-p(1)*trend_x(1:seg_pts,1)-p(2); % Subtract from ch 1.
        end
        if trend_chan_2                               % Linear trend removal.
            p=polyfit(trend_x(1:seg_pts,1),rd2(1:seg_pts,ind),1);                 % Fit 1st order polynomial.
            rd2(1:seg_pts,ind)=rd2(1:seg_pts,ind)-p(1)*trend_x(1:seg_pts,1)-p(2); % Subtract from ch 2.
        end
        if norm_chan
            rd1(1:seg_pts,ind)=rd1(1:seg_pts,ind)/std(rd1(1:seg_pts,ind));
            rd2(1:seg_pts,ind)=rd2(1:seg_pts,ind)/std(rd2(1:seg_pts,ind));
        end
    end
    
    % Call sp2_fn2a() Weighted periodogram based spectral estimation routine.
    if (nargout>3)
        [f(:,:,ind_off),t(:,:,ind_off),cl(ind_off),sc(:,:,ind_off)] = sp2_fn2a(rd1,rd2,seg_samp_no,samp_rate,flags);
    else
        [f(:,:,ind_off),t(:,:,ind_off),cl(ind_off)]                 = sp2_fn2a(rd1,rd2,seg_samp_no,samp_rate,flags);
    end
    
    % Display No of segments & resolution
    if (sp_type ==2)
        disp(['Segments: ',num2str(seg_tot),', Segment length: ',num2str(seg_size/samp_rate),...
            ' sec,  Offset: ',num2str(offset(ind_off)),' pts,   Resolution: ',num2str(cl(ind_off).df),' Hz.']);
    else
        disp(['Segments: ',num2str(seg_tot),', Segment length: ',num2str(seg_size/samp_rate),...
            ' sec,  Resolution: ',num2str(cl(ind_off).df),' Hz.']);
    end
end

% Set additional data and type dependent elements in cl structure.
for ind_off=1:length(offset)
    switch (sp_type)
        case 0                                % Type 0.
        case 1                                % Type 1.
        case 2                                % Type 2.
            cl(ind_off).seg_pts=seg_pts;        % Data points in segment.
            cl(ind_off).offset=offset(ind_off); % Offset.
    end
    cl(ind_off).N1=0;                       % N1, No of events in ch 1. (zero for TS data)
    cl(ind_off).N2=0;                       % N2, No of events in ch 2.          "
    cl(ind_off).P1=0;                       % P1, mean intensity ch1.            "
    cl(ind_off).P2=0;                       % P2, mean intensity ch2.            "
    cl(ind_off).opt_str=opt_str;            % Copy of options string.
    cl(ind_off).what='';                    % Field for plot label.
end

function [f,t,cl,sc] = sp2_fn2a(d1,d2,seg_samp_no,samp_rate,flags);
% function [f,t,cl,sc] = sp2_fn2a(d1,d2,seg_samp_no,samp_rate,flags);
% Function with core routines for periodogram based spectral estimates.
% Implements a weighted periodogram analysis with:
%  Variable number of data points in each segment, T
%  Fixed DFT segment length, S.
% see NOTE below about calling this routine.
%
% Copyright (C) 2008, David M. Halliday.
% This file is part of NeuroSpec.
%
%    NeuroSpec is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    NeuroSpec is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with NeuroSpec; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%    NeuroSpec is available at:  http://www.neurospec.org/
%    Contact:  contact@neurospec.org
%
%  Inputs are two matrices containing pre processed time series or point process data.
%  Matrices have L columns with S rows.
%   L = Number of segments: seg_tot.
%   S = DFT segment length: seg_size.
%
% Input arguments
%  d1          Channel 1 data matrix.
%  d2          Channel 2 data matrix.
%  seg_samp_no Vector with number of data points (T) in each segment.
%  samp_rate   Sampling rate (samples/sec)
%  flags       Structure with flags to control processing options.
%              Flags supported:
%                line: Apply simple filter to suppress mains/line frequency (0:No; 1:Yes).
%                inv:  Invert channel reference for phase & cumulant        (0:No; 1:Yes).
%                han:  Additional smoothing to spectra using hanning filter (0:No; 1:Yes).
%                wind: Apply data window (split cosine taper) to raw data.
%                        (0:No; 1,2,3,...: Type).  1:10%; 2:20%; 3:50%; 4:100% taper.
%                gain: Estimate additional system identification parameters (0:No; 1:Yes).
%                        Log{|Gain|} included in f matrix.
%                        Impulse response included in t matrix.
%                cov:  Calculates Periodogram COV test for both channels.
%
% Output arguments
%  f    column matrix with frequency domain parameter estimates.
%  t    column matrix with      time domain parameter estimates.
%  cl   single structure with scalar values related to analysis.
%  sc   optional matrix for spectral coefficients.
%
% Output parameters
%  f column 1       frequency in Hz.
%  f column 2       Log  input/d1 spectrum.
%  f column 3       Log output/d2 spectrum.
%  f column 4       Coherence.
%  f column 5       Phase.
%  f column 6       Log of gain magnitude (with s option).
%  f column 6 or 7  Periodogram COV test  input channel (with c option).
%  f column 7 or 8  Periodogram COV test output channel (with c option).
%
%  t column 1       Lag in ms.
%  t column 2       Cumulant density.
%  t column 3       Impulse response (with s option).
%
%  cl.type          Analysis type (0, 1, 2)
%  cl.seg_size      Segment length.
%  cl.seg_tot       Number of segments.
%  cl.seg_tot_var   Effective no of segments, used to calculate confidence limits.
%  cl.samp_tot      Number of samples analysed.
%  cl.samp_rate     Sampling rate of data (samps/sec).
%  cl.dt            Time domain bin width (ms).
%  cl.df            Frequency domain bin width (Hz).
%  cl.f_c95         95% confidence limit for Log spectral estimates.
%  cl.ch_c95        95% confidence limit for coherence.
%  cl.q_c95         95% confidence limit for cumulant density.
%  cl.a_c95         95% Confidence limits for impulse response (with s option).
%  cl.col_g         Column containing log gain (with s option).
%  cl.col_a         Column containing impulse response (with s option).
%  col_cova         Column containing Periodogram COV test on  input channel (with c option).
%  col_covb         Column containing Periodogram COV test on output channel (with c option).
%
%  sc column 1      f11.
%  sc column 2      f22.
%  sc column 3      f21 (Complex).
%
% References:
% 1. Halliday D.M., Rosenberg J.R., Amjad A.M., Breeze P., Conway B.A. & Farmer S.F.
%    Progress in Biophysics and molecular Biology, 64, 237-278, 1995.
% 2. Bloomfield, P. Fourier Analysis of Time Series: An Introduction.
%    2nd edition. Wiley, New York, 2000.
% 3. Nielsen J.B., Conway B.A., Halliday D.M., Perreault M-C. and Hultborn H.
%    Journal of Physiology, 569, 291-304, 2005.
%
% function [f,t,cl,sc] = sp2_fn2a(d1,d2,seg_samp_no,samp_rate,flags);
%
% NOTE: This routine is not intended to support analysis of raw data.
% It is intended as a support routine for the 2 channel spectral analysis functions:
% sp2_m1.m, sp2a_m1.m, sp2a2_m1.m. Refer to these functions for further details.

% PBMB     refers to above Progress in Biophysics article.
% JPHYSIOL refers to above  Journal of Physiology article.

[seg_size,seg_tot]=size(d1); % Calculate seg_size & seg_tot from data matrices.
samp_tot=sum(seg_samp_no);   % Total no of samples, Sum(T).

% Check if matrix of periodogram and cross periodogram coefficients required
flags.matrix=0;   % Default - No
if (flags.cov)    % Required in Periodogram COV test.
    flags.matrix=1;
end

% Apply data window, if required.
if (flags.wind)
    switch (flags.wind)
        case 1       % 10%  split cosine taper applied  (5% at each end).
            p=0.1;
            u2=0.9377; % U2 factor (Bloomfield, P181).
        case 2       % 20%  split cosine taper applied (10% at each end).
            p=0.2;
            u2=0.8751;
        case 3       % 50%  split cosine taper applied (25% at each end).
            p=0.5;
            u2=0.6877;
        case 4       % 100%  full cosine taper applied (50% at each end).
            p=1;
            u2=0.3752;
    end
    
    % Calculate length of default taper.
    % Type 0, 1: data segment is length seg_size.
    % Type    2: data segment is length flags.seg_pts.
    if (flags.sp_type==2)
        tap_pts_tot=flags.seg_pts;
    else
        tap_pts_tot=seg_size;
    end
    
    % Generate lower taper section for default taper.
    % No of points to taper at each end.
    tap_pts=round(tap_pts_tot*p/2);
    % Generate lower half of cosine over this number of points.
    taper(1:tap_pts,1)=0.5-0.5*cos(pi*(1:tap_pts)'/tap_pts);  % Bloomfield, (6.9).
    % Indices for upper data section to taper, in reverse order.
    tap_ind_up=tap_pts_tot:-1:tap_pts_tot-tap_pts+1;
    
    % Taper each segment
    for ind=1:seg_tot
        if (seg_samp_no(ind)==tap_pts_tot)
            % Use default taper.
            d1(1:tap_pts,ind) =d1(1:tap_pts,ind) .*taper; % Lower section.
            d1(tap_ind_up,ind)=d1(tap_ind_up,ind).*taper; % Upper section.
            d2(1:tap_pts,ind) =d2(1:tap_pts,ind) .*taper; % Lower section.
            d2(tap_ind_up,ind)=d2(tap_ind_up,ind).*taper; % Upper section.
            tap_pts_no(ind,1)=tap_pts;
        else
            % Have smaller number of data points - generate custom taper.
            % No of points to taper at each end.
            tap_pts_1=round(seg_samp_no(ind)*p/2);
            % Generate lower half of cosine over this number of points.
            taper_1=[];  % Clear any previous taper.
            taper_1(1:tap_pts_1,1)=0.5-0.5*cos(pi*(1:tap_pts_1)'/tap_pts_1); % Bloomfield, (6.9).
            % Indices for upper data section to taper in reverse order.
            tap_ind_up_1=seg_samp_no(ind):-1:seg_samp_no(ind)-tap_pts_1+1;
            d1(1:tap_pts_1,ind) =d1(1:tap_pts_1,ind) .*taper_1; % Lower section.
            d1(tap_ind_up_1,ind)=d1(tap_ind_up_1,ind).*taper_1; % Upper section.
            d2(1:tap_pts_1,ind) =d2(1:tap_pts_1,ind) .*taper_1; % Lower section.
            d2(tap_ind_up_1,ind)=d2(tap_ind_up_1,ind).*taper_1; % Upper section.
            tap_pts_no(ind,1)=tap_pts_1;
        end
    end
    
    % Normalisation factor for DFT to take account of tapering.
    u2_norm=sqrt(u2);            % Bloomfield, (9.11).
    fd1=fft(d1)/u2_norm;         % Take DFT across columns/segments ch 1, PBMB (4.1)/(4.2).
    fd2=fft(d2)/u2_norm;         % Take DFT across columns/segments ch 2, PBMB (4.1)/(4.2).
    
else
    % No tapering
    fd1=fft(d1);                 % Take DFT across columns/segments ch 1, PBMB (4.1)/(4.2).
    fd2=fft(d2);                 % Take DFT across columns/segments ch 2, PBMB (4.1)/(4.2).
end
t_fac=2*pi*samp_tot;           % Normalization for weighted periodogram spectral estimates.
% JPHYSIOL (5), (6).

% Constuct 2D matrix of periodograms for individual segments.
if flags.matrix
    % Definition of indices to include - DC to f_n.
    per_ind=(1:seg_size/2+1)';
    for ind=1:seg_tot
        if flags.inv % Channel reference invert:  Channel 2 is input, Channel 1 is output.
            per_2d(:,1,ind)=abs(fd2(per_ind,ind).*fd2(per_ind,ind)/t_fac);   % Periodogram ch1
            per_2d(:,2,ind)=abs(fd1(per_ind,ind).*fd1(per_ind,ind)/t_fac);   % Periodogram ch2
            per_2d(:,3,ind)=fd1(per_ind,ind).*conj(fd2(per_ind,ind)/t_fac);  % Cross periodogram (complex)
        else         % Channels in default order: Channel 1 is input, Channel 2 is output.
            per_2d(:,1,ind)=abs(fd1(per_ind,ind).*fd1(per_ind,ind)/t_fac);   % Periodogram ch1
            per_2d(:,2,ind)=abs(fd2(per_ind,ind).*fd2(per_ind,ind)/t_fac);   % Periodogram ch2
            per_2d(:,3,ind)=fd2(per_ind,ind).*conj(fd1(per_ind,ind)/t_fac);  % Cross periodogram (complex)
        end
    end
end

% Construct spectra based on average across segments
if flags.inv
    % Channel reference invert: Channel 2 is input, Channel 1 is output.
    f11=sum(abs(fd2.*fd2)/t_fac,2);   % Spectrum 1, PBMB (5.2), Mag squared for  input auto spectra.
    f22=sum(abs(fd1.*fd1)/t_fac,2);   % Spectrum 2, PBMB (5.2), Mag squared for output auto spectra.
    f21=sum(fd1.*conj(fd2)/t_fac,2);  % Cross spectrum (complex valued), PBMB (5.2).
else
    % Default order - Channel 1 is input, Channel 2 is output.
    f11=sum(abs(fd1.*fd1)/t_fac,2);   % Spectrum 1, PBMB (5.2), Mag squared for  input auto spectra.
    f22=sum(abs(fd2.*fd2)/t_fac,2);   % Spectrum 2, PBMB (5.2), Mag squared for output auto spectra.
    f21=sum(fd2.*conj(fd1)/t_fac,2);  % Cross spectrum (complex valued), PBMB (5.2).
end

deltaf=samp_rate/seg_size;        % Spacing of Fourier frequencies in Hz.

% Set line frequency in Hz.
line_freq=50;
% line_freq=60;  % Uncomment this line to set the line frequency to 60 Hz.

% Suppression of mains/line frequency - smooth out using adjacent values.
if flags.line
    line_ind=round(line_freq/deltaf)+1;        % NB Index 1 is DC.
    f11(line_ind)=0.5*(f11(line_ind-2)+f11(line_ind+2));    % Spectrum ch 1.
    f11(line_ind-1)=0.5*(f11(line_ind-2)+f11(line_ind-3));
    f11(line_ind+1)=0.5*(f11(line_ind+2)+f11(line_ind+3));
    f22(line_ind)=0.5*(f22(line_ind-2)+f22(line_ind+2));    % Spectrum ch 2.
    f22(line_ind-1)=0.5*(f22(line_ind-2)+f22(line_ind-3));
    f22(line_ind+1)=0.5*(f22(line_ind+2)+f22(line_ind+3));
    f21(line_ind)=0.5*(f21(line_ind-2)+f21(line_ind+2));    % Cross spectrum.
    f21(line_ind-1)=0.5*(f21(line_ind-2)+f21(line_ind-3));
    f21(line_ind+1)=0.5*(f21(line_ind+2)+f21(line_ind+3));
    % Smooth elements in upper hermetian section of cross spectral estimate.
    % This data used in ifft() to generate cumulant. Data is complex conjugate.
    f21(seg_size-line_ind+2)=conj(f21(line_ind));
    f21(seg_size-line_ind+3)=conj(f21(line_ind-1));
    f21(seg_size-line_ind+1)=conj(f21(line_ind+1));
end

% Estimate cumulant density using inverse DFT of cross spectrum.
deltat=1000.0/samp_rate;    % dt in msec.

cov=ifft(f21);              % Inverse DFT.

% Construct output time domain matrix t.
% Column 1 of t matrix is time in msec. Range (-S/2)*dt to (S/2-1)*dt.
t(:,1)=((1:seg_size)'-seg_size/2-1)*deltat;

% Column 2 of t matrix is cumulant, shifted by S/2 so that time zero is in centre.
% 2pi/S factor is 2*pi, since ifft routine includes 1/S term.
t([seg_size/2+1:seg_size,1:seg_size/2],2)=real(cov(1:seg_size))*2*pi; % PBMB (5.9).

% Estimate variance of cumulant density estimate.
var_fac=4*pi*pi/(seg_size*samp_tot);                           % Factor (2pi/S)(2pi/R).
q_var=var_fac*2*sum(f11(1:seg_size/2+1).*f22(1:seg_size/2+1)); % PBMB (6.10).

% Optional impulse response calculation.
if flags.gain
    % Estimate impulse response using inverse DFT.
    f21_gain(1,1)=0;    % Disregard DC value, may well be zero.
    % Calculate gain, PBMB (11.2), f11 is input spectrum.
    f21_gain(2:seg_size,1)=f21(2:seg_size)./f11(2:seg_size);
    
    a=ifft(f21_gain);
    
    % Append impulse response to t matrix, shifted by T/2 so that time zero is in centre.
    % 1/T factor not included, since ifft routine includes 1/T term.
    col_a=size(t,2)+1;
    t([seg_size/2+1:seg_size,1:seg_size/2],col_a)=real(a(1:seg_size)); % PBMB (11.5).
    
    % Calculate variance of impulse response estimate.
    a_var_fac=1/(seg_size*samp_tot);                                 % Factor (1/T)(1/R).
    a_var=a_var_fac*2*sum(f22(2:seg_size/2+1)./f11(2:seg_size/2+1)); % PBMB (12.4).
end

% Reduce arrays to lower hermitian section only
f11=f11(1:seg_size/2+1);
f22=f22(1:seg_size/2+1);
f21=f21(1:seg_size/2+1);

% Apply smoothing using hanning filter.
var_smooth=1.0;     % Default - no smoothing correction.
if (flags.han==1)
    f11=han(f11);
    f22=han(f22);
    f21=han(f21);
    if flags.gain
        % Recalculate gain from smoothed spectra. PBMB (11.2), f11 is input spectrum.
        % Calculated only over lower hermitian section.
        f21_gain(2:seg_size/2+1,1)=f21(2:seg_size/2+1)./f11(2:seg_size/2+1);
    end
    var_smooth=0.375; % Correction for variance of smoothed estimates, Bloomfield (9.7).
end

% Construct output spectral matrix f.
f_index=(2:seg_size/2+1)';               % Indexing for output, DC component not output.
% Set minimum frequency output value for type two analysis, which takes account of zero padding.
if (flags.sp_type ==2)
    f_res=samp_rate/flags.seg_pts;         % Frequency resolution for type 2 analysis.
    % Calculate starting index for f matrix, with resolution deltaf.
    % Gives first frequency ordinate above f_res, +1 to allow for DC vlaue in first bin.
    f_index_start=ceil(f_res/deltaf)+1;
    f_index=(f_index_start:seg_size/2+1)'; % Indexing for type 2 analysis.
end
f(:,1)=(f_index-1)*deltaf;    % Column 1 - frequencies in Hz.
f(:,2)=log10(f11(f_index));   % Column 2 - Log spectrum ch 1.
f(:,3)=log10(f22(f_index));   % Column 3 - Log spectrum ch 2.
% Column 4 - Coherence, PBMB (5.5).
f(:,4)=abs(f21(f_index)).*abs(f21(f_index))./(f11(f_index).*f22(f_index));
f(:,5)=angle(f21(f_index));   % Column 5 - Phase, PBMB (5.7).

% Optional gain calculation.
if flags.gain
    % Append log magnitude of gain to f matrix, PBMB (11.3)
    col_g=size(f,2)+1;
    f(:,col_g)=log10(abs(f21_gain(f_index)));
end

% Optional periodogram COV estimates
if flags.cov
    col_cova=size(f,2)+1;
    col_covb=size(f,2)+2;
    f(:,col_cova)=std(per_2d(f_index,1,:),0,3)./mean(per_2d(f_index,1,:),3);  % COV of periodogram Ch 1
    f(:,col_covb)=std(per_2d(f_index,2,:),0,3)./mean(per_2d(f_index,2,:),3);  % COV of periodogram Ch 2
end

% Construct optional output spectral matrix sc, DC value is output.
if (nargout>3)
    sc(:,1)=f11(1:seg_size/2+1);  % Column 1 - f11
    sc(:,2)=f22(1:seg_size/2+1);  % Column 2 - f22
    sc(:,3)=f21(1:seg_size/2+1);  % Column 3 - f21 (Complex valued)
end

% Calculate effective number of segments, L', used to approximate variance of estimates.
seg_tot_var=1/sum((seg_samp_no(:,1)/samp_tot).*(seg_samp_no(:,1)/samp_tot)); % JPHYSIOL (8).
% Correct for smoothing.
seg_tot_var=seg_tot_var/var_smooth;

% Construct cl structure, confidence limits for parameter estimates.
cl.type=flags.sp_type;       % Analysis type.
cl.seg_size=seg_size;        % S.
cl.seg_tot=seg_tot;          % L.
cl.seg_tot_var=seg_tot_var;  % Effective no of segments (L').
cl.samp_tot=samp_tot;        % R.
cl.samp_rate=samp_rate;      % Sampling rate.
cl.dt=deltat;                % Delta t.
cl.df=deltaf;                % Delta f.
cl.f_c95=0.8512*sqrt(1/seg_tot_var);  % 95% Confidence limit for spectral estimates, PBMB (6.2).
% N.B. Confidence interval for log plot of spectra is TWICE this value.
cl.ch_c95=1-0.05^(1/(seg_tot_var-1)); % 95% Confidence limit for coherence, PBMB (6.6).
cl.q_c95=1.96*sqrt(q_var);            % 95% Confidence limits for cumulant, PBMB (6.11).

% Extra values for optional gain and impulse response estimates.
if flags.gain
    cl.a_c95=1.96*sqrt(a_var); % 95% Confidence limits for impulse response, PBMB (12.5).
    cl.col_g=col_g;            % Column in f matrix containing gain.
    cl.col_a=col_a;            % Column in t matrix containing impulse response.
end

if flags.cov
    cl.cov_c95=1.96*sqrt(3/(2*seg_tot));
    cl.col_cova=col_cova;
    cl.col_covb=col_covb;
end

return

%------------------------------------------------------------------------------
function [dat_out] = han(dat_in);
% function [dat_out] = han(dat_in);
%
% Function to smooth data vector using Hanning filter with coefficients: (1/4, 1/2, 1/4).
% End points smoothed using two point filter: (1/2, 1/2).

a=1;
b=[0.25;0.5;0.25];
dat_pts=length(dat_in);

% Smooth using function: filter.
dat_out=filter(b,a,dat_in);

% Shift one place to cancel out one sample delay in filter function.
dat_out=[dat_out(2:dat_pts);0];

% Smooth end points using two point filter: (1/2, 1/2).
dat_out(1)=0.5*(dat_in(1)+dat_in(2));
dat_out(dat_pts)=0.5*(dat_in(dat_pts)+dat_in(dat_pts-1));
return

