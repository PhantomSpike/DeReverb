function [X_ft, t, params] = cochleagram_msc_power(x, fs, dt, type, varargin)
% function [X_ft, t, params] = cochleagram(x, fs, dt, type, varargin)
%
% calculate log-frequency-spaced or cat-erb-spaced cochleagram
% of input sound.
%
% input params:
% x -- the sound
% fs -- sample rate in Hz
% dt -- desired time bin size in ms
% type -- 'log' or 'cat-erb'
% varargin -- additional parameters, currently f_min, f_max, n_f
%
% output:
% X_ft -- the cochleagram
% t -- times at which cochleagram is measured
% params -- parameters used to make the cochleagram
%
% e.g.
% 1/2-octave spacing between 500 and 1600Hz, 10ms time window:
% [X_ft, t, params] = cochleagram(rand(10000,1), 44100, 10, 'log', 500, 16000, 11)

n = 1;
c = 2.5e-3;
SAT = 1;

params.fs = fs;
params.threshold = -40;

if strcmp(type, 'log')
  if isempty(varargin)
    params.f_min = 200;
    params.f_max = 16000;
    params.n_f = 30;
  else
	[params.f_min, params.f_max, params.n_f] = varargin{:};
  end

  params.nfft_mult = 4;
  params.meltype = 'lusc';

elseif strcmp(type, 'cat-erb')
  if isempty(varargin)
    params.f_min = 1000;
    params.f_max = 32000;
    params.n_f = 23;
  else
	[params.f_min, params.f_max, params.n_f] = varargin{:};
  end
  %bank filters still have no coefficients; giving up\

  params.nfft_mult = 1;
  params.meltype = 'kusc';

end

% get actual dt (which is an integer number of samples)
dt_sec_nominal = dt/1000;
dt_bins = round(dt_sec_nominal*params.fs);
params.dt_sec = dt_bins/params.fs;

% get window, overlap sizes
t_window_bins = dt_bins * 2;
params.t_window_sec = t_window_bins/params.fs;
t_overlap_bins = t_window_bins - dt_bins;

[melbank.x, melbank.mc, melbank.na, melbank.nb] = ...
    melbankbw(params.n_f, t_window_bins*params.nfft_mult, params.fs,...
        params.f_min/params.fs, params.f_max/params.fs, params.meltype);        
if any(sum(melbank.x')==0)
	fprintf('some melbank filters have no coefficients; increasing nfft_mult\n');
  params.nfft_mult = 8;
  [melbank.x, melbank.mc, melbank.na, melbank.nb] = ...
    melbankbw(params.n_f, t_window_bins*params.nfft_mult, params.fs,...
        params.f_min/params.fs, params.f_max/params.fs, params.meltype); 
  if any(sum(melbank.x')==0) 
    error('some melbank filters still have no coefficients; giving up\n');       
  end
end

params.melbank = melbank;

params.freqs = 10.^params.melbank.mc;

[spec, freqs, t] = spectrogram(x, t_window_bins, t_overlap_bins, ...
    t_window_bins*params.nfft_mult, params.fs);

X = nan(size(params.melbank.x,1), size(spec,2));
for tt = 1:size(spec,2)
    X_ft(:,tt)=params.melbank.x*(abs(spec(params.melbank.na:params.melbank.nb,tt)).^2);
end


end