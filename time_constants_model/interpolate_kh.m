function [int_k_h,int_t_ms] = interpolate_kh(k_h,t_ms,upsamp_factor,method)
%function [int_wf,int_time_ms] = interpolate_wf(wf,wf_time,upsamp_factor,method)
%This function does interpolation on spike waveforms to obtain a more
%continuous version of the original waveform. This aids in the computaton
%of clustering dimensions parameters used later
%>>INPUT>>
%wf - The original waveform values
%wf_time - The time points where the wf values were sampled
%upsamp_factor - Upsampling factor. How many times to increase the number
%of interpolated points
%method - Which method to use for the interpolation. Possible options are:
%'linear' - Linear interpolation. The interpolated value at a query point is based on linear interpolation of the values at neighboring grid points in each respective dimension. This is the default interpolation method.
%'nearest' - Nearest neighbor interpolation. The interpolated value at a query point is the value at the nearest sample grid point.
%'next' - Next neighbor interpolation. The interpolated value at a query point is the value at the next sample grid point.
%'previous' - Previous neighbor interpolation. The interpolated value at a query point is the value at the previous sample grid point.
%'pchip' - Shape-preserving piecewise cubic interpolation. The interpolated value at a query point is based on a shape-preserving piecewise cubic interpolation of the values at neighboring grid points.
%'v5cubic' - Cubic convolution used in MATLAB® 5.
%'makima' - Modified Akima cubic Hermite interpolation. The interpolated value at a query point is based on a piecewise function of polynomials with degree at most three. The Akima formula is modified to avoid overshoots.
%'spline' - Spline interpolation using not-a-knot end conditions. The interpolated value at a query point is based on a cubic interpolation of the values at neighboring grid points in each respective dimension.
%<<OUTPUT<<
%int_wf - The interpolated waveform
%int_time_ms - The new time vector with the novel itnerpolation points 

if ~exist('upsamp_factor','var')
    upsamp_factor = 10;
end
if ~exist('method', 'var')
    method = 'makima';
end

d_ms = t_ms(2)-t_ms(1); % Find the difference between two adjacent samples 
int_t_ms = [t_ms(1):d_ms/upsamp_factor:t_ms(end)]; %Define new time points where you want to interpolate
k_h = k_h'; %Transpose input because of interp1 requirement 
int_k_h = interp1(t_ms,k_h,int_t_ms,method);
int_k_h = int_k_h'; %Transpose back the output