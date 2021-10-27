kh_fits_name{1} = '/mnt/40086D4C086D41D0/Reverb_analysis/dexp_fits/open/40pos/2.5ms/2ms_open manypos spechill sep none kernel_ub_A_B_25/kh_fits.mat';
% kh_fits_name{1} = '/mnt/40086D4C086D41D0/Reverb_analysis/dexp_fits/open/40pos/2.5ms/2ms_open manypos spechill sep none kernel_ub_A_B_25/kh_fits.mat';
% kh_fits_name{2} = '/mnt/40086D4C086D41D0/Reverb_analysis/dexp_fits/open/40pos/5ms/open 40pos 5ms 150_spechill sep none kernel_ub_A_B_25/kh_fits.mat';
% kh_fits_name{3} = '/mnt/40086D4C086D41D0/Reverb_analysis/dexp_fits/open/40pos/5ms/open 40pos 5ms 300_spechill sep none kernel_ub_A_B_25/kh_fits.mat';
% kh_fits_name{4} = '/mnt/40086D4C086D41D0/Reverb_analysis/dexp_fits/open/40pos/10ms/open 40pos 10ms_spechill sep none kernel_ub_A_B_25/kh_fits.mat';
% kh_fits_name{5} = '/mnt/40086D4C086D41D0/Reverb_analysis/dexp_fits/open/single_pos/5ms/open onepos 5ms 150_spechill sep none kernel_ub_A_B_25/kh_fits.mat';
% kh_fits_name{6} = '/mnt/40086D4C086D41D0/Reverb_analysis/dexp_fits/open/single_pos/5ms/open onepos 5ms 300_spechill sep none kernel_ub_A_B_25/kh_fits.mat';
% kh_fits_name{7} = '/mnt/40086D4C086D41D0/Reverb_analysis/dexp_fits/open/single_pos/10ms/open onepos 10ms 300_spechill sep none kernel_ub_A_B_25/kh_fits.mat';
% kh_fits_name{8} = '/mnt/40086D4C086D41D0/Reverb_analysis/dexp_fits/closed/chopped/40pos/2ms/2ms_closed manypos spechill sep none kernel_ub_A_B_25/kh_fits.mat';
% kh_fits_name{9} = '/mnt/40086D4C086D41D0/Reverb_analysis/dexp_fits/closed/chopped/40pos/5ms/closed 40pos 5ms 150_spechill sep none kernel_ub_A_B_25/kh_fits.mat';
% kh_fits_name{10} = '/mnt/40086D4C086D41D0/Reverb_analysis/dexp_fits/closed/chopped/40pos/5ms/closed 40pos 5ms 300_spechill sep none kernel_ub_A_B_25/kh_fits.mat';
% kh_fits_name{11} = '/mnt/40086D4C086D41D0/Reverb_analysis/dexp_fits/closed/chopped/40pos/10ms/closed 40pos 10ms 300_spechill sep none kernel_ub_A_B_25/kh_fits.mat';
% kh_fits_name{12} = '/mnt/40086D4C086D41D0/Reverb_analysis/dexp_fits/closed/chopped/single_pos/5ms/onepos chopped 5ms_spechill sep none kernel_ub_A_B_25/kh_fits.mat';

load('/mnt/40086D4C086D41D0/Reverb_analysis/Model_fits/open/single_pos/10ms/open_onepos_10ms_300/spechill_sep_none_kernel.mat');

r_type{1} = 'small';
r_type{2} = 'med';
r_type{3} = 'big';
n_rooms = length(r_type);
ub_A_B = 25;
int_name = [];
n_ker = 30; %The number of kernels

params.r_type = r_type;
params.int_name = int_name;
params.freqs = kernels.freqs;
params.n_rooms = n_rooms;
params.n_ker = n_ker;
params.ub_A_B = ub_A_B;

for k = 1:length(kh_fits_name)
    load(kh_fits_name{k});
    [save_dir_full,~] = fileparts(kh_fits_name{k});
    params.save_dir_full = save_dir_full;
    plot_tconst(kh_fits,params);
end



