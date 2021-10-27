kernel_name{1} = '/mnt/40086D4C086D41D0/Reverb_analysis/Model_fits/open/many/2ms/open_manypos_spechill_sep_none_kernel.mat';
coch_name{1} = '/mnt/40086D4C086D41D0/Reverb_analysis/Cochleagrams/open/many/2.5ms/coch_all_conditions_open_manypos_spechill.mat';

% kernel_name{2} = '/mnt/40086D4C086D41D0/Reverb_analysis/Model_fits/closed_40pos_5ms_300/spechill_sep_none_kernel.mat';
% coch_name{2} = '/mnt/40086D4C086D41D0/Reverb_analysis/Cochleagrams/closed_manypos/5ms_tbin/spechill/coch_all_conditions_spechill.mat';
% 
% kernel_name{3} = '/mnt/40086D4C086D41D0/Reverb_analysis/Model_fits/closed_40pos_10ms_300/spechill_sep_none_kernel.mat';
% coch_name{3} = '/mnt/40086D4C086D41D0/Reverb_analysis/Cochleagrams/closed_manypos/10ms_tbin/spechill/coch_all_conditions_spechill.mat';



for j = 1:length(kernel_name)
    fit_tconst_model(kernel_name{j},coch_name{j});
end