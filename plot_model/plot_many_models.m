model_dir = '/mnt/40086D4C086D41D0/Reverb_analysis/Model_fits/open_40pos_5ms_150';
coch_dir = '/mnt/40086D4C086D41D0/Reverb_analysis/Cochleagrams/open_manypos/5ms_tbin/norm_pos';
%Try different Cochleagrams
coch_type{1} = 'spechill';
% coch_type{2} = 'specpower';
% coch_type{3} = 'speclog';
num_coch = length(coch_type);

%Try different models
model{1} = 'sep';
% model{2} = 'elnet_lasso';
% model{3} = 'elnet_ridge';
num_model = length(model);

%Try different normalizations
normz{1} = 'none';
% normz{2} = 'perfreq';
% normz{2} = 'global';

num_normz = length(normz);

j=0;
for c = 1:num_coch
    for m = 1:num_model
        for n = 1:num_normz
            try
                plot_norm_model(coch_type{c},coch_dir,model{m},normz{n},model_dir);
            catch 
                j=j+1;
                er{j} = ['Error with ',strjoin({coch_type{c},model{m},normz{n}},'_')];
            end
        end
    end
end
