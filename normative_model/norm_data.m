function [reverb_norm_X_ft,anech_norm_X_ft] = norm_data(reverb_X_ft,anech_X_ft,type)
%[reverb_norm_X_ft,anech_norm_X_ft] = norm_data(reverb_X_ft,anech_X_ft,type)
%This functon normalizes (standardise) the data before model fitting. The normalization
%is either for each freqeuncy separately or across all the data

reverb_norm_X_ft = zeros(size(reverb_X_ft));
anech_norm_X_ft = zeros(size(anech_X_ft));


switch type
    case 'perfreq'
        reverb_norm_X_ft = (reverb_X_ft - mean(reverb_X_ft,2))./std(reverb_X_ft,[],2);
        anech_norm_X_ft = (anech_X_ft - mean(anech_X_ft,2))./std(anech_X_ft,[],2);
    case 'global'
        reverb_norm_X_ft = (reverb_X_ft - mean(reverb_X_ft(:)))/std(reverb_X_ft(:));
        anech_norm_X_ft = (anech_X_ft - mean(anech_X_ft,2))./std(anech_X_ft,[],2);
end

