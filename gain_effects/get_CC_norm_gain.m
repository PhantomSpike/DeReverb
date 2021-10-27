function [CCnorm_small, CCnorm_big] = get_CC_norm_gain(coch_in, y_in, kernel, normz, animal_name)

%% Setup and params
coch = coch_in.coch;

if ismember(animal_name, ["Ronnie", "PLP", "Cork", "Kilkenny", "Derry"])
    small1_X_ft = coch(3).X_ft;
    small2_X_ft = coch(4).X_ft;
    
elseif ismember(animal_name, ["Noah", "Derekah"])
    small1_X_ft = coch(1).X_ft;
    small2_X_ft = coch(2).X_ft;
    
else
    error('Unrecognized animal');
    
end

big1_X_ft = coch(5).X_ft;
big2_X_ft = coch(6).X_ft;


n_h = kernel.n_h; %Find the number of history steps necessary given the bin size and max history 
t_edges_s = coch(1).t; %Get the edges in sec as used in the cochleagram

%% Make psths

% Make histograms for the 2 reverb conditions: stim
%For "Ronnie", "PLP", "Cork", "Kilkenny", "Derry" this is:
% [1,2]->anech, [3,4]->small reverb, [5,6]->large reverb
%For "Derekah" and "Noah" this is:
% [1,2]->small reverb, [3,4] ->medium reverb,[5,6]->large reverb

n_stim = length(y_in.stim);
n_rep = length(y_in.stim(1).repeat);

for s = 1:n_stim
    psth_temp = [];
    for r = 1:n_rep
        psth_temp(r,:) = histc(y_in.stim(s).repeat(r).spiketimes, t_edges_s);
    end
    y_temp{s} = psth_temp;
end

switch animal_name
    
    case {"Ronnie","PLP","Cork","Kilkenny","Derry"}
        y_t_small = [y_temp{3}, y_temp{4}];        
        
    case {"Noah","Derekah"}
        y_t_small = [y_temp{1}, y_temp{2}];
end

y_t_big = [y_temp{5}, y_temp{6}];

%% Optionally normalize the data
switch normz
    
    case {'perfreq','perfreq_noneuro'}
        %Calculate the mean and std for each frequency across all stimuli of the
        %same kind
        mean_reverb = mean([small1_X_ft, small2_X_ft, big1_X_ft, big2_X_ft],2); std_reverb = std([small1_X_ft, small2_X_ft, big1_X_ft, big2_X_ft],[],2);
        small1_X_ft = (small1_X_ft - mean_reverb)./std_reverb; small2_X_ft = (small2_X_ft - mean_reverb)./std_reverb;
        big1_X_ft = (big1_X_ft - mean_reverb)./std_reverb; big2_X_ft = (big2_X_ft - mean_reverb)./std_reverb;
        
    case 'none'
        
end

%% Tensorize
fprintf('== Tensorizing the cochleagrams ==\n');tic;
small1_X_fht = tensorize(small1_X_ft,n_h);
small2_X_fht = tensorize(small2_X_ft,n_h);
big1_X_fht = tensorize(big1_X_ft,n_h);
big2_X_fht = tensorize(big2_X_ft,n_h);

%% Concatenate together
small_X_fht = cat(3,small1_X_fht,small2_X_fht);
big_X_fht = cat(3,big1_X_fht,big2_X_fht);
n_t = size(small_X_fht,3);


%% Optionally Normalize neuronal firing too
switch normz
    case 'perfreq'
        y_t_small = (y_t_small - mean(y_t_small))./std(y_t_small);
        y_t_big = (y_t_big - mean(y_t_big))./std(y_t_big);
        
    case 'perfreq_noneuro'
        
    case 'none'
        
end


%% Compute CCnorm for the different folds and conditions
kernel_small = kernel.small.allKernels;
kernel_big = kernel.big.allKernels;

%Small
CCnorm_small = kfold_CCnorm(small_X_fht, y_t_small, kernel_small);

%Large
CCnorm_big = kfold_CCnorm(big_X_fht, y_t_big, kernel_big);


end