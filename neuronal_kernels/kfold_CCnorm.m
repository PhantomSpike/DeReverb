function res = kfold_CCnorm(X_fht, y_nt, kernel)

kfold = length(kernel) - 1; %Get the number of folds

%Params for validation
n_t = length(y_nt);
chunk_sz = floor(n_t/kfold);

%Run the model k times
for k = 1:kfold
    start_ix = (k-1)*chunk_sz + 1;
    end_ix = start_ix + chunk_sz -1;
    val_ix = [start_ix:end_ix]; %Get the indices corresponding to the validation set for the current fold
    
    curr_kernel = kernel{k}(1); %Get the kernel that was fit on the training set only
    curr_X_fht = X_fht(:,:,val_ix); %Get the cochleagram corresponding to the validation set for the current fold
    R = y_nt(:,val_ix); %Get the psths corresponding to the validation set for the current fold
    
    %Predict the firing rate in the current fold
    y_hat_t = kernelconv(curr_X_fht, curr_kernel);
    
    %Calculate CCnorm
    [CCnorm, CCabs, CCmax] = calc_CCnorm(R, y_hat_t);
    
    res.CCnorm(k) = CCnorm;
    res.CCabs(k) = CCabs;
    res.CCmax(k) = CCmax;
end

res.CCnorm_mean = nanmean(res.CCnorm);
res.CCabs_mean = nanmean(res.CCabs);
res.CCmax_mean = nanmean(res.CCmax);
res.hasNaN = any(isnan(res.CCnorm));