function ker = run_model_kfold(reverb_X_fht,y_t,model)

kfold = 10;
fprintf('== Fitting %s model ==\n',model);tic

%Params for validation
n_t = length(y_t);
chunk_sz = floor(n_t/kfold);

%Run the model k times
for k = 1:kfold
    fprintf('== Fold %0.f/%0.f ==\n',k,kfold);
    start_ix = (k-1)*chunk_sz + 1;
    end_ix = start_ix + chunk_sz -1;
    val_ix = [start_ix:end_ix];
    [ker.main{k}, ker.allKernels{k}, ker.res{k}] = alexnet_fht(reverb_X_fht, y_t, model, val_ix);
end
fprintf('== Done! ==\n');


fprintf('== Fitting best model ==\n');
%Find best lambda
for k = 1:kfold
    lambdas(k) = ker.main{k}.lambda;
end

%Because the spacing of lambdas is log find the log weighted mean
best_lambda = 10.^(mean(log10(lambdas)));

%Run the model with the best lambda
options = glmnetSet;

switch model
    case 'lasso'
        alpha = 1;
    case 'ridge'
        alpha = 0.01;
end

options.alpha = alpha;
options.lambda = best_lambda;
[ker.main{kfold+1}, ker.allKernels{kfold+1}, ker.res{kfold+1}] = alexnet_fht(reverb_X_fht, y_t, options);
fprintf('== Done! ==\n');