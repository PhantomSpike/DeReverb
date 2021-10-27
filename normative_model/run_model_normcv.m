function ker = run_model_normcv(reverb_X_fht,anechoic_X_ft,model,room_sz)

kfold = 10; %How many folds for the model

if isnumeric(model)
    model_name = 'manual_alpha';
else
    model_name = model;
end
fprintf('== Fitting %s model for %s room ==\n',model_name,room_sz);tic
n_f = size(anechoic_X_ft,1);
ker = cell(n_f,1);
parfor k = 1:n_f
    fprintf('== Fitting frequency band %0.f/%0.f ==\n',k,n_f);
    if ischar(model) || isnumeric(model)
        [ker{k}.main, ker{k}.allKernels, ker{k}.res] = kfold_model(reverb_X_fht, anechoic_X_ft(k,:), model, kfold);
    else
        error('Unrecognised model type')
    end
end
fprintf('== Done! This took %0.fs ==\n',toc);
end