function ker = run_model_norm_test(reverb_X_fht,anechoic_X_ft,model_type,room_sz)

fprintf('== Fitting %s model for %s room ==\n',model_type,room_sz);tic
n_f = 1;
for k = 1:n_f
    fprintf('== Fitting frequency band %0.f/%0.f ==\n',k,n_f);
    switch model_type
        case 'sep'
            ker{k} = sepkerneltensor2(reverb_X_fht, anechoic_X_ft(k,:));
        case 'elnet_lasso'
            [ker{k}.main, ker{k}.allKernels, ker{k}.res] = elnet_fht(reverb_X_fht, anechoic_X_ft(k,:), 'lasso');
        case 'elnet_ridge'
            [ker{k}.main, ker{k}.allKernels, ker{k}.res] = elnet_fht(reverb_X_fht, anechoic_X_ft(k,:), 'ridge');
    end
end
fprintf('== Done! This took %0.fs ==\n',toc);
end