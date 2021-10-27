function ker = run_model(reverb_X_fht,y_t,model_type,k_f)
fprintf('== Fitting %s model ==\n',model_type);tic

%Params for validation
n_t = length(y_t); 
val_per = 0.1; %What percentage of the data to use for validation
val_sz = floor(n_t*val_per); %Convert to samples
start_ix = randi([1 n_t-val_sz],1,1); %Define random starting index
val_ix = [start_ix:start_ix+val_sz]; %Define the val indices

switch model_type
    case 'sep'
        ker = sepkerneltensor2(reverb_X_fht, y_t);
    case {'sep_kh_all','sep_kh'}
        if ~exist('k_f','var') || isempty(k_f)
            error('k_f not supplied');
        end
        ker =  get_k_h(reverb_X_fht, y_t, k_f);
    case {'lasso','ridge','elastic'}
        [ker.main, ker.allKernels, ker.res] = alexnet_fht(reverb_X_fht, y_t, model_type,val_ix);
    case 'dnet'
        lambda = 1e-7;
        n_hidden = 20;
        n_pass = 40;
        n_minibatch = 20;
        sz_minibatch = floor(length(y_t)/n_minibatch);
        for b = 1:n_minibatch
            ix_start = (b-1)*sz_minibatch + 1;
            ix_end = ix_start + sz_minibatch;
            train_nfht(b,:,:,:) =  reverb_X_fht(:,:,ix_start:ix_end);
            train_nt(b,:) = y_t(ix_start:ix_end);
        end
        [ker.theta,ker.train_err]=fit_DNet_model(train_nfht,train_nt,'sq','abs',lambda,{n_hidden 1},'sigmoid',n_pass,'mDNet');
end
fprintf('== Done! This took %0.fs ==\n',toc);
end