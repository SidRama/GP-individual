function compute_KL(probeset_file, file_save)

    %   Code to perform GP fittings and compute KL scores.
    %   Input:
    %       probeset_file - .mat file containing probe-set measurements amd 
    %                       associated covariates. This file also contains
    %                       the disease age (sero_age)
    %       file_save - results file name (ending with .mat)
    %   Output:
    %        .mat file with KL-scores for each probe-set
    % 


    load(probeset_file);
    age_case = age_case*4;
    age_control = age_control*4;
    age_combined = [age_case; age_control];  
    [age_combined, sortIndex] = sort(age_combined);
    total_probes = size(case_data, 1)
    ratios = cell(total_probes, 1);
    ratios_bf = cell(total_probes, 1);

    window_start = sero_age - 26;
    window_end = sero_age;
    
    min_age = min([window_start, age_combined(1)]); 
    x_test = (min_age:age_combined(end))';

    iter_track = 1;
    for gene_num = 1:total_probes

        % obtain case, control and combined covariates.
        sample_case = case_data(gene_num, :)';
        sample_control = control_data(gene_num, :)';
        
        age_case = age_case_original;
        age_control = age_control_original;
        
        sample_combined = [sample_case; sample_control];
        age_combined = [age_case; age_control];  
        [age_combined, sortIndex] = sort(age_combined);
        x_test = (age_combined(1):age_combined(end))';
        % sort with respect to age
        sample_combined = sample_combined(sortIndex);

        % mean centring
        yMean_case = mean(sample_case);
        sample_centred_case = sample_case - yMean_case;

        yMean_control = mean(sample_control);
        sample_centred_control = sample_control - yMean_control;

        yMean_combined = mean(sample_combined);
        sample_centred_combined = sample_combined - yMean_combined;

        age_list = unique(age_combined);

        % case GP
        lik = lik_gaussian('sigma2',sinvchi2rand(1, 0.01));
        pn = prior_sinvchi2('s2', 0.01, 'nu', 1);
        lik = lik_gaussian(lik,'sigma2_prior', pn);

        pl = prior_gaussian('mu', 40, 's2', 4);
        pm = prior_sqrtt('s2',1,'nu',20);
        gpcf = gpcf_sexp('lengthScale',normrnd(40,2), 'magnSigma2', 0.2^2,...
            'lengthScale_prior', pl, 'magnSigma2_prior', pm);
        gp = gp_set('lik', lik, 'cf', gpcf);

        [gp_case_arr, p_th_case, th_case] = gp_ia(gp, age_case, sample_centred_case);

        [Ef_case, Varf_case] = gpia_pred(gp_case_arr, age_case, sample_centred_case, x_test);
        Eff_case = Ef_case + yMean_case;
        std_ft_case = sqrt(Varf_case);

        % control GP
        lik = lik_gaussian('sigma2',sinvchi2rand(1, 0.01));
        pn = prior_sinvchi2('s2', 0.01, 'nu', 1);
        lik = lik_gaussian(lik,'sigma2_prior', pn);

        pl = prior_gaussian('mu', 40, 's2', 4);
        pm = prior_sqrtt('s2',1,'nu',20);
        gpcf = gpcf_sexp('lengthScale',normrnd(40,2), 'magnSigma2', 0.2^2,...
            'lengthScale_prior', pl, 'magnSigma2_prior', pm);
        gp = gp_set('lik', lik, 'cf', gpcf);

        [gp_control_arr, p_th_control, th_control] = gp_ia(gp, age_control, sample_centred_control);
        [Ef_control, Varf_control] = gpia_pred(gp_control_arr, age_control, sample_centred_control, x_test);
        Eff_control = Ef_control + yMean_control;
        std_ft_control = sqrt(Varf_control);

        % combined GP
        lik = lik_gaussian('sigma2',sinvchi2rand(1, 0.01));
        pn = prior_sinvchi2('s2', 0.01, 'nu', 1);
        lik = lik_gaussian(lik,'sigma2_prior', pn);

        pl = prior_gaussian('mu', 40, 's2', 4);
        pm = prior_sqrtt('s2',1,'nu',20);
        gpcf = gpcf_sexp('lengthScale',normrnd(40,2), 'magnSigma2', 0.2^2,...
            'lengthScale_prior', pl, 'magnSigma2_prior', pm);
        gp = gp_set('lik', lik, 'cf', gpcf);

        [gp_combined_arr, p_th_combined, th_combined] = gp_ia(gp, age_combined, sample_centred_combined);
        [Ef_combined, Varf_combined] = gpia_pred(gp_combined_arr, age_combined, sample_centred_combined, x_test);
        Eff_combined = Ef_combined + yMean_combined;
        std_ft_combined = sqrt(Varf_combined);

        % Bayes Factor computation - marginal likelihood

        % separate model
        % case 
        ml_case = 0;
        for cur = 1:length(gp_case_arr)
            w = gp_pak(gp_case_arr{cur});
            [e, ~,~ ] = gp_e(w, gp_case_arr{cur}, age_case, sample_centred_case);
            ml_case = ml_case + (exp(-e) * gp_case_arr{cur}.ia_weight) ;
        end

        % control
        ml_control = 0;
        for cur = 1:length(gp_control_arr)
            w = gp_pak(gp_control_arr{cur});
            [e, ~, ~] = gp_e(w, gp_control_arr{cur}, age_control, sample_centred_control);
            ml_control = ml_control + (exp(-e) * gp_control_arr{cur}.ia_weight);
        end

        % joint model
        % combined
        ml_combined = 0;
        for cur = 1:length(gp_combined_arr)
            w = gp_pak(gp_combined_arr{cur});
            [e, ~, ~] = gp_e(w, gp_combined_arr{cur}, age_combined, sample_centred_combined);
            ml_combined = ml_combined + (exp(-e) * gp_combined_arr{cur}.ia_weight);
        end

        % BF-score
        ml_ratio = log(ml_case) + log(ml_control) - log(ml_combined);
        ratios_bf{iter_track, 1} = ml_ratio;
        
    
        % Symmetric Kullback-Leibler computation for KL-scores
        index_start = find(x_test == window_start);
        index_end = find(x_test == window_end);
        % window computation
        mu_win_1 = [Eff_case(index_start:index_end); Eff_control(index_start:index_end)];
        vars = [Varf_case(index_start:index_end); Varf_control(index_start:index_end)];
        cov_win_1 = diag(vars);

        mu_win_2 = [Eff_combined(index_start:index_end); Eff_combined(index_start:index_end)];
        vars = [Varf_combined(index_start:index_end); Varf_combined(index_start:index_end)];
        cov_win_2 = diag(vars);

        kl_win_1 = 0.5 * (log(det(cov_win_2)/det(cov_win_1)) - length(mu_win_1) + trace(pinv(cov_win_2) * cov_win_1) ...
                + ((mu_win_2 - mu_win_1)' * pinv(cov_win_2) * (mu_win_2 - mu_win_1)));


        mu_win_1 = [Eff_combined(index_start:index_end); Eff_combined(index_start:index_end)];
        vars = [Varf_combined(index_start:index_end); Varf_combined(index_start:index_end)];
        cov_win_1 = diag(vars);

        mu_win_2 = [Eff_case(index_start:index_end); Eff_control(index_start:index_end)];
        vars = [Varf_case(index_start:index_end); Varf_control(index_start:index_end)];
        cov_win_2 = diag(vars);

        kl_win_2 = 0.5 * (log(det(cov_win_2)/det(cov_win_1)) - length(mu_win_1) + trace(pinv(cov_win_2) * cov_win_1) ...
                + ((mu_win_2 - mu_win_1)' * pinv(cov_win_2) * (mu_win_2 - mu_win_1)));


        kl_win_score = 0.5 * (kl_win_1 + kl_win_2); 
        
        ratios{iter_track, 1} = kl_win_score;
      
        if (mod(gene_num, 200) == 0)
            fprintf('\n\n\n %i probe-sets complete\n\n\n', gene_num)
        end
        iter_track = iter_track + 1;
    end
    
    file_save = strcat(file_save, '.mat');
    save(file_save)
end