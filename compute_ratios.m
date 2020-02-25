function compute_ratios(probeset_file, file_save)

    %   Code to perform GP fittings and compute BF scores.
    %   Input:
    %       probeset_file - .mat file containing probe-set (or gene) 
    %                       measurements and associated covariates.
    %       file_save - results file name (ending with .mat)
    %   Output:
    %        .mat file with BF-scores for each feature
    % 
    

    load(probeset_file);    
    age_case_original = age_case;
    age_control_original = age_control;
    total_probes = size(case_data, 1)

    ratios = cell(total_probes, 1);
    
    conv_fail= [];
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
        ratios{iter_track, 1} = ml_ratio;
        
	
        if (mod(gene_num, 200) == 0)
            fprintf('\n\n\n %i probe-sets complete\n\n\n', gene_num)
        end
        iter_track = iter_track + 1;
    end
    
    file_save = strcat(file_save, '.mat');
    save(file_save)
end
