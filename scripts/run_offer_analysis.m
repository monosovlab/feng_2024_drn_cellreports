psnames = {ps_analysis_periods{1}.epochs.ps};
offi=1;

xname = intersect(glm_analysis_periods{1}(1).efits{10,1,1,3}.xname,glm_analysis_periods{2}(1).efits{314-224,1,1,3}.xname,'stable'); % get min set of attr across monkeys
xnamelabels = get_xnamelabels(xname);
halfi=3; 
bothoffi_relevant_bs = {};
bothoffi_relevant_ps = {};
bothoffi_relevant_ses = {};
[bothoffi_relevant_bs{1}, bothoffi_relevant_ps{1}, bothoffi_relevant_ses{1}] = get_relevant_bs_ps(glm_analysis_periods,1,halfi,monkeyabbr,xname,allu_ok);

savedir = neuralsavedir;
for mi=0:length(monkeyabbr)
    this_allu_ok = allu_ok;
    this_allu_ok_manex = allu_ok_manex;
    this_allu_ok_mintrials = allu_ok_mintrials;
    this_allu_ok_initin = allu_ok_initin;
    switch mi
        case 1
            this_allu_ok(225:end) = false;
            this_allu_ok_manex(225:end) = false; 
            this_allu_ok_mintrials(225:end) = false;
            this_allu_ok_initin(225:end) = false;
        case 2
            this_allu_ok(1:224) = false;
            this_allu_ok_manex(1:224) = false;
            this_allu_ok_mintrials(1:224) = false;
            this_allu_ok_initin(1:224) = false;
    end
    if mi
        monktxt = char(monkeyabbr(mi));
    else
        monktxt = 'both monks';
    end 

    relevant_bs = bothoffi_relevant_bs{offi};
    relevant_ps = bothoffi_relevant_ps{offi};
    relevant_ses = bothoffi_relevant_ses{offi};

    %%
    ERTRconds
    
    %% figure 3c
    if mi==0
        neural_proportion_bar_plots
    end
    
    %% figure 4; figure 5a; supp figure 2g; supp figure 4c
    RvsTout_scatter

    if mi==0
        %% figure 7b; supp figure 8c
        epochs
        
        %% figure 3b; figure 5b; figure 6a; supp figure 10c
        selective_subpop_psths 
        
        %% figure 7c; figure 7d; supp figure 8a; supp figure 8b; supp figure 8d; supp figure 8e; supp figure 9
        rpe_analysis
        
        %% figure 6b; figure 6c; supp figure 5; supp figure 6; supp figure 10a; supp figure 10b
        value_analyses
        
        %% supp figure 3a; supp figure 3b; supp table 2
        intrinsicproperties
        
        %% supp figure 3c, supp table 3
        DistanceFromAC
    end
end
