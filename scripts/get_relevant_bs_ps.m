function [relevant_bs, relevant_ps, relevant_ses] = get_relevant_bs_ps(glm_analysis_periods,offi,halfi,monkeyabbr,xname,allu_ok)
    relevant_bs = NaN(length(allu_ok),length(xname));
    relevant_ps = relevant_bs;
    relevant_ses = relevant_bs;

    relevant_bs_ = [];
    relevant_ps_ = [];
    relevant_ses_ = [];

    for gi=1:length(glm_analysis_periods)
        glm_efit = [glm_analysis_periods{gi}(offi).efits{:,1,1,halfi}];
        bs = [glm_efit.b]';
        ss = [glm_efit.stats];
        ps = [ss.p]';
        ses = [ss.se]';
        xnamei = find(ismember(glm_efit(1).xname,xname)); % for off2, xname will be unique

        tmpbs = [];
        tmpps = [];
        tmpses = [];
        for tmpoffi=1:offi
            tmpi=(length(xname)-1)*(tmpoffi-1)+(1:(length(xname)-1));
            if strcmp(monkeyabbr(gi),'P') && offi==2 && tmpoffi==1 % Monkey P matched info - so a caveat is while this is indeed the effect of offer2 being info, it could also be explained by off1
                tmpi=[tmpi(1:2) find(strcmp(glm_efit(1).xname,'I')) tmpi(3:end-1)]; % remove last one because it will be first reg of off2
            end
            tmpi = [tmpi length(xnamei)]; % add bias back in, will be repeated for off2
            tmpbs = cat(3,tmpbs,bs(:,xnamei(tmpi)));
            tmpps = cat(3,tmpps,ps(:,xnamei(tmpi)));
            tmpses = cat(3,tmpses,ses(:,xnamei(tmpi)));
        end

        relevant_bs_ = [relevant_bs_; tmpbs];
        relevant_ps_ = [relevant_ps_; tmpps];
        relevant_ses_ = [relevant_ses_; tmpses];
    end

    relevant_bs(allu_ok,:) = relevant_bs_;
    relevant_ps(allu_ok,:) = relevant_ps_;
    relevant_ses(allu_ok,:) = relevant_ses_;
end