% rpe_analysis.m
scriptname = 'rpe_analysis';

%% calc stats 
psapnames = arrayfun(@(x)x.ps,ps_analysis_periods{1}.epochs,'UniformOutput',false);
epochparams = struct('ps',{'cue','rev'},'info',{1,0},'anawind',{[250 750],[250 1050]});

cicolors = {[1 0.5 0],[1 0.75 0.5];[0.5 0 1],[0.75 0.5 1]};
lcolors = [0 0 1; 0 0 0; 1 0 0];
statmode = {'diff(-,+)','diff(0,+)','diff(-,0)','roc(-,+)','roc(0,+)','roc(-,0)','mean(roc(0,+),roc(-,0))'};

for do255025only = 0:1
    if do255025only
        trialconds = struct(...
            'name',{'Info RPE +','Info RPE 0','Info RPE -';'Noinfo RPE +','Noinfo RPE 0','Noinfo RPE -'}, ...
            'attributes',{{'rpe','info','entropy'},{'rpe','info','entropy'},{'rpe','info','entropy'};{'rpe','info','entropy'},{'rpe','info','entropy'},{'rpe','info','entropy'}}, ...
            'levels',{[3 2 3],[2 2 3],[1 2 3];[3 1 3],[2 1 3],[1 1 3]}, ...
            'color',{cicolors{1,1},cicolors{1,2},cicolors{1,1};cicolors{2,1},cicolors{2,2},cicolors{2,1}},...
            'style',{cistyles{1},cistyles{1},cistyles{2};cistyles{1},cistyles{1},cistyles{2}});
    else
        trialconds = struct(...
            'name',{'Info RPE +','Info RPE 0','Info RPE -';'Noinfo RPE +','Noinfo RPE 0','Noinfo RPE -'}, ...
            'attributes',{{'rpe','info'},{'rpe','info'},{'rpe','info'};{'rpe','info'},{'rpe','info'},{'rpe','info'}}, ...
            'levels',{[3 2],[2 2],[1 2];[3 1],[2 1],[1 1]}, ...
            'color',{cicolors{1,1},cicolors{1,2},cicolors{1,1};cicolors{2,1},cicolors{2,2},cicolors{2,1}},...
            'style',{cistyles{1},cistyles{1},cistyles{2};cistyles{1},cistyles{1},cistyles{2}});
    end
    
    trialconds = flip(trialconds,2);
    
    anastat = NaN(length(allu_ok),length(epochparams),length(statmode));
    anastatp = anastat;
    
    %%
    for epi=1:size(anastat,2)
        %%
        api = find(strcmp(psnames,epochparams(epi).ps));
        tmp = max(size(ps_analysis_periods{1}.epochs(api).ys,1),size(ps_analysis_periods{2}.epochs(api).ys,1));
        level = {};
        ys = {};
        for mi_=1:2
            level{mi_} = ps_analysis_periods{mi_}.epochs(api).level;
            ys{mi_} = ps_analysis_periods{mi_}.epochs(api).ys;
            if size(ps_analysis_periods{mi_}.epochs(api).ys,1)<tmp
                tmpdiff = tmp-size(ps_analysis_periods{mi_}.epochs(api).level,1);
                
                tmp_level = size(ps_analysis_periods{mi_}.epochs(api).level);
                tmp_level(1) = tmpdiff;
                level{mi_} = vertcat(level{mi_},NaN(tmp_level));
                tmp_ys = size(ps_analysis_periods{mi_}.epochs(api).ys);
                tmp_ys(1) = tmpdiff;
                ys{mi_} = vertcat(ys{mi_},NaN(tmp_ys));
            end
        end
        t = ps_analysis_periods{1}.epochs(api).t;
        attributes = ps_analysis_periods{1}.epochs(api).attributes;
        level = cat(2,level{1},level{2});
        ys = cat(3,ys{1},ys{2});
        noktrall = cat(2,sum(ps_analysis_periods{1}.oktrall),sum(ps_analysis_periods{2}.oktrall));

        for ui=1:length(allu_ok)
            if epochparams(epi).info
                these_trialconds = trialconds(1,:);
            else
                these_trialconds = trialconds(2,:);
            end
            rates = cell(size(these_trialconds));
            for tci=1:length(these_trialconds)
                tmpys = squeeze(ys(1:noktrall(ui),:,ui));
                tmplevel = squeeze(level(1:noktrall(ui),ui,ismember(attributes,these_trialconds(tci).attributes),3));
                tmp = tmpys(all(tmplevel==these_trialconds(tci).levels,2),:);
                rates{tci} = nanmean(tmp(:,inbounds(t,epochparams(epi).anawind)),2);
            end
            if any(cellfun(@(x)length(x),rates)<2)
                continue;
            end
            if 1
                for stmi=1:size(anastat,3)
                    switch statmode{stmi}
                        case 'diff(-,+)'
                            anastat(ui,epi,stmi) = mean(rates{3})-mean(rates{1});
                            anastatp(ui,epi,stmi) = ranksum(rates{1},rates{3});
                        case 'diff(0,+)'
                            anastat(ui,epi,stmi) = mean(rates{3})-mean(rates{2});
                            anastatp(ui,epi,stmi) = ranksum(rates{2},rates{3});
                        case 'diff(-,0)'
                            anastat(ui,epi,stmi) = mean(rates{2})-mean(rates{1});
                            anastatp(ui,epi,stmi) = ranksum(rates{1},rates{2});
                        case 'roc(-,+)'
                            [anastat(ui,epi,stmi),anastatp(ui,epi,stmi)] = rocarea3(rates{1},rates{3});
                        case 'roc(0,+)'
                            [anastat(ui,epi,stmi),anastatp(ui,epi,stmi)] = rocarea3(rates{2},rates{3});
                        case 'roc(-,0)'
                            [anastat(ui,epi,stmi),anastatp(ui,epi,stmi)] = rocarea3(rates{1},rates{2});
                        case 'mean(roc(0,+),roc(-,0))'
                            [neg,negp] = rocarea3(rates{1},rates{2});
                            [pos,posp] = rocarea3(rates{2},rates{3});
                            anastat(ui,epi,stmi) = mean([neg pos]);
                            anastatp(ui,epi,stmi) = combine_pvalues([negp posp]);
                        case 'diff+(2,3)'
                            rpelevel = 3;
                            distlevel = [2 3];
                            signallevel = rpelevel+3*(distlevel(2)-1);
                            noiselevel = rpelevel+3*(distlevel(1)-1);
                            if ~isempty(rates{signallevel}) && ~isempty(rates{noiselevel})
                                anastat(ui,epi,stmi) = mean(rates{signallevel})-mean(rates{noiselevel});
                                anastatp(ui,epi,stmi) = ranksum(rates{signallevel},rates{noiselevel});
                            end
                        case 'diff-(2,3)'
                            rpelevel = 1;
                            distlevel = [2 3];
                            signallevel = rpelevel+3*(distlevel(2)-1);
                            noiselevel = rpelevel+3*(distlevel(1)-1);
                            if ~isempty(rates{signallevel}) && ~isempty(rates{noiselevel})
                                anastat(ui,epi,stmi) = mean(rates{signallevel})-mean(rates{noiselevel});
                                anastatp(ui,epi,stmi) = ranksum(rates{signallevel},rates{noiselevel});
                            end
                        case 'diff0(1,2)'
                            rpelevel = 2;
                            distlevel = [1 2];
                            signallevel = rpelevel+3*(distlevel(2)-1);
                            noiselevel = rpelevel+3*(distlevel(1)-1);
                            if ~isempty(rates{signallevel}) && ~isempty(rates{noiselevel})
                                anastat(ui,epi,stmi) = mean(rates{signallevel})-mean(rates{noiselevel});
                                anastatp(ui,epi,stmi) = ranksum(rates{signallevel},rates{noiselevel});
                            end
                        case 'roc+(2,3)'
                            rpelevel = 3;
                            distlevel = [2 3];
                            signallevel = rpelevel+3*(distlevel(2)-1);
                            noiselevel = rpelevel+3*(distlevel(1)-1);
                            [anastat(ui,epi,stmi),anastatp(ui,epi,stmi)] = rocarea3(rates{noiselevel},rates{signallevel});
                        case 'roc-(2,3)'
                            rpelevel = 1;
                            distlevel = [2 3];
                            signallevel = rpelevel+3*(distlevel(2)-1);
                            noiselevel = rpelevel+3*(distlevel(1)-1);
                            [anastat(ui,epi,stmi),anastatp(ui,epi,stmi)] = rocarea3(rates{noiselevel},rates{signallevel});
                        case 'roc0(1,2)'
                            rpelevel = 2;
                            distlevel = [1 2];
                            signallevel = rpelevel+3*(distlevel(2)-1);
                            noiselevel = rpelevel+3*(distlevel(1)-1);
                            [anastat(ui,epi,stmi),anastatp(ui,epi,stmi)] = rocarea3(rates{noiselevel},rates{signallevel});
                        case 'mean(roc+(2,3),roc-(2,3),roc0(1,2))'
                            rpelevel = 3;
                            distlevel = [2 3];
                            signallevel = rpelevel+3*(distlevel(2)-1);
                            noiselevel = rpelevel+3*(distlevel(1)-1);
                            [pos,posp] = rocarea3(rates{noiselevel},rates{signallevel});
    
                            rpelevel = 1;
                            distlevel = [2 3];
                            signallevel = rpelevel+3*(distlevel(2)-1);
                            noiselevel = rpelevel+3*(distlevel(1)-1);
                            [neg,negp] = rocarea3(rates{noiselevel},rates{signallevel});
    
                            rpelevel = 2;
                            distlevel = [1 2];
                            signallevel = rpelevel+3*(distlevel(2)-1);
                            noiselevel = rpelevel+3*(distlevel(1)-1);
                            [zer,zerp] = rocarea3(rates{noiselevel},rates{signallevel});
                            
                            anastat(ui,epi,stmi) = mean([neg zer pos]);
                            anastatp(ui,epi,stmi) = combine_pvalues([negp zerp posp]);
                    end
                end
            end
        end
    end
    if do255025only==0
        anastatfull = anastat;
        anastatpfull = anastatp;
    end

    %% plots
    do_plot = 1;
    conds = ertr_conds;
    bargraphconds = [3 5 2 4];
    histconds = [7 6];
    meanconds = flip(histconds);
    
    stmis = [1 4 7];
    epises = {1,2,1:2};
    for stmi = 7
        x0 = ternary(~isempty(regexp(statmode{stmi},'roc', 'once')),0.5,0); y0 = x0;
        for ei=3
            switch ei
                case {1,2}
                    etag = [ternary(epochparams(epises{ei}).info,'Info','Noinfo') ' ' epochparams(epises{ei}).ps];
                case 3
                    etag = 'combined';
            end
            indextag = [statmode{stmi} ', ' etag];
            
            %% for different conds, compare rpe dists
            thisanastat = mean(anastat(:,epises{ei},stmi),2);
            thisanastatp = combine_pvalues(anastatp(:,epises{ei},stmi),2);
            
            x_s = {};
            sig_s = {};
            for condi=1:length(conds)
                x_ = thisanastat(conds(condi).thisalluok);
                x_s{condi} = x_;
                sig_ = thisanastatp(conds(condi).thisalluok)<sigthresh;
                sig_s{condi} = sig_;
            end
            rsps = NaN(length(conds));
            avgdiffs = rsps;
            srps = NaN(1,length(conds));
            avgs = srps;
            ses = srps;
            nsig = srps;
            for condi_=1:length(conds)
                if isempty(x_s{condi_})
                    continue
                end
                avgs(condi_) = mean(x_s{condi_});
                ses(condi_) = std(x_s{condi_})/sqrt(length(x_s{condi_}));
                srps(condi_)=signrank(x_s{condi_}-x0);
                nsig(condi_) = sum(sig_s{condi_}); 
                for condj_=1:length(conds)
                    if condi_>condj_
                        if isempty(x_s{condj_})
                            continue
                        end
                        avgdiffs(condi_,condj_)=avgs(condj_)-avgs(condi_);
                        rsps(condi_,condj_)=ranksum(x_s{condi_},x_s{condj_});
                    end
                end
            end
            T = table({conds.name}',condn',nsig',avgs',ses',srps');
            
            if do_plot 
                %% rpe vs offer scatter
                attri = 1;
                basettl = ['rpe_vs_' xn{attri,2}]; 
                basettl = [monktxt ', Offer ' num2str(offi) ' conds, ' basettl];
                
                x = relevant_bs(this_allu_ok,strcmp(xname,xn{attri,1})); xp = relevant_ps(this_allu_ok,strcmp(xname,xn{attri,1}));
                 y = thisanastat(this_allu_ok); yp = thisanastatp(this_allu_ok);
                y0 = ternary(~isempty(regexp(statmode{stmi},'roc', 'once')),0.5,0);
                maxlim = {[-0.45 0.45],max(abs(y))*[-1 1]};
                if y0==0.5
                    maxlim{2} = [0 1];
                end
                
                % with hist
                figure('Visible','on');
                nans_to_remove = any(isnan([x y]),2); % remove NaN entries here where y could not be computed
                myScatterWithHistogram([xn{attri,2} ' GLM weight'],escape(indextag),x(~nans_to_remove),y(~nans_to_remove),xp(~nans_to_remove),yp(~nans_to_remove),[],[],[],maxlim,0,y0,[],sigthresh);
    
                if do255025only
                    figname = 'supp_figure_8de';
                else
                    figname = 'figure_7cd';
                end

                savepath=fullfile(neuralsavedir,figname);
                saveas(gcf,[savepath '.fig']);
                saveas(gcf,[savepath '.png']);
                close(gcf);
                
                if ~do255025only
                    %% mostly bidirectional (positively signed) rpes, justify combining + and - component
                    basettl = 'rpe_+vs-';
                    basettl = [monktxt ', Offer ' num2str(offi) ' conds, ' basettl];
        
                    if ~isempty(regexp(statmode{stmi},'mean', 'once'))
                        posstmi=4+1;
                    else
                        posstmi=stmi+1;
                    end
                    xylabels = escape({[statmode{posstmi} ', ' etag],[statmode{posstmi+1} ', ' etag]});
                    x = mean(anastat(this_allu_ok,:,posstmi),2);
                    y = mean(anastat(this_allu_ok,:,posstmi+1),2);
                    xp = combine_pvalues(anastatp(this_allu_ok,:,posstmi),2);
                    yp = combine_pvalues(anastatp(this_allu_ok,:,posstmi+1),2);
                    if ~isempty(regexp(statmode{stmi},'roc', 'once'))
                        maxlim = {[0 1],[0 1]};
                    else
                        tmp = [max(abs(x))*[-1 1],max(abs(y))*[-1 1]];
                        maxlim = {tmp,tmp};
                    end
        
                    figure('Visible','on');
                    myScatterWithHistogram(xylabels{1},xylabels{2},x,y,xp,yp,[],[],[],maxlim,x0,y0,[],sigthresh,1);
        
                    figname = 'supp_figure_8a';

                    savepath=fullfile(neuralsavedir,figname);
                    saveas(gcf,[savepath '.fig']);
                    saveas(gcf,[savepath '.png']);
                    close(gcf);
                end
            end
        end
        
        if do_plot
            if ~do255025only
                %% info and noinfo RPE is the same, justify combining epochs
                basettl = 'rpe_info_vs_noinfo';
                basettl = [monktxt ', Offer ' num2str(offi) ' conds, ' basettl];
        
                x = anastat(this_allu_ok,1,stmi);
                y = anastat(this_allu_ok,2,stmi);
                xp = anastatp(this_allu_ok,1,stmi);
                yp = anastatp(this_allu_ok,2,stmi);
                if ~isempty(regexp(statmode{stmi},'roc', 'once'))
                    maxlim = {[0 1],[0 1]};
                else
                    tmp = {[median([min(x),max(x)])+1.05*range(x)*[-1/2 1/2]],[median([min(y),max(y)])+1.05*range(y)*[-1/2 1/2]]};
                    tmp = [min(tmp{1}(1),tmp{2}(1)),max(tmp{1}(2),tmp{2}(2))];
                    maxlim = {tmp,tmp};
                end
        
                figure('Visible','on');
                myScatterWithHistogram([statmode{stmi} ', Info cue'],[statmode{stmi} ', Noinfo rev'],x,y,xp,yp,[],[],[],maxlim,x0,y0,[],sigthresh,1);
        
                figname = 'supp_figure_8b';
                savepath=fullfile(neuralsavedir,figname);
                saveas(gcf,[savepath '.fig']);
                saveas(gcf,[savepath '.png']);
                close(gcf);
            end
        end
    end
end
%%
if 1
    %% formal glm 
    rpe_efits = {}; 
    apis = 3:4; model_i = 3; half_i = 3; wi=1; rpe_bis = 1:2; 

    for apii=1:length(apis)
        rpe_efits{apii}={};
        for mi_=1:2
            rpe_efits{apii} = vertcat(rpe_efits{apii},glm_analysis_periods{mi_}(apis(apii)).efits{:,wi,model_i,half_i});
        end
        rpe_efits{apii}=vertcat(rpe_efits{apii}{:});
    end

    savedir = neuralsavedir;
    x = {}; y = x; xp = x; yp = x; x_rpe = x; xp_rpe = x; y_rpe = x; yp_rpe = x;
    for apii=1:2
        psname = glm_analysis_periods{1}(apis(apii)).ps;
        bs = [rpe_efits{apii}.b]';
        ss = [rpe_efits{apii}.stats];
        ps = [ss.p]';
        bnames = rpe_efits{apii}(1).xname;
        rpe_bnames = bnames(rpe_bis);
        rpe_bs = bs(:,rpe_bis);
        rpe_ps = ps(:,rpe_bis);
        xname_ = rpe_bnames{1}; yname_ = rpe_bnames{2};
        
        x{apii} = NaN(sum(allu_ok),1); y{apii} = x{apii};
        xp{apii} = x{apii}; yp{apii} = x{apii};

        x_rpe{apii} = NaN(size(allu_ok)); y_rpe{apii} = x_rpe{apii};
        xp_rpe{apii} = x_rpe{apii}; yp_rpe{apii} = x_rpe{apii};
        
        x{apii} = rpe_bs(:,1); y{apii} = rpe_bs(:,2);
        xp{apii} = rpe_ps(:,1); yp{apii} = rpe_ps(:,2);

        x_rpe{apii}(allu_ok) = rpe_bs(:,1); y_rpe{apii}(allu_ok) = rpe_bs(:,2);
        xp_rpe{apii}(allu_ok) = rpe_ps(:,1); yp_rpe{apii}(allu_ok) = rpe_ps(:,2);
    end
    
    figure(); 
    myScatterWithHistogram(escape(xname_),escape(yname_),mean([x{:}],2),mean([y{:}],2),combine_pvalues([xp{:}],2),combine_pvalues([yp{:}],2),[],[],[],'seplims',[],[],[],sigthresh,1);
    set(gcf,'Position',[1 1 800 800]);
    figname = 'supp_figure_9a';

    savepath=fullfile(neuralsavedir,figname);
    saveas(gcf,[savepath '.fig']);
    saveas(gcf,[savepath '.png']);
    close(gcf);

    %% trial by trial ER/RPE noise corr
    off1modeli = 1; rpemodeli = 3;
    off1_rpe_noisecorr = NaN(2,2,length(allu_ok),2);
    ns = NaN(size(allu_ok));
    for wimode = 1:2
        switch wimode
            case 1
                off1wi = 1; inforpewi = 1; noinforpewi = 1; 
            case 2
                off1wi = 2; inforpewi = 2; noinforpewi = 3; 
        end
    
        off1efits = vertcat(glm_analysis_periods{1}(1).efits(:,off1wi,off1modeli,3), glm_analysis_periods{2}(1).efits(:,off1wi,off1modeli,3));
        inforpeefits = vertcat(glm_analysis_periods{1}(3).efits(:,inforpewi,rpemodeli,3), glm_analysis_periods{2}(3).efits(:,inforpewi,rpemodeli,3));
        noinforpeefits = vertcat(glm_analysis_periods{1}(4).efits(:,noinforpewi,rpemodeli,3), glm_analysis_periods{2}(4).efits(:,noinforpewi,rpemodeli,3));
        
        for ui=1:length(allu_ok)
            if ~allu_ok(ui)
                continue;
            end
            if ui<=224
                u = ui;
                mi_ = 1;
            else
                u = ui-224;
                mi_ = 2;
            end

            oktrall = ps_analysis_periods{mi_}.oktrall(:,u);
            
            switch off1wi
                case 1
                    off1resid = off1efits{ui}.stats.resid;
                case 2
                    off1resid = off1efits{ui}.y;
            end
    
            rperesid = NaN(sum(oktrall),1);
            rperesid(chose_info{mi_}(oktrall)) = inforpeefits{ui}.stats.resid;
            rperesid(~chose_info{mi_}(oktrall)) = noinforpeefits{ui}.stats.resid;
        
            [off1_rpe_noisecorr(1,1,ui,wimode),off1_rpe_noisecorr(1,2,ui,wimode)] = corr(rperesid(chose_off1{mi_}(oktrall)),off1resid(chose_off1{mi_}(oktrall)),'type','Spearman');
            [off1_rpe_noisecorr(2,1,ui,wimode),off1_rpe_noisecorr(2,2,ui,wimode)] = corr(rperesid(chose_off1{mi_}(oktrall)),off1resid(chose_off1{mi_}(oktrall)),'type','Pearson');
            ns(ui) = sum(chose_off1{mi_}(oktrall));
        end
    end
    
    % see https://onlinelibrary.wiley.com/doi/10.1002/0471667196.ess5050
    % Treating the Spearman coefficients as though they were Pearson coefficients and using the standard Fisher's z-transformation and subsequent comparison was more robust with respect to Type I error than either ignoring the nonnormality and computing Pearson coefficients or converting the Spearman coefficients to Pearson equivalents prior to transformation.
    p_compare_cc = NaN(length(allu_ok),2);
    for corrmode=1:2
        for ui=1:length(allu_ok)
            p_compare_cc(ui,corrmode) = compare_correlation_coefficients(squeeze(off1_rpe_noisecorr(corrmode,1,ui,1)),squeeze(off1_rpe_noisecorr(corrmode,1,ui,2)),ns(ui),ns(ui));
        end
    end
    
    figure; set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.6 0.9]); tcl = tiledlayout(1,2);
    for corrmode=1:2
        switch corrmode
            case 1
                cmt = 'Spearman';
            case 2
                cmt = 'Pearson';
        end
        condi=1;
        condname = ertr_conds(condi).name;
        nexttile();
        mySigHistogram(diff(squeeze(off1_rpe_noisecorr(corrmode,1,ertr_conds(condi).thisalluok,:)),1,2),p_compare_cc(ertr_conds(condi).thisalluok,corrmode),sigthresh,[],[],[],[],[],{'Difference in corr(Offer 1 E[R] resid, RPE resid)', 'for analysis window (post-event) vs control window (pre-event)'},'neurons');
        axis square;
        title([cmt ' correlation , ' condname ' neurons']);
    end
    
    figname = 'supp_figure_9e';

    savepath=fullfile(neuralsavedir,figname);
    saveas(gcf,[savepath '.fig']);
    saveas(gcf,[savepath '.png']);
    close(gcf);
end

if 1
    %%
    combinedR = mean([x{:}],2); combinedRp = combine_pvalues([xp{:}],2); sigcombinedR = combinedR; sigcombinedR(combinedRp>=sigthresh) = 0;
    combinedrew_ul = mean([y{:}],2); combinedrew_ulp = combine_pvalues([yp{:}],2); sigcombinedrew_ul = combinedrew_ul; sigcombinedrew_ul(combinedrew_ulp>=sigthresh) = 0;
    
    % how many in each quadrant?
    % how many significant in each quadrant?
    
    allT = [sum(combinedR < 0 & combinedrew_ul > 0) sum(combinedR > 0 & combinedrew_ul > 0); sum(combinedR < 0 & combinedrew_ul < 0) sum(combinedR > 0 & combinedrew_ul < 0)];
    sigT = [sum(sigcombinedR < 0 & sigcombinedrew_ul > 0) sum(sigcombinedR == 0 & sigcombinedrew_ul > 0) sum(sigcombinedR > 0 & sigcombinedrew_ul > 0); sum(sigcombinedR < 0 & sigcombinedrew_ul == 0) sum(sigcombinedR == 0 & sigcombinedrew_ul == 0) sum(sigcombinedR > 0 & sigcombinedrew_ul == 0); sum(sigcombinedR < 0 & sigcombinedrew_ul < 0) sum(sigcombinedR == 0 & sigcombinedrew_ul < 0) sum(sigcombinedR > 0 & sigcombinedrew_ul < 0)];
    RsigT = cat(3,[sum(sigcombinedR < 0 & sigcombinedrew_ul > 0 & relevant_ps(allu_ok,1)>=sigthresh) sum(sigcombinedR == 0 & sigcombinedrew_ul > 0 & relevant_ps(allu_ok,1)>=sigthresh) sum(sigcombinedR > 0 & sigcombinedrew_ul > 0 & relevant_ps(allu_ok,1)>=sigthresh); sum(sigcombinedR < 0 & sigcombinedrew_ul == 0 & relevant_ps(allu_ok,1)>=sigthresh) sum(sigcombinedR == 0 & sigcombinedrew_ul == 0 & relevant_ps(allu_ok,1)>=sigthresh) sum(sigcombinedR > 0 & sigcombinedrew_ul == 0 & relevant_ps(allu_ok,1)>=sigthresh); sum(sigcombinedR < 0 & sigcombinedrew_ul < 0 & relevant_ps(allu_ok,1)>=sigthresh) sum(sigcombinedR == 0 & sigcombinedrew_ul < 0 & relevant_ps(allu_ok,1)>=sigthresh) sum(sigcombinedR > 0 & sigcombinedrew_ul < 0 & relevant_ps(allu_ok,1)>=sigthresh)],[sum(sigcombinedR < 0 & sigcombinedrew_ul > 0 & relevant_ps(allu_ok,1)<sigthresh) sum(sigcombinedR == 0 & sigcombinedrew_ul > 0 & relevant_ps(allu_ok,1)<sigthresh) sum(sigcombinedR > 0 & sigcombinedrew_ul > 0 & relevant_ps(allu_ok,1)<sigthresh); sum(sigcombinedR < 0 & sigcombinedrew_ul == 0 & relevant_ps(allu_ok,1)<sigthresh) sum(sigcombinedR == 0 & sigcombinedrew_ul == 0 & relevant_ps(allu_ok,1)<sigthresh) sum(sigcombinedR > 0 & sigcombinedrew_ul == 0 & relevant_ps(allu_ok,1)<sigthresh); sum(sigcombinedR < 0 & sigcombinedrew_ul < 0 & relevant_ps(allu_ok,1)<sigthresh) sum(sigcombinedR == 0 & sigcombinedrew_ul < 0 & relevant_ps(allu_ok,1)<sigthresh) sum(sigcombinedR > 0 & sigcombinedrew_ul < 0 & relevant_ps(allu_ok,1)<sigthresh)]);
    
    ORx = NaN(size(sigT));
    fisherp = ORx;
    for ii=1:size(ORx,1)
        for jj=1:size(ORx,2)
            a = [RsigT(ii,jj,2) (sum(RsigT(:,:,2),'all')-RsigT(ii,jj,2));RsigT(ii,jj,1) (sum(RsigT(:,:,1),'all')-RsigT(ii,jj,1))];
            ORx(ii,jj) = (a(1,1)/a(1,2))/(a(2,1)/a(2,2));
            [~,fisherp(ii,jj)] = fishertest(a);
        end
    end
    
    %%
    figure; set(gcf,'Position',get(0,'ScreenSize').*[1 1 1 1]);
    tiledlayout(2,2);
    nexttile;
    myheatmap(RsigT(:,:,2)/sum(RsigT(:,:,2),'all'));
    xticks(1:3); yticks(1:3);
    xticklabels({'Expected GLM weight < 0','Expected GLM weight not sig.','Expected GLM weight > 0'}); yticklabels({'Delivered GLM weight > 0','Delivered GLM weight not sig.','Delivered GLM weight < 0'});
    title('Offer E[R]-signaling neurons'); 
    nexttile;
    myheatmap(RsigT(:,:,1)/sum(RsigT(:,:,1),'all'));
    xticks(1:3); yticks(1:3);
    xticklabels({'Expected GLM weight < 0','Expected GLM weight not sig.','Expected GLM weight > 0'}); yticklabels({'Delivered GLM weight > 0','Delivered GLM weight not sig.','Delivered GLM weight < 0'});
    title('Not Offer E[R]-signaling neurons'); 
    nexttile;
    myheatmap(log10(ORx),fisherp,[-1,1],20);
    xticks(1:3); yticks(1:3);
    xticklabels({'Expected GLM weight < 0','Expected GLM weight not sig.','Expected GLM weight > 0'}); yticklabels({'Delivered GLM weight > 0','Delivered GLM weight not sig.','Delivered GLM weight < 0'});
    title('Log10 Odds Ratio'); 

    figname = 'supp_figure_9cd';

    savepath=fullfile(neuralsavedir,figname);
    saveas(gcf,[savepath '.fig']);
    saveas(gcf,[savepath '.png']);
    close(gcf);
    
    %%
    combinedR = mean([x_rpe{:}],2); combinedRp = combine_pvalues([xp_rpe{:}],2); sigcombinedR = combinedR; sigcombinedR(combinedRp>=sigthresh) = 0;
    combinedrew_ul = mean([y_rpe{:}],2); combinedrew_ulp = combine_pvalues([yp_rpe{:}],2); sigcombinedrew_ul = combinedrew_ul; sigcombinedrew_ul(combinedrew_ulp>=sigthresh) = 0;
    fontsize = 5;
    barheightp = 0.03;
    formalGLMRPE = sigcombinedR < 0 & sigcombinedrew_ul > 0; 
    formalGLMoutcome = sigcombinedR == 0 & sigcombinedrew_ul > 0; 
    formalGLMoutcome_veryweakpred_p = combinedRp>0.3 & sigcombinedrew_ul > 0;
    rew_ul_R_theta = atan(combinedR./combinedrew_ul)/pi*180;
    formalGLMoutcome_veryweakpred_ratio = sigcombinedR == 0 & sigcombinedrew_ul > 0 & abs(rew_ul_R_theta)<10;
    formalGLMnothing = sigcombinedR == 0 & sigcombinedrew_ul == 0; 
    ERneuron = relevant_ps(:,strcmp(xname,'R'))<sigthresh;
    
    %%
    uconds = struct('name',{'Both +Delivered and -Predicted','Mainly +Delivered'},'thisalluok',{formalGLMRPE & allu_ok,formalGLMoutcome_veryweakpred_ratio & allu_ok});

    h = {};
    maxylims = zeros(length(uconds),2);

    figure; tiledlayout(length(uconds),length(revoutplotconds),'Padding','compact','TileSpacing','compact'); set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.7 1]);
    for uci=1:length(uconds)
        for pci=1:length(revoutplotconds)
            meanfullps = mean(revoutplotconds(pci).fullps(uconds(uci).thisalluok,:,:),1);
            sefullps = std(revoutplotconds(pci).fullps(uconds(uci).thisalluok,:,:),[],1)/sqrt(sum(uconds(uci).thisalluok));
            
            h{uci,pci} = nexttile();
            for pcii=revoutplotconds(pci).condi
                shadedErrorBar(revoutplotconds(pci).t,meanfullps(1,:,pcii),sefullps(1,:,pcii),'lineProps',{'color',revoutplotconds(pci).condcolors{pcii},'linewidth',3,'linestyle','-'});
            end
            
            maxylims(uci,:) = [min([maxylims(uci,1) min(ylim)]) max([maxylims(uci,2) max(ylim)])];
        end
    end
    for uci=1:length(uconds)
        for pci=1:length(revoutplotconds)
            axes(h{uci,pci});
            hold on;
            ylim(maxylims(uci,:));
            xlim(revoutplotconds(pci).xlims);
            axis square;
            linex(0,'k','HandleVisibility','off');

            if pci==1
                text(max(xlim)-0.95*range(xlim),max(ylim)-0.05*range(ylim),{uconds(uci).name,['(n=' num2str(sum(uconds(uci).thisalluok)) ')']},'FontSize',12,'HorizontalAlignment','left','VerticalAlignment','top');
                ylabel('Normalized activity');
            else
                yticklabels([]);
            end
            if uci==1
                legend(revoutplotconds(pci).condnames(revoutplotconds(pci).condi),'Location','northeast');
                title(revoutplotconds(pci).name);
                xticklabels([]);
            else
                xticklabels(xticks/1000);
                xlabel(['Time from ' revoutplotconds(pci).xmarkname ' onset (s)']);
            end
        end
    end

    figname = 'supp_figure_9b';

    savepath=fullfile(neuralsavedir,figname);
    saveas(gcf,[savepath '.fig']);
    saveas(gcf,[savepath '.png']);
    close(gcf);
end
