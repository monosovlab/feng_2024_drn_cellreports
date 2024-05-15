%% extra_behavior.m
scriptname = 'extrabehavior';

altmodis = [1:6];
nboot = 2000;
if 1
    %% compare halves
    plotxnames = {'R','O','S'}; plotregcolors = [1 0 0; 0 0 1; 0 0.7 0;];
    
    figure; set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.8 0.8]);
    tiledlayout(2,length(plotxnames)); zp = NaN(2,length(plotxnames));
    for mi_=1:2
        for pri=1:length(plotxnames)
            nexttile;
            for hi=1:2
                xsel = strcmp(halfbefits{mi_,hi}.xname,plotxnames{pri});
                bar(hi,halfbefits{mi_,hi}.b(xsel),'FaceColor',plotregcolors(pri,:)/hi,'BarWidth',0.7);
                hold on;
                errorbar(hi,halfbefits{mi_,hi}.b(xsel),halfbefits{mi_,hi}.stats.se(xsel),'HandleVisibility','off','color','k','CapSize',0,'LineWidth',2);
            end
            z = abs(halfbefits{mi_,1}.b(xsel)-halfbefits{mi_,2}.b(xsel))/norm([halfbefits{mi_,1}.stats.se(xsel) halfbefits{mi_,2}.stats.se(xsel)]);
            zp(mi_,pri) = 2*normcdf(-z);
            
            xlim([0 3]); xticks(1:2); xticklabels({'1st','2nd'}); 

            ylim([-1.2 1.2]*max(abs(ylim)));
            axis square;

            xlabel('Halves of sessions');
            ylabel('{\beta}_{preference} GLM weight');
            title(['Monkey ' monkeyabbr(mi_) ' ' get_xnamelabels(plotxnames{pri})]);
        end
    end
    tmp = 1;
    for mi_=1:2
        for pri=1:length(plotxnames)
            nexttile(tmp);
            [~, ~, ~, adj_p] = fdr_bh(zp);
            sigstar({[1,2]},adj_p(mi_,pri),[],[],[],[],1);
            axis square;
            tmp = tmp+1;
        end
    end

    ttl = 'comparehalvesofsession';
    figname = 'supp_figure_2f';
    savepath=fullfile(sigthreshsavedir,figname);
    saveas(gcf,[savepath '.fig']);
    saveas(gcf,[savepath '.png']);
    close(gcf);
    
    %% can do per session behavior trends here, just get ses from above loop

    thissessionbs = sessionbs;
    thissessionses = sessionses;
    
    figure; set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.8 0.8]); tiledlayout(1,2);
    for mi_=1:2
        [uniquethissessionbs, tmp] = unique(thissessionbs{mi_},'rows','stable');
        uniquethissessionses = thissessionses{mi_}(tmp,:);
        nexttile;
        for pri=1:length(plotxnames)
            xsel = strcmp(befits{1,mi_}.xname,plotxnames{pri});
            plot(uniquethissessionbs(:,xsel),'color',plotregcolors(pri,:));
            hold on;
            errorbar(uniquethissessionbs(:,xsel),uniquethissessionses(:,xsel),'HandleVisibility','off','color',plotregcolors(pri,:));
        end
        liney(0,'k');
        axis square;
        xticks(1:size(uniquethissessionbs,1));
        xlabel('Session #');
        ylabel('{\beta}_{preference} GLM weight');
        title(['Monkey ' monkeyabbr(mi_)]);
        legend(get_xnamelabels(plotxnames),'Location','northwest');
    end
    ttl = 'sessionbysessiontrend';
    figname = 'supp_figure_2e';

    savepath=fullfile(sigthreshsavedir, figname);
    saveas(gcf,[savepath '.fig']);
    saveas(gcf,[savepath '.png']);
    close(gcf);
end

if 1
    %% behavior model comparison
    for ami=1:length(altmodis)
        for mi_=1:length(monkeyabbr)
            bt = altbefits{mi_,ami};
            betas{mi_,ami} = bt.b;
            nparams(mi_,ami) = length(betas{mi_,ami});
            logliks(mi_,ami) = bt.loglik;
            cvlogliks(mi_,ami) = bt.cvloglik;
            bootlogliks{mi_,ami} = bt.bootloglik;
            bootcvlogliks{mi_,ami} = bt.bootcvloglik;
            outputs{mi_,ami} = bt.output;
        end
    end
    tmp = [befits{1,:}];
    aic = 2*(-logliks + nparams);
    bic = -2*logliks + nparams.*log([tmp.n]');
    keyweights = NaN(size(nparams));
    for ami=1:length(altmodis)
        for mi_=1:length(monkeyabbr)
            switch ami
                case 1
                    switch mi_
                        case 1
                            tmp = 10;
                        case 2
                            tmp = 6;
                    end
                case {2,3}
                    tmp = 1;
                case 4
                    tmp = 2;
                case {5,6}
                    tmp = 3;
            end
            keyweights(mi_,ami) = betas{mi_,ami}(tmp);
        end
    end
    
    Tdef = struct(...
        'name',{'AIC','BIC','# of parameters'},'mat',{aic,bic,nparams});
    
    vars = {'Linear','Hyp. discounting','Exp. discounting','Mag:Power, Prob:Prelec-1','Mag:Power, Prob:Prelec-2','Mag:Power, Prob:Linear-in-log-odds'};
    rownames = {['Monkey ' monkeyabbr(1)],['Monkey ' monkeyabbr(2)]};
    
    for ti=1:length(Tdef)
        T = array2table(round(Tdef(ti).mat,2));
        T.Properties.VariableNames = vars;
        T.Properties.RowNames = rownames;
        disp(Tdef(ti).name);
        disp(T);
        figname = ['supp_table_' num2str(ti+3)];

        savepath=fullfile(sigthreshsavedir, [figname '.txt']);
        writetable(T,savepath,'Delimiter','\t','WriteRowNames',true);
    end
    
    T = array2table(keyweights);
    T.Properties.VariableNames = {'Linear (Tout beta)','Hyp. disc. (k)','Exp. disc. (delta)','Prelec-1 (mag power)','Prelec-2 (mag power)','Linear-in-log-odds (mag power)'};
    T.Properties.RowNames = {monkeyabbr(1),monkeyabbr(2)};
    disp('key weights');
    disp(T);
    
    %% bootstrap model comparioson histograms
    figure; set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.6 0.9]); tiledlayout(2,2);
    for mci=1:2
        for mi_=1:2
            x = [bootlogliks{mi_,1} bootlogliks{mi_,5}];
            switch mci
                case 1
                    mc = 'AIC'; 
                    z = 2*(-x+nparams(mi_,[1 5]));
                    z_real = diff(aic(mi_,[1 5]));
                case 2
                    mc = 'BIC'; 
                    z = -2*x + nparams(mi_,[1 5]).*log(befits{1,mi_}.n);
                    z_real = diff(bic(mi_,[1 5]));
            end
            z = diff(z,1,2);
            nexttile();
            histogram(z,'NumBins',20,'FaceColor',0.5*[1 1 1]); axis square;
            hold on;
            xline(z_real,'k--');
            if prod(xlim)<0
                xline(0,'k-');
            end
    
            ttl = ['Comparison of Prelec-2 vs Linear value models by ' mc];
            title(['Monkey ' monkeyabbr(mi_) ' ' ttl]);
            xlabel(['Prelec-2 - Linear, ' mc]);
            ylabel('# of bootstraps');
            pseudop = 2*min(mean([z -z]<0));
            ciexcluding0 = 100.*(1 - pseudop);
            text(max(xlim)-0.95*range(xlim),max(ylim)-0.05*range(ylim),sprintf('%4.1f%% CI excludes 0',ciexcluding0),'FontSize',10,'HorizontalAlignment','left','VerticalAlignment','top');
            text(max(xlim)-0.95*range(xlim),max(ylim)-0.15*range(ylim),[num2str(numel(z)) ' bootstraps'],'FontSize',10,'HorizontalAlignment','left','VerticalAlignment','top');
        end
    end
    ttl = 'prelec2vslinearbootstraps';
    figname = 'supp_figure_2a';

    savepath=fullfile(sigthreshsavedir,figname);
    saveas(gcf,[savepath '.fig']);
    saveas(gcf,[savepath '.png']);
    close(gcf);
    
    %%
    figure; set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.9 0.9]);
    tiledlayout(1,2);
    for mi_ = 1:2
        nexttile();
        switch mi_
            case 0
                ttl = 'both monks';
                x = vertcat([altbefits{1,1}.valpred(:)],[altbefits{2,1}.valpred(:)]);
                y = vertcat([altbefits{1,5}.valpred(:)],[altbefits{2,5}.valpred(:)]);
            case {1,2}
                ttl = ['Monkey ' monkeyabbr(mi_)];
                x = vertcat([altbefits{mi_,1}.valpred(:)]);
                y = vertcat([altbefits{mi_,5}.valpred(:)]);
        end
        myScatterWithHistogram('Linear value (log-odds)','Prelec-2 value (log-odds)',x,y,[],[],[],[],[],[],[],[],[],sigthresh,1);
        h = gca;
        h.Children(end-2).SizeData = 1;
        axis square;
        legend('off');
        title(ttl);
    end
    ttl = 'Value of offer based on Prelec-2 vs Linear model';
    figname = 'supp_figure_2b_value';

    savepath=fullfile(sigthreshsavedir,figname);
    saveas(gcf,[savepath '.fig']);
    saveas(gcf,[savepath '.png']);
    close(gcf);
    
    %%
    for mi_=1:2
        figure; set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.9 0.9]);
        stdrpe = out_rew_ul{mi_} - altbefits{mi_,1}.xreg{1}.terms{1};
        prelecrpe = getPrelecRPE(out_rew_ul{mi_},altbefits{mi_,1},altbefits{mi_,5},mi_==2);
        myScatterWithHistogram('RPE_{std}','RPE_{subjective}',stdrpe(~isnan(prelecrpe)),prelecrpe(~isnan(prelecrpe)),[],[],[],[],[],'seplims',[],[],[],sigthresh,1,[],[],[],[],[],[],1);
        h = gca;
        h.Children(end-2).SizeData = 0.5;
        axis square;
        ttl = ['Monkey ' monkeyabbr(mi_)];
        title(ttl);
        ttl = [ttl '_stdvsprelec2rpescatter'];
        figname = ['supp_figure_2b_rpe_' monkeyabbr(mi_)];

        savepath=fullfile(sigthreshsavedir,figname);
        saveas(gcf,[savepath '.fig']);
        saveas(gcf,[savepath '.png']);
        close(gcf);
    end
end


%% helper functions
function rpe = getPrelecRPE(out_rew_ul,linearbefit,prelec2befit,t_cue_is_fixed)
    nuisanceregs = {'Lft','Rgt','2'};
    ev = linearbefit.xreg{strcmp(linearbefit.xname,'R')}.terms{1};
    info = linearbefit.xreg{strcmp(linearbefit.xname,'IxS')}.terms{1};
    t_rev_s = linearbefit.xreg{strcmp(linearbefit.xname,'O')}.terms{1};
    if t_cue_is_fixed
        t_cue_s = 0.5*ones(size(info));
    else
        t_cue_s = t_rev_s - linearbefit.xreg{strcmp(linearbefit.xname,'TadvxS')}.terms{1};
    end
    
    twaited = t_cue_s.*info + t_rev_s.*info;
    rpe = out_rew_ul;
    pred = NaN(size(rpe)); 
    out = pred;

    isT = strcmp(prelec2befit.xname,'O');
    isnuisance = ismember(prelec2befit.xname,nuisanceregs);
    xb = prelec2befit.b((end-length(prelec2befit.xname)+1):end);
    tcorrection = xb(isT)*twaited;
    nuisancecorrection = squeeze(sum(xb(isnuisance)'.*prelec2befit.x_sep(:,isnuisance,:),2));
    pred = prelec2befit.valpred - tcorrection - nuisancecorrection;

    isT = strcmp(prelec2befit.xname,'O');
    isnuisance = ismember(prelec2befit.xname,nuisanceregs);
    xb = prelec2befit.b((end-length(prelec2befit.xname)+1):end);
    tcorrection = xb(isT)*twaited;
    nuisancecorrection = squeeze(sum(xb(isnuisance)'.*prelec2befit.x_sep(:,isnuisance,:),2));
    tmprpe = rpe - ev;
    tmprpe(abs(tmprpe)<100) = 0;
    tmp = sign(tmprpe)+2; tmp=permute(repmat(tmp,1,1,3),[1 3 2]);
    tmp2 = repmat([1 2 3],sum(prelec2befit.oktr),1,2);
    feedbackprobs = double(tmp==tmp2);
    distcorrection = squeeze(prelec2befit.nonlinutilfun(prelec2befit.probs,prelec2befit.mags,prelec2befit.b(1),prelec2befit.b(2),prelec2befit.b(3))-prelec2befit.nonlinutilfun(feedbackprobs,prelec2befit.mags,prelec2befit.b(1),prelec2befit.b(2),prelec2befit.b(3)));
    out = prelec2befit.valpred - tcorrection - distcorrection - nuisancecorrection;

    rpe = out - pred;
end
