%% analog_analysis.m
analog_types = {'eyex','eyey','pupil','licking'};
minmaxthresh = 4900; 
t = -2000:2000;
nanalog = 4;
ps = struct('name',{'fix','off1','off2','cstart','cue','rev','out','stimoff'},'code',num2cell(14000+[2 4 5 8 8 8 2000-14000 9]));
first_trial_search_code_seq = [3800 4001 4010 4021];
infonesses = [1 0 2];
dofake = 0; 
psnames = {'cue','rev','out'}; psis = find(ismember({ps.name},psnames));
rpesigns = [1 0 -1];
rpecolors = [1 0 0; 0.5 0.5 0.5; 0 0 1];
nqs = length(rpesigns);
do255025only = 1;

makePlots(monkeydef,infonesses,do255025only,rpesigns,rpecolors,dofake,nqs,psis,ps,t,sigthreshsavedir);

%% helper functions
function makePlots(monkeydef,infonesses,do255025only,rpesigns,rpecolors,dofake,nqs,psis,ps,t,savedir)
    mi_ = 0;
    monktxt = 'both_monkeys';
        
    %% sac_aligned

    thesestats = cateachcell(1,monkeydef(1).sac_aligned,monkeydef(2).sac_aligned);
    pupylimscale = 2;
    stattypei=2;
    statwini=1;

    z = cell(1,1,3);
    for rpei=1:3
        z{:,:,rpei} = vertcat(thesestats{2,stattypei+2*(statwini-1),rpei},thesestats{1,stattypei+2*(statwini-1),rpei});
    end
    z = squeeze(z);
    zn = cellfun(@(x)sum(~isnan(x)),z);

    statttl = ['pupil before ' ternary(statwini==1,'1st','2nd') ' sac'];
    ylbl = {'zscored and base subtr', '(100ms prior to event)'};
    ylims = pupylimscale*1*[-1 1];

    infotxt = 'info and noinfo';


    figure('Visible','on'); set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.8 0.8]);
    rsps = NaN(nqs,1); rspgroups = cell(size(rsps));
    for qi=1:nqs
        qii = mod(qi,nqs)+1;
        rsps(qi) = ranksum(z{qi},z{qii});
        rspgroups{qi} = [qi qii];
    end

    zm = cellfun(@(x)nanmedian(x),z);
    zerr = 1.25/1.35*cellfun(@(x)diff(prctile(x,[25 75])),z)./sqrt(zn);
    figtypetxt = 'dotmedian';

    plot(1:length(zm),zm,'-','Color','k');
    hold on;
    errorbar(1:length(zm),zm,zerr,'HandleVisibility','off','color','k','CapSize',0,'LineWidth',2,'LineStyle','none');
    for qi=1:nqs
        plot(qi,zm(qi),"Marker",'o','MarkerFaceColor',rpecolors(qi,:),'MarkerEdgeColor','none','MarkerSize',10);
    end

    sigstar(rspgroups,rsps,[],[],[],[],1);
    xticks(1:nqs);
    xticklabels({['RPE+ (n=' num2str(zn(1)) ')'],['RPE0 (n=' num2str(zn(2)) ')'],['RPE- (n=' num2str(zn(3)) ')']})
    set( gca, 'xdir', 'reverse' )
    ylabel(ylbl);
    title({infotxt,statttl});
    axis square;
    ttl = [monktxt infotxt statttl ' feedback_and_sac_triggered' ternary(do255025only,' do255025only','') figtypetxt];
    title(ttl);
    figname = 'supp_figure_7a_left';
    savepath=fullfile(savedir,figname);
    saveas(gcf,[savepath '.fig']);
    saveas(gcf,[savepath '.png']);
    close(gcf);

    %% ps_aligned

    thesestats = cateachcell(1,monkeydef(1).ps_aligned,monkeydef(2).ps_aligned);
    monkn = {cellfun(@(x)size(x,1),monkeydef(1).ps_aligned),cellfun(@(x)size(x,1),monkeydef(2).ps_aligned)};

    boxplotstats = cell(size(1:4)); boxplotgroups = boxplotstats;
    for pltmodei=3:4
        switch pltmodei
            case 3
                ttl = 'p gaze at chosen';
                figname = 'supp_figure_7a_middle';
                ylbl = ttl;
                ylims = [];
                    boxplotanawind = [300 800];

                    boxpsis = [6 5]; %noinfo info
                doany = false;
            case 4
                ttl = 'p gaze at unchosen';
                figname = 'supp_figure_7a_right';
                ylims = [];
                ylbl = ttl;
                boxplotanawind = [-150 650];
                boxpsis = [7 7]; %noinfo info
                doany = true;
        end

        ttl = [ttl ternary(do255025only,' do255025only','')];

        for infoness_=1:2
            infoness = infonesses(infoness_);
            switch infoness
                case 0
                    infotxt = 'noinfo';
                case 1
                    if dofake
                        infotxt = 'noinfo';
                    else
                        infotxt = 'info';
                    end
                case 2
                    infotxt = 'info and noinfo';
            end
            for psi_=1:length(psis) 
                psi = psis(psi_);
                for rpei=1:length(rpesigns)
                    z = thesestats{pltmodei,infoness_,psi_,rpei};
                    if ~isempty(boxplotanawind) && psi==boxpsis(infoness+1)
                        if doany

                            tmpcpcutoffs = vertcat(monkeydef(1).cpcutoff*ones(monkn{1}(pltmodei,infoness_,psi_,rpei),1),monkeydef(2).cpcutoff*ones(monkn{2}(pltmodei,infoness_,psi_,rpei),1));

                            boxplotstat = nanmean(z(:,inbounds(t,boxplotanawind)),2) > tmpcpcutoffs;

                        else
                            boxplotstat = nanmean(z(:,inbounds(t,boxplotanawind)),2); 
                        end
                        boxplotstats{pltmodei} = vertcat(boxplotstats{pltmodei},boxplotstat);
                        boxplotgroup = rpei*ones(size(boxplotstat));
                        boxplotgroups{pltmodei} = vertcat(boxplotgroups{pltmodei},boxplotgroup);
                    end
                end
            end
        end

        %% box/bar plots
        if ~isempty(boxplotanawind)
            ttl = [ttl num2str(boxplotanawind)];
            h2 = figure('Visible','on'); set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.8 0.8]);
            axis square;
            hold on;
            rsps = NaN(nqs,1); rspgroups = cell(size(rsps)); zn = rsps;
            tmpmeans = [];
            for qi=1:nqs
                qii = mod(qi,3)+1;
                zn(qi) = sum(~isnan(boxplotgroups{pltmodei}) & boxplotgroups{pltmodei} == qi);
                if doany
                    pd = fitdist(boxplotstats{pltmodei}(boxplotgroups{pltmodei} == qi),'Binomial');
                    tmp = paramci(pd,'Alpha',0.32);
                    ci68 = tmp(:,2);
                    tmpmean = pd.mean;
                    negerr = tmpmean-min(ci68);
                    poserr = max(ci68)-tmpmean;

                    rsps(qi) = permutation_exact_test_of_proportions(nansum(boxplotstats{pltmodei}(boxplotgroups{pltmodei} == qi)),sum(~isnan(boxplotstats{pltmodei}(boxplotgroups{pltmodei} == qi))),nansum(boxplotstats{pltmodei}(boxplotgroups{pltmodei} == qii)),sum(~isnan(boxplotstats{pltmodei}(boxplotgroups{pltmodei} == qii))));
                    rspgroups{qi} = [qi qii];
                    figtypetxt = 'dot';
                else
                    tmpy = boxplotstats{pltmodei}(boxplotgroups{pltmodei} == qi); tmpyy = boxplotstats{pltmodei}(boxplotgroups{pltmodei} == qii);
                    rsps(qi) = ranksum(tmpy,tmpyy);
                    rspgroups{qi} = [qi qii];

                    figtypetxt = 'dotmedian';
                    tmpmean = nanmedian(tmpy);% there was a weird bug here??? nanstd
                    tmperr = 1.25/1.35*diff(prctile(tmpy,[25 75]))./sqrt(sum(~isnan(tmpy)));
                    negerr = tmperr; poserr = tmperr;

                end
                tmpmeans(end+1) = tmpmean;

                errorbar(qi,tmpmean,negerr,poserr,'HandleVisibility','off','color','k','CapSize',0,'LineWidth',2,'LineStyle','none');
                plot(qi,tmpmean,"Marker",'o','MarkerFaceColor',rpecolors(qi,:),'MarkerEdgeColor','none','MarkerSize',10);

            end

            tmp = plot(1:nqs,tmpmeans,'-','Color','k');
            uistack(tmp,'bottom');

            if ~isempty(ylims)
                ylim(ylims);
            end
            xticks(1:nqs);
            xticklabels({['RPE+ (n=' num2str(zn(1)) ')'],['RPE0 (n=' num2str(zn(2)) ')'],['RPE- (n=' num2str(zn(3)) ')']})
            sigstar(rspgroups,rsps,[],[],[],[],1);
            set( gca, 'xdir', 'reverse' )
            thisttl = [monktxt ' ' infotxt ' ' ttl ' ' figtypetxt];
            title(thisttl);
            axis square;

            savepath=fullfile(savedir,figname);
            saveas(gcf,[savepath '.fig']);
            saveas(gcf,[savepath '.png']);
            close(gcf);
        end
    end
end