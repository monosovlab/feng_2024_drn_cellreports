savedir = sigthreshsavedir;
basettl='choice behavior GLM from recording';
soi=1;
tmpmodes = {'reduced','IxS'};
for gi=1:length(monkeyabbr)
    for tmpi=1:length(tmpmodes)
        ttl = [basettl ', ' tmpmodes{tmpi}];
        % select regs to show
        selxnames = {'R','O','S','IxS'};
        if strcmp(tmpmodes{tmpi},'reduced')
            selxnames = selxnames([1:3]);
            figname = ['figure_1c_' char(monkeyabbr(gi))];
        else
            figname = ['supp_figure_4b_' char(monkeyabbr(gi))];
        end
    
        selxnamelabels = get_xnamelabels(selxnames);
        [~,~,attridx]=intersect(selxnames,befits{soi,gi}.xname,'stable');
        attridx = attridx';
    
        ttl=[ttl ', ' ternary(soi==1,'offer_diff','sep_offers')];
        ttl=[char(monkeyabbr(gi)) ' ' ttl];
    
        figure('Visible','on');
        set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.4 0.6]);
        if soi==1
            bar(befits{soi,gi}.b(attridx),'FaceColor',0.5*[1 1 1],'LineWidth',2);
            hold on;
            errorbar(befits{soi,gi}.b(attridx),befits{soi,gi}.stats.se(attridx),'LineStyle','none','Color','k','CapSize',0,'LineWidth',2);
            xts = 1:length(attridx);
            xticks(xts);
            xticklabels(selxnamelabels);
            sigstar(xts,befits{soi,gi}.stats.p(attridx),[],sigthresh);
        end
    
        ylim([-1.5 1.5]);
    
        set(gca,'FontSize',20);
        title(escape(ttl));
        ylabel({['Effect on ' ternary(soi==1,'subjective value','choosing offer 2')],'(units: log odds of biasing choice)'});
        xlabel('Regressor');
        text(max(xlim)-0.95*range(xlim),max(ylim)-0.05*range(ylim),['Monkey ' char(monkeyabbr(gi))],'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','top');
        text(max(xlim)-0.05*range(xlim),max(ylim)-0.95*range(ylim),[num2str(befits{soi,gi}.n) ' trials'],'FontSize',20,'HorizontalAlignment','right','VerticalAlignment','bottom');

        savepath=fullfile(savedir, figname);
        saveas(gcf,[savepath '.fig']);
        saveas(gcf,[savepath '.png']);
        close(gcf);
    end
end
