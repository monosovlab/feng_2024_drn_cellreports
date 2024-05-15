function plot_pref_comparison(allu_ok,x,color,comparisonstat,comparisonstatwin,vstext,anawind,comparisonstatmode,psname,gauss_smoothing_sigma,xlims)
    u_ok = allu_ok & any(~isnan(comparisonstat),2);
    flippedcomparisonstat = comparisonstat(u_ok,:);
    flippedcomparisonstatwin = comparisonstatwin(u_ok);

    if strcmp(comparisonstatmode,'ROC')
        flippedcomparisonstatwin = flippedcomparisonstatwin-0.5;
        y0 = 0.5;
    else
        y0 = 0;
    end

    psn = flippedcomparisonstat;
    if size(psn,1)
        plot_psn(x,psn,color,gauss_smoothing_sigma);
        p = signrank(flippedcomparisonstatwin);
    else
        p = NaN;
        warning([' no cells']);
    end
    linex(0,'k'); liney(y0,'k');
    xlim(xlims);
    set(gca,'FontSize',20);
    xlabel(['Time from ' psname ' onset (ms)']);
    ylabel(comparisonstatmode);
    text(max(xlim)-0.05*range(xlim),max(ylim)-0.95*range(ylim),['n=' num2str(round(size(psn,1),2))],'FontSize',20,'HorizontalAlignment','right','VerticalAlignment','bottom');
    text(max(xlim)-0.95*range(xlim),max(ylim)-0.05*range(ylim),vstext,'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','top');
    text(200,min(ylim)+0.15*range(ylim),['p=' num2str(round(p,3))],'FontSize',20);
    rectangle('Position',[min(anawind) min(ylim) range(anawind) 0.05*range(ylim)],'FaceColor',[0 0 0]);
end

function flippedstat = flipstat(stat,statmode)
    if strcmp(statmode,'ROC')
        flippedstat = 1-stat;
    elseif strcmp(statmode,'diff')
        flippedstat = -stat;
    else
        error('undefined');
    end
end