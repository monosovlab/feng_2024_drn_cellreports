function plot_pref_vs_nonpref(vstext,allu_ok,psns,x,colors,psname,gauss_smoothing_sigma,doerrorbar,xlims)
    if ~exist('doerrorbar','var') || isempty(doerrorbar)
        doerrorbar = 1;
    end
    u_ok = allu_ok & all(any(~isnan(psns),2),3);
    psnnonpref = psns(u_ok,:,1);
    psnpref = psns(u_ok,:,2);

    % plot nonpref
    psn = psnnonpref;
%     lcolor = 0.5*[1 1 1];
    plot_psn(x,psn,colors{2},gauss_smoothing_sigma,doerrorbar);
    hold on;
    
    %plot pref
    psn = psnpref;
%     lcolor = 0*[1 1 1];
    plot_psn(x,psn,colors{1},gauss_smoothing_sigma,doerrorbar);
    linex(0,'k'); liney(0,'k');
    xlim(xlims);

    set(gca,'FontSize',20);
    text(max(xlim)-0.05*range(xlim),max(ylim)-0.95*range(ylim),['n=' num2str(round(size(psn,1),2))],'FontSize',20,'HorizontalAlignment','right','VerticalAlignment','bottom');
    text(max(xlim)-0.95*range(xlim),max(ylim)-0.05*range(ylim),vstext,'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','top');

    xlabel(['Time from ' psname ' onset (ms)']);
    ylabel('Normalized activity');
end