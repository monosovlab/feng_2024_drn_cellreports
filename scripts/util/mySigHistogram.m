function [h, stats] = mySigHistogram(xy,xyp,sigthresh,colors,thislim,nbin,fontsize,h,xyname,numname)
    if ~exist('xyp','var') || isempty(xyp)
        xyp = NaN(size(xy));
    end
    if ~exist('h','var') || isempty(h)
        h = gca;
    end
    if ~exist('colors','var') || isempty(colors)
        colors = [0.7; 0].*[1 1 1];
    end
    if ~exist('thislim','var') || isempty(thislim)
        thislim = max(abs(xy))*[-1 1];
    end
    if ~exist('nbin','var') || isempty(nbin)
        nbin = 10;
    end
    if ~exist('fontsize','var') || isempty(fontsize)
        fontsize = 10;
    end
    if ~exist('xyname','var') || isempty(xyname)
        xyname = 'x';
    end
    if ~exist('numname','var') || isempty(numname)
        numname = 'y';
    end
    plotHist(xy,xyp,sigthresh,colors,thislim,nbin,h);
    [h, stats] = labelHist(xy,xyp,sigthresh,colors,fontsize,h,xyname,numname);
end

function h = plotHist(xy,xyp,sigthresh,colors,thislim,nbin,h)
%     nbin = 20;
    hbw = range(thislim)/nbin;
    histogram(xy,'FaceColor',colors(1,:),'BinWidth',hbw);
    hold on;
    histogram(xy(xyp < sigthresh),'FaceColor',colors(2,:),'BinWidth',hbw);
    
    xlim(thislim);
end

function [h, stats] = labelHist(xy,xyp,sigthresh,colors,fontsize,h,xyname,numname)
    negi=0;
    for tmpj=1:2
        neg = 2*tmpj-3;
        text([negi*0.05+0.5+neg*0.43]*range(xlim)+min(xlim),0.95*range(ylim)+min(ylim),[char(string(round(100*mean(neg*xy > 0 & xyp<sigthresh),1))),'%'],'FontSize',fontsize,'VerticalAlignment','top','HorizontalAlignment',ternary(tmpj==1,'left','right'));
        text([negi*0.05+0.5+neg*0.43]*range(xlim)+min(xlim),0.9*range(ylim)+min(ylim),['p = ', char(string(round(myBinomTest(sum(neg*xy > 0 & xyp<sigthresh),length(xy),sigthresh/2),4)))],'FontSize',fontsize,'VerticalAlignment','top','HorizontalAlignment',ternary(tmpj==1,'left','right'));
    end

    all_mean = mean(xy);
    all_median_p = signrank(xy);
    l = line([all_mean all_mean],ylim,'LineStyle',':','Color',colors(1,:),'LineWidth',3); uistack(l,'top');
    text(all_mean+0.01*diff(xlim),max(ylim)-0.05*diff(ylim),['p = ' num2str(all_median_p)],'FontSize',fontsize,'VerticalAlignment','top','HorizontalAlignment','left');
    
    sig_mean = nan; sig_median_p = sig_mean;
    if any(xyp < sigthresh)
        sig_mean = mean(xy(xyp < sigthresh));
        sig_median_p = signrank(xy(xyp < sigthresh));
    end
    l = line([sig_mean sig_mean],ylim,'LineStyle',':','Color',colors(2,:),'LineWidth',3); uistack(l,'bottom');
    text(sig_mean+0.01*diff(xlim),max(ylim)-0.1*diff(ylim),['p = ' num2str(sig_median_p)],'FontSize',fontsize,'VerticalAlignment','top','HorizontalAlignment','left');

    l = line([0 0],ylim,'LineStyle','-','Color',[0 0 0],'LineWidth',1/2); uistack(l,'bottom');

    text([negi*0.05+0.5-0.43]*range(xlim)+min(xlim),0.05*range(ylim)+min(ylim),['(n = ' num2str(length(xy)) ')'],'FontSize',fontsize,'VerticalAlignment','bottom','HorizontalAlignment','left');
%     xticks([]);
    ylabel(['# of ' numname],'FontSize',fontsize,'FontWeight','bold');
    xlabel(xyname,'FontSize',fontsize,'FontWeight','bold');

    stats.all_mean = all_mean;
    stats.all_median_p = all_median_p;

    stats.sig_mean = sig_mean;
    stats.sig_median_p = sig_median_p;
end