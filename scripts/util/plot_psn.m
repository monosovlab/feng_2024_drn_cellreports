function h = plot_psn(x,psn,lcolor,gauss_smoothing_sigma,doerrorbar)
    if ~exist('doerrorbar','var') || isempty(doerrorbar)
        doerrorbar = 1;
    end
    smoothpsn = efilternan(squeeze(psn),'gauss',gauss_smoothing_sigma);
    smooththensepsn = std(smoothpsn)/sqrt(size(smoothpsn,1));
    smooththenmeanpsn = mean(smoothpsn,1);
    h = plot(x,smooththenmeanpsn,'Color',lcolor,'LineWidth',3);
    hold on;
    if doerrorbar && size(psn,1) > 1
        shadedErrorBar(x,smooththenmeanpsn,smooththensepsn,'lineProps',{'color',lcolor});
    end
end