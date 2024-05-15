%% neural proportion bar plots
scriptname = 'neural_proportion_bar_plots';

basettl = ['proportions of neurons with Offer ' num2str(offi) ' coding by GLM'];
basettl = [monktxt ', ' basettl];
ttl = basettl;
figname = 'figure_3c';

selxnames = {'R','O','S'};
selxnamelabels = get_xnamelabels(selxnames);

[~,~,attridx]=intersect(selxnames,xname,'stable');
sig = relevant_ps(this_allu_ok,attridx,offi) < sigthresh;

% clopper pearson confidence intervals
if 1
    ci95 = NaN(2,length(selxnames));
    ci68 = ci95;
    for pdi=1:length(selxnames)
        pd = fitdist(double(sig(:,pdi)),'Binomial');
        tmp = paramci(pd,'Alpha',sigthresh);
        ci95(:,pdi) = tmp(:,2);
        tmp = paramci(pd,'Alpha',0.32);
        ci68(:,pdi) = tmp(:,2);
    end
    negerr = mean(sig)-min(ci68);
    poserr = max(ci68)-mean(sig);
end

p = myBinomTest(sum(sig),length(sig),sigthresh,'one');
RToff12stats(offi).p = p;
RToff12stats(offi).prop = mean(sig);
RToff12stats(offi).err = [negerr; poserr];
RToff12stats(offi).n = length(sig);

% start plotting
figure('Visible','on');
liney(sigthresh,'k--');
hold on;

errorbar(1:length(selxnames),mean(sig),negerr,poserr,"LineStyle","none","Marker","none","Color",0*[1 1 1]);
bar(mean(sig),'FaceColor',0.5*[1 1 1]);

xticks(1:length(attridx));
sigstar(xticks,p,[],sigthresh);

RO_chars = {'R','O'};
[~,lesser_of_RO] = min(sum(sig(:,ismember(selxnames,RO_chars))));
lesser_of_RO = RO_chars{lesser_of_RO};
tmp = find(ismember(selxnames,{lesser_of_RO,'S'}));

set(gca,'FontSize',20);
xticklabels(selxnamelabels);
ylabel('fraction of cells with significant coding');
text(max(xlim)-0.05*range(xlim),max(ylim)-0.05*range(ylim),['n=' num2str(sum(this_allu_ok))],'FontSize',20,'HorizontalAlignment','right','VerticalAlignment','top');

title(ttl);
savepath=fullfile(neuralsavedir, figname);
saveas(gcf,[savepath '.fig']);
saveas(gcf,[savepath '.png']);
close(gcf);
