%% intrinsicproperties.m
scriptname = 'intrinsicproperties';

%% table
basettl = 'intrinsic_properties_table';
ttl = [monktxt ', Offer ' num2str(offi) ' conds, ' basettl];
propnames = {'Mean firing rate (95% CI), baseline','Spike width','IR, baseline','CV, baseline'};
props = [baserate,SaveDur,basemedianirreg,basecv];
colorps = [relevant_ps(:,strcmp(xname,'R'),end),relevant_ps(:,strcmp(xname,'O'),end)];
T = {};
conds = struct(...
    'name',...
        {'All','Task-responsive','Not task-responsive'},...
    'thisalluok',...
        {this_allu_ok_manex & this_allu_ok_mintrials,...
        this_allu_ok_manex & this_allu_ok_mintrials & this_allu_ok_initin,...
        this_allu_ok_manex & this_allu_ok_mintrials & ~this_allu_ok_initin...
        });
for condi=1:3
    T{condi,1} = [conds(condi).name ' (n=' num2str(sum(conds(condi).thisalluok)) ')'];
    for propi = 1:length(propnames)
        thisprop = props(conds(condi).thisalluok,propi);
        T{condi,propi+1} = [num2str(round(mean(thisprop),2)) char(177) num2str(round(std(thisprop)/sqrt(length(thisprop)),2)) ' (' num2str(round(prctile(thisprop,2.5),2)) '-' num2str(round(prctile(thisprop,97.5),2)) ')'];
    end
end
T = cell2table(T);
T.Properties.VariableNames = ['Neurons',propnames];
ttl = [strjoin(propnames,'_') '_scatters'];
figname = 'supp_table_2';
savepath=fullfile(neuralsavedir, figname);
writetable(T,[savepath '.txt'],'Delimiter','\t');

%% main fig density, hist, scatter
basettl = ['main_scatter'];
thisalluok = this_allu_ok_manex & this_allu_ok_mintrials;
propnames = {'Mean firing rate, baseline','Spike width','IR, baseline','IR, all epochs','CV, baseline'};
props = [baserate,SaveDur,basemedianirreg,medianirreg,basecv]; props = props(thisalluok,:);
colorps = [relevant_ps(thisalluok,strcmp(xname,'R')),relevant_ps(thisalluok,strcmp(xname,'O'))];
fillp = ~this_allu_ok_initin(thisalluok);
colorps(~fillp<sigthresh,:) = NaN;
xi = 1; yi = 2;
x = props(:,xi); y = props(:,yi);

figure('Visible','on');
tiledlayout(7,1,'TileSpacing','none');
set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.5 0.9]);

nexttile(2,[6 1]);
myScatterWithHistogram(propnames{xi},propnames{yi},x,y,[],[],[],fillp,[],{[0 1.03*max(x)],[0 1.03*max(y)]},[],[],[],sigthresh,1,lgtxt); 
h = legend;
h.Visible = 'off';
thisxlim = xlim;

withist = 1;
if withist
    nexttile(1,[1 1]);
    h = histogram(x,'FaceAlpha',0.5,'FaceColor',0.5*[1 1 1],'BinWidth',1);
    hold on;
    plot(mean(x),0.9*max(ylim),'Marker','v','MarkerFaceColor',0.5*[1 1 1],'MarkerEdgeColor',0.5*[1 1 1],'MarkerSize',15);

    xlim(thisxlim);
    h=gca; h.YAxis.TickLength = [0 0];
    xticklabels([]);
end

set(findall(gcf,'type','axes'),'TickDir','out');
set(findall(gcf,'-property','FontSize'),'FontSize',20);
ttl = [monktxt ', Offer ' num2str(offi) ' conds, ' basettl];
ttl = [ttl ' ' regexprep([propnames{xi} 'vs' propnames{yi} ternary(withist,'withhist','')],' ','_')];
figname = 'supp_figure_3a';
savepath=fullfile(neuralsavedir, figname);
saveas(gcf,[savepath '.png']);
saveas(gcf,[savepath '.fig']);
close(gcf);

propnames = {'firing rate, baseline','spike width','IR (median), baseline'};
props = [baserate,SaveDur,basemedianirreg];
attrs = {'R'};
attrlabels = {'E[R]'};

figure('Visible','on'); set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.95 0.95]);
tiledlayout(length(attrs),length(propnames));
for ai=1
    colorxpypxy_ = cell(1,4);
    [colorxpypxy_{:}] = deal(ones(size(colorxpypxy{1})));
    colorxpypxy_{ai} = colorxpypxy{ai};
    for propi=1:length(propnames)
        nexttile();
        x = props(this_allu_ok,propi);
        y = relevant_bs(this_allu_ok,strcmp(xname,attrs{ai}));

        maxlim = {[median([min(x),max(x)])+1.05*range(x)*[-1/2 1/2]],[median([min(y),max(y)])+1.05*range(y)*[-1/2 1/2]]};

        x0 = 0; y0 = 0; markersize = 20; nohist = 1;

        myScatterWithHistogram(propnames{propi}, attrlabels{ai}, x, y, [], [], [], [], [], maxlim, x0, y0, markersize, sigthresh, nohist, lgtxt, colorxpypxy_,'none');
    end
end
ttl = regexprep([strjoin(propnames,'_') 'vs' strjoin(attrlabels,'_') num2str(sum(this_allu_ok))],' ','_');
figname = 'supp_figure_3b';
savepath=fullfile(neuralsavedir, figname);
saveas(gcf,[savepath '.png']);
saveas(gcf,[savepath '.fig']);

close(gcf);
