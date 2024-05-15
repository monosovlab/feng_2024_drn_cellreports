%% DistanceFromAC.m
scriptname = 'DistanceFromAC';

colorps = [relevant_ps(:,strcmp(xname,'R'),end),relevant_ps(:,strcmp(xname,'O'),end)];

%% table
tmp = combinedAtlasAligned(this_allu_ok_manex & this_allu_ok_mintrials,:);
T = [mean(tmp); std(tmp)];
T = T(:,2:-1:1);
T = array2table(T);
T.Properties.RowNames = {'Mean','SD'};
T.Properties.VariableNames = {'AP','LM'};
ttl = 'AP_ML';
figname = 'supp_table_3';
savepath=fullfile(neuralsavedir,figname);
writetable(T,[savepath '.txt'],'Delimiter','\t');

%% GLM vs ap/ml scatters
ai=1;
colorxpypxy_ = cell(1,4);
[colorxpypxy_{:}] = deal(ones(size(colorxpypxy{1})));
colorxpypxy_{ai} = colorxpypxy{ai};
x_ = combinedAtlasAligned(this_allu_ok,:);
y = colorbs(this_allu_ok,ai); 
yp = colorps(this_allu_ok,ai); 
xp = NaN(size(y));
figure; set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.9 0.9]); tiledlayout(1,2);
nexttile;
myScatterWithHistogram('AP','E[R]',x_(:,2),y,xp,yp,[],[],[],{[2 4],0.05*range(y)*[-1 1]+[min(y) max(y)]},0,0,[],sigthresh,1,[],colorxpypxy_,'none');
nexttile;
myScatterWithHistogram('LM','E[R]',x_(:,1),y,xp,yp,[],[],[],{[-2 2],0.05*range(y)*[-1 1]+[min(y) max(y)]},0,0,[],sigthresh,1,[],colorxpypxy_,'none');

ttl = ['ERvsAP_LM'];
figname = 'supp_figure_3c';
savepath=fullfile(neuralsavedir,figname);
saveas(gcf,[savepath '.fig']);
saveas(gcf,[savepath '.png']);
close(gcf);
