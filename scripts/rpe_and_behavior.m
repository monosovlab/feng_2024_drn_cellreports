%% rpe_and_behavior.m
scriptname = 'rpe_and_behavior';

chr_names = {'R','O','S','Lft','Cnt','Rgt','2'};
chr_namelabels = get_xnamelabels(chr_names);
propensities = {NaN(befits{1,1}.n,length(chr_names)),NaN(befits{1,2}.n,length(chr_names))};

off1value = rpes;

nuisanceregs = {'Lft','Rgt','2'};

for mi_=1:2
    propensities{mi_} = getPropensities(befits{1,mi_},chr_names); 
    off1value{mi_} = getValues(befits{1,mi_},1,{'Lft','Rgt','2'});
end

tclttl = 'Prev RPE effect on behavior';
onlydo255025 = 1; dosurprise = 0; plus_minus_only = 1;
yplotlim = 1.3*[-1000 1000];
figure; set(gcf,'Position',get(0,'ScreenSize').*[1 1 1 0.5]);
tlsize = [length(monkeyabbr),length(chr_names)+1]; 
corrstats = NaN([2 2 tlsize]);
tcl = tiledlayout(tlsize(1),tlsize(2)); kwstats = NaN([2 tlsize]);
axs = {};
for mi_=1:2
    for chri=1:tlsize(2)
        axs{mi_,chri} = nexttile;
        switch chri
            case length(chr_names)+1
                xname_ = ['Fix RT'];
                ttl = xname_;
                x = fixRTs{mi_};
                if mi_<2
                    xl = [0.9 1.1];
                else
                    xl = [0.75 0.95];
                end
            otherwise
                xname_ = ['% stay with prev trial ' chr_namelabels{chri} ' pref'];
                ttl = [chr_namelabels{chri}];
                tmp_cur = sign(propensities{mi_}(:,chri)); tmp_prev = [0; tmp_cur(1:end-1)]; 
                x = tmp_prev.*tmp_cur;
                xl = [0.3 0.7];
        end
         yname_ = 'prev. RPE';
         y = [0; rpes{mi_}(1:end-1)]; % trial i's behavioral change (relative to trial i-1), and trial i-1 rpe
        
         if dosurprise
            y = abs(y);
            theseedges = [-inf 20 350 inf]; 
            yname_ = regexprep(yname_,'RPE','Surprise');

        else
            theseedges = [-inf -350 -280 20 350 inf]; 
            if plus_minus_only
                theseedges([2 end-1]) = [];
            end
        end

        removethese = x==0;
        if onlydo255025
            tmp = squeeze(altbefits{mi_,end}.probs(:,2,:));
            tmptmp = sub2ind(size(tmp),[1:size(tmp,1)]',altbefits{mi_,end}.yreg+1);
            tmp = tmp(tmptmp);
            removethese = removethese | (2*[0; tmp(1:end-1)]~=1);
        end
        x(removethese) = []; y(removethese) = []; 
        x = (x+1)/2;

        binmedians = NaN(length(theseedges)-1,1);
        x_mean = binmedians;
        x_CI = NaN(length(theseedges)-1,2);
        for bi=1:length(binmedians)
            inbin = inbounds(y,theseedges([0:1]+bi));
            binmedians(bi) = median(y(inbin));
            y(inbin) = binmedians(bi);
            x_mean(bi) = mean(x(inbin));
            if chri<(length(chr_names)+1)
                pd = fitdist(x(inbin),'Binomial');
                tmp = paramci(pd,'Alpha',0.32);
                ci68 = tmp(:,2);
                x_CI(bi,1) = x_mean(bi)-min(ci68);
                x_CI(bi,2) = max(ci68)-x_mean(bi);
            else
                x_CI(bi,:) = [1 1]*(std(x(inbin)/sqrt(sum(inbin))));
            end
        end
        
        [c,p] = corr(x,y,'type','Spearman');
        corrstats(1,:,mi_,chri) = [c p];
        [c,p] = corr(x,y,'type','Pearson');
        corrstats(2,:,mi_,chri) = [c p];
        
        [ksp,kstbl] = kruskalwallis(x,y,'off');
        kwstats(1,mi_,chri) = kstbl{2,5};
        kwstats(2,mi_,chri) = ksp;
        
        yl = 1.2*max(abs(binmedians))*[-1 1];
        plot(x_mean,binmedians,'Marker','o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',2);
        hold on;
        errorbar(x_mean,binmedians,0,0,x_CI(:,1),x_CI(:,2),'HandleVisibility','off','color','k','CapSize',0,'LineWidth',2);
        
        if chri<(length(chr_names)+1)
            linex(0.5,'k');
        end
        liney(0,'k');
        xlim(xl); ylim(yl);
        
        xlabel(xname_);
        ylabel(yname_);
        if chri<(length(chr_names)+1)
            xticklabels(arrayfun(@(x)num2str(x*100),xticks,'UniformOutput',false));
        end
    end
    if chri==1
        ylabel({monkeyabbr(mi_),get(gca,'YLabel').String},'FontSize',20);
    else
        ylabel([]);
    end
    
    title(ttl);
    axis square;
end


allp = cat(3,squeeze(corrstats(1,2,:,:)),squeeze(kwstats(2,:,:)));
[~, ~, ~, adj_p] = fdr_bh(allp);
for mi_=1:2
    for chri=1:tlsize(2)
        axes(axs{mi_,chri});
        text(max(xlim)-0.95*range(xlim),max(ylim)-0.6*range(ylim),['KW \chi^2 = ', char(string(round(kwstats(1,mi_,chri),2)))],'FontSize',10,'HorizontalAlignment','left','VerticalAlignment','bottom');
        text(max(xlim)-0.95*range(xlim),max(ylim)-0.7*range(ylim),['p_{adj} = ', char(string(round(adj_p(mi_,chri,2),2)))],'FontSize',10,'HorizontalAlignment','left','VerticalAlignment','bottom');
        text(max(xlim)-0.95*range(xlim),max(ylim)-0.85*range(ylim),['rho = ', char(string(round(corrstats(1,1,mi_,chri),2)))],'FontSize',10,'HorizontalAlignment','left','VerticalAlignment','bottom');
        text(max(xlim)-0.95*range(xlim),max(ylim)-0.95*range(ylim),['p_{adj} = ', char(string(round(adj_p(mi_,chri,1),3)))],'FontSize',10,'HorizontalAlignment','left','VerticalAlignment','bottom');
    end
end

if onlydo255025
    tclttl = [tclttl ', only 25/50/25'];
end
if dosurprise
    tclttl = regexprep(tclttl,'RPE','Surprise');
end
if plus_minus_only
    tclttl = [tclttl ', plus_minus_only'];
end
title(tcl,tclttl);

figname = 'supp_figure_7b';

savepath=fullfile(sigthreshsavedir,figname);
saveas(gcf,[savepath '.fig']);
saveas(gcf,[savepath '.png']);
close(gcf);

nqs=5;
mis_=1:2;
figure; set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.9 0.6]);
tcl=tiledlayout(length(mis_),1);
for mi_=mis_
    RTs_ = off1RTs;
    values_ = off1value;
    xname_ = 'SV1';
    yname_ = 'Offer 1 init fix reciprocal RT';

    RTs_{1}(RTs_{1}==0) = NaN;
    RTs_{2}(RTs_{2}==0) = NaN;

    y = [];
    x = [];
    chose2 = true([]);
    thisisinfo = [];
    if mi_
        mis__=mi_;
        monktxt = monkeyabbr(mi_);
    else
        mis__=1:2;
        monktxt = 'Both monks';
    end
    for mi__=mis__
        chose2 = vertcat(chose2,logical(befits{1,mi__}.y));
        thisisinfo = vertcat(thisisinfo,isinfo{mi__});

        y = vertcat(y,1./RTs_{mi__});
        stattxt = 'raw';

        x = vertcat(x,values_{mi__});
    end

    xname_ = xname_;
    maxlim = {[-1 1]*max(abs(x)), [0 1]*max(abs(y))};

    nexttile;

    tmpsel = true(size(chose2));

    if nqs>2
        tmp = nqs-1;
    else
        tmp = 0.5;
    end
    tmpq = quantile(x(tmpsel),tmp);

    tmpq = [-inf tmpq inf]; 
    xmedians = NaN(size(x)); uxmedians = NaN(nqs,1); tmpmeans = uxmedians; negerrs = uxmedians; poserrs = uxmedians;
    rsps = NaN(nqs,1); rspgroups = {};
    for qi=1:nqs
        hold on;
        qii = mod(qi,nqs)+1;
        uxmedians(qi) = median(x(inbounds(x,tmpq(qi+(0:1)))));
        uxmedians(qii) = median(x(inbounds(x,tmpq(qii+(0:1)))));
        xmedians(inbounds(x,tmpq(qi+(0:1)))) = uxmedians(qi);
        xmedians(inbounds(x,tmpq(qii+(0:1)))) = uxmedians(qii);
        
        if nqs>2 || qi==1
            rsps(qi) = ranksum(y(xmedians==uxmedians(qi) & tmpsel),y(xmedians==uxmedians(qii) & tmpsel));
            rspgroups{qi} = [uxmedians(qi) uxmedians(qii)];
        end

        tmpy = y(xmedians==uxmedians(qi) & tmpsel);
        tmpmean = nanmean(tmpy);
        tmperr = nanstd(tmpy)./sqrt(sum(~isnan(tmpy))); negerrs(qi) = tmperr; poserrs(qi) = tmperr;
        tmpmeans(qi) = tmpmean;
    end

    plot(uxmedians,tmpmeans,"Marker",'o','MarkerFaceColor',0*[1 1 1],'MarkerEdgeColor','none','MarkerSize',8,'Color',0*[1 1 1],'LineWidth',2);
    errorbar(uxmedians,tmpmeans,negerrs,poserrs,'HandleVisibility','off','color','k','CapSize',0,'LineWidth',2,'LineStyle','none');

    tmp = all(~isnan([x y]),2);
    [c,p] = corr(x(tmp),y(tmp),'type','Spearman');
    text(0.1*range(xlim)+min(xlim),0.9*range(ylim)+min(ylim),['rho = ' num2str(c)]);
    text(0.1*range(xlim)+min(xlim),0.86*range(ylim)+min(ylim),['p = ', num2str(p)]);
    
    [c,p] = corr(x(tmp),y(tmp),'type','Pearson');
    text(0.1*range(xlim)+min(xlim),0.75*range(ylim)+min(ylim),['r = ' num2str(c)]);
    text(0.1*range(xlim)+min(xlim),0.71*range(ylim)+min(ylim),['p = ', num2str(p)]);
    sigstar(rspgroups(1:(nqs-1)),rsps(1:(nqs-1)),[],[],[],[],1);
    if (nqs>=3)
        sigstar(rspgroups(end),rsps(end),[],[],[],[],1);
    end

    xlabel([xname_ ' bin median']);
    ylabel(yname_);
    axis square;

    ylabel({monktxt,get(gca,'YLabel').String},'FontSize',10);

    text(max(xlim)-0.05*range(xlim),max(ylim)-0.05*range(ylim),['n = ', num2str(sum(~isnan(y(tmpsel))))],'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','bottom');
end

figname = 'supp_figure_2c';

savepath=fullfile(sigthreshsavedir,figname);
saveas(gcf,[savepath '.fig']);
saveas(gcf,[savepath '.png']);
close(gcf);

%% value vs attr model comparison off1 RT

monknames = {}; xnames = {}; notnuisance = {};
attrb = {}; attrstats = {};

for mi_=1:2
    monknames{mi_} = char(monkeyabbr(mi_));
    y = off1RTs{mi_}; isvalid = y>0;
    [~,attrx,~,xnames{mi_},notnuisance{mi_}] = getValues(befits{1,mi_},1,{'Lft','Rgt','2'});
    valx = off1value{mi_};
    
    y = 1./y(isvalid);
    attrx = attrx(isvalid,:);

    theseinds = 1:length(y);
    [attrb_,~,attrstats_] = glmfit(attrx(theseinds,:),y(theseinds),"normal",'link','identity','constant','on');
    attrb{mi_} = attrb_;
    attrstats{mi_} = attrstats_;
end

%% RT vs choice beta scatters

figure; set(gcf,'Position',get(0,'ScreenSize')); 
xname_ = '\beta_{Offer\_1\_reciprocal\_latency}';
yname_ = '\beta_{preference}';

mis__ = 1:2;
monktxt = 'both monkeys';

x = []; y = [];
xerr = []; yerr = [];
xp = []; yp = [];
for mi__=mis__
    x = vertcat(x,attrb{mi__}(2:end));
    xerr = vertcat(xerr,attrstats{mi__}.se(2:end));
    xp = vertcat(xp,attrstats{mi__}.p(2:end));
    y = vertcat(y,befits{1,mi__}.b(notnuisance{mi__}));
    yerr = vertcat(yerr,befits{1,mi__}.stats.se(notnuisance{mi__}));
    yp = vertcat(yp,befits{1,mi__}.stats.p(notnuisance{mi__}));
end

myScatterWithHistogram(yname_,xname_,y,x,yp,xp,{yerr,xerr},[],[],'seplims',0,0,[],sigthresh,1);
axis square;
title(monktxt);

ttl = 'off1RTandchoiceweightsscatter';
figname = 'supp_figure_2d';

savepath=fullfile(sigthreshsavedir,figname);
saveas(gcf,[savepath '.fig']);
saveas(gcf,[savepath '.png']);
close(gcf);

%% helper function
function propensities = getPropensities(befit,chr_names)
    xreg = {};
    for cni=1:length(chr_names)
        if strcmp(chr_names{cni},'Cnt')
            tmp1 = befit.xreg{strcmp(befit.xname,'Lft')};
            tmp2 = befit.xreg{strcmp(befit.xname,'Rgt')};
            xreg{cni} = tmp1;
            xreg{cni}.terms{1} = ~(tmp1.terms{1} | tmp2.terms{1});
        else
            xreg{cni} = befit.xreg{strcmp(befit.xname,chr_names{cni})};
        end
    end
    yreg = befit.yreg;
    propensities = (2*yreg-1).*eglm_make_regression_inputs(xreg,yreg);
end

function [values,x,b,xname,notnuisance] = getValues(befit,valuemode,nuisanceregs) % nuisance regs put at top?
    notnuisance = ~ismember(befit.xname,nuisanceregs);
    b = befit.b(notnuisance);
    switch valuemode
        case 0
            x = befit.x(:,notnuisance);
        case {1,2}
            x = eglm_make_regression_inputs(befit.xreg,befit.yreg,true(size(befit.yreg)),1);
            x = squeeze(x(:,notnuisance,valuemode));
    end
    values = x*b;
    xname = befit.xname(notnuisance);
end

function lldiff = getLLdiff(y,attrx,valx,theseinds)
    [attrb,~,attrstats] = glmfit(attrx(theseinds,:),y(theseinds),"normal",'link','identity','constant','on');
    attrypred = glmval(attrb,attrx(theseinds,:),'identity','constant','on');
    attrlltr = log(normpdf(y(theseinds),attrypred,attrstats.sfit));
    attrll = sum(attrlltr);

    [valb,~,valstats] = glmfit(valx(theseinds,:),y(theseinds),"normal",'link','identity','constant','on');
    valypred = glmval(valb,valx(theseinds,:),'identity','constant','on');
    vallltr = log(normpdf(y(theseinds),valypred,valstats.sfit));
    valll = sum(vallltr);
    
    lldiff = attrll - valll;
end

