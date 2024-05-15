function [xp,yp] = plot_offer_scatter(befits,xs,xps,xses,ys,yps,yses,xname,xnames,xplotnames,allu_ok,monkeyabbr,ii,jj,doflip,modestr,sigthresh,xylims,nohist,basecolors)
    if ~exist('nohist','var') || isempty(nohist)
        nohist=0;
    end 
    if strcmp(modestr,'ROC')
        x = roc_auc_and_p(allu_ok,1,ii); xp = roc_auc_and_p(allu_ok,2,ii);
        y = roc_auc_and_p(allu_ok,1,jj); yp = roc_auc_and_p(allu_ok,2,jj);
        xylims = [0,1];
    elseif strcmp(modestr,'GLM')
        x = xs(allu_ok,strcmp(xname,xnames{ii})); xp = xps(allu_ok,strcmp(xname,xnames{ii}));
        y = ys(allu_ok,strcmp(xname,xnames{jj})); yp = yps(allu_ok,strcmp(xname,xnames{jj}));
        if exist("xses",'var') && ~isempty(xses)
            xse = 1.96*xses(allu_ok,strcmp(xname,xnames{ii}));
        else
            xse = [];
        end
        if exist("yses",'var') && ~isempty(yses)
            yse = 1.96*yses(allu_ok,strcmp(xname,xnames{jj}));
        else
            yse = [];
        end
        if ~exist('xylims','var')
            xylims = 1.1*max(abs([x;y]))*[-1 1];
        end
    else
        error('modestr');
    end
    if ~exist('basecolors','var') || isempty(basecolors)
        basecolors = {[1 0 0],[0 0 1]};
    end
    if doflip
        for xyi=1:2
            ubs = [];
            switch xyi
                case 1
                    ij = ii;
                case 2
                    ij = jj;
            end
            for gi=1:size(befits,2)
                if any(strcmp(befits{1,gi}.xname,xnames{ij}))
                    tmp = befits{1,gi}.b(strcmp(befits{1,gi}.xname,xnames{ij}));
                else
                    tmp = NaN;
                end
                ubs = [ubs; tmp*ones(ternary(gi==1,224,length(allu_ok)-224),1)];
            end
            switch xyi
                case 1
                    x = sign(ubs(allu_ok)).*x;
                case 2
                    y = sign(ubs(allu_ok)).*y;
            end
        end
    end
    
    if ii==jj
        warning('implicit hardcoding off1=x vs off2=y');
        lg = {['neither (p > ' num2str(sigthresh) ')'],'Offer 1','Offer 2','Both Offers'};
    else
        lg = {['neither (p > ' num2str(sigthresh) ')'],xplotnames{ii},xplotnames{jj},[xplotnames{ii} ' and ' xplotnames{jj}]};
    end

    myScatterWithHistogram(lg{2},lg{3},x,y,xp,yp,{xse,yse},[],[],xylims,[],[],[],sigthresh,nohist,[],[],[],[],basecolors);
end