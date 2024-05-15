function [h, stats] = myScatterWithHistogram(xname,yname,x,y,xp,yp,xyerr,fillp,colorconds,maxlim,x0,y0,markersize,sigthresh,nohist,lgtxt,colorxpypxy,lglocation,usecurfig,basecolors,fontsize,dojitter,fillname)
    if ~exist('xp','var') || isempty(xp)
        xp = NaN(size(x));
    end
    if ~exist('yp','var') || isempty(yp)
        yp = NaN(size(x));
    end
    removethese = any(isnan([x y]),2);
    x(removethese) = []; y(removethese) = []; 
    xp(removethese) = []; yp(removethese) = [];
    if ~exist('colorxpypxy','var') || isempty(colorxpypxy) % doesn't affect hists
        colorxpypxy = {xp, yp, x, y};
    else
        for ii=1:length(colorxpypxy)
            if ~isempty(colorxpypxy{ii})
                colorxpypxy{ii} =  colorxpypxy{ii}(~removethese); 
            end
        end
    end
    if any(removethese)
        warning('removing nans (%d/%d data points)',sum(removethese(:)),numel(removethese));
    end
    if ~exist('sigthresh','var') || isempty(sigthresh)
        sigthresh = 0.05;
    end
    if ~exist('fillp','var') || isempty(fillp)
        fill_ok = true(size(x));
    else
        fill_ok = fillp < sigthresh;
    end
    if ~exist('maxlim','var') || isempty(maxlim)
        maxlim = max(max(abs(x)),max(abs(y)))*[-1.05 1.05];
    end
    if ischar(maxlim)
        assert(strcmp(maxlim,'seplims'),'should be seplims');
        maxlim = {};
        maxlim{1} = max(abs(x))*[-1.05 1.05];
        maxlim{2} = max(abs(y))*[-1.05 1.05];
    end
    if iscell(maxlim)
        xlims = maxlim{1};
        ylims = maxlim{2};
    else
        xlims = maxlim;
        ylims = maxlim;
    end
    if exist('xyerr','var') && ~isempty(xyerr) && ~any(cellfun(@(x)isempty(x),xyerr))
        tmpx = max(xyerr{:,1});
        tmpy = max(xyerr{:,2});
        xlims = xlims + tmpx*[-1 1];
        ylims = ylims + tmpy*[-1 1];
    end
    if ~exist('basecolors','var') || isempty(basecolors)
        basecolors = {[1 0 0],[0 0 1]};
    end
    
    if ~exist('colorconds','var') || isempty(colorconds)
        tmp = sum(vertcat(basecolors{:}));
        colorconds = struct(...
            'is_ok',{...
                @(xp,yp,x,y) ~(xp < sigthresh | yp < sigthresh),...
                @(xp,yp,x,y) xp < sigthresh & ~(yp < sigthresh),...
                @(xp,yp,x,y) ~(xp < sigthresh) & yp < sigthresh,...
                @(xp,yp,x,y) xp < sigthresh & yp < sigthresh...
            },...
            'color',{[0 0 0],basecolors{1},basecolors{2},tmp/max(tmp)}...
        );
    end
    if ~exist('xyerr','var') 
        xyerr = [];
    end
    if ~exist('x0','var') || isempty(x0)
        x0 = 0;
    end
    if ~exist('y0','var') || isempty(y0)
        y0 = 0;
    end
    if ~exist('markersize','var') || isempty(markersize)
        markersize=20;
    end
    if ~exist('nohist','var') || isempty(nohist)
        nohist=0;
    end 
    if ~exist('lgtxt','var') || isempty(lgtxt)
        lgtxt = {['neither (p>' num2str(sigthresh) ')'],xname,yname,'both'};
    end
    if ~exist('lglocation','var') || isempty(lglocation)
        lglocation = 'southeast';
    end
    if ~exist('usecurfig','var') || isempty(usecurfig)
        usecurfig = true;
    end
    if ~exist('fillname','var') || isempty(fillname)
        fillname = 'filled';
    end
    if ~nohist
        if ~usecurfig
            h = figure;
        end
        h = tiledlayout(5,5, 'Padding', 'none', 'TileSpacing', 'none');
        nexttile(6,[4 4]);
    %     subplot(5,5,[1*5+[1:4] 2*5+[1:4] 3*5+[1:4] 4*5+[1:4]]);
    else
        h = [];
    end
    if ~exist('fontsize','var') || isempty(fontsize)
        fontsize = markersize*0.5;
    end
    if ~exist('dojitter','var') || isempty(dojitter)
        dojitter = false;
    end
    uselsline = 0;
    if dojitter
        jitterxy = normrnd(0,0.01,[length(x),2]).*[range(xlims) range(ylims)];
    else
        jitterxy = [0 0];
    end
    jitteredx = x + jitterxy(:,1);
    jitteredy = y + jitterxy(:,2);

    stats = struct('corr',NaN(2,2));

    if numel(x)>1
        scatter(x,y,'.','MarkerEdgeColor','none','HandleVisibility','off');
        xlim(xlims);
        ylim(ylims);
        if uselsline
            l = lsline();
        else
            [tmpm,tmpb] = lsqfitma(x,y);
            l = line(xlims,tmpm*xlims+tmpb);
        end
        l.Color = 0.7*[1 1 1];
        l.LineWidth = 2;
        l.HandleVisibility = 'off';
    end
    hold on;
    
    if ~isempty(xyerr)
        xpos = xyerr{1,1};
        ypos = xyerr{1,2};
        if size(xyerr,1)>1
            xneg = xyerr{2,1};
            yneg = xyerr{2,2};
        else
            xneg = xpos;
            yneg = ypos;
        end
        errorbar(jitteredx,jitteredy,yneg,ypos,xneg,xpos,'Color',0.7*[1 1 1],'LineStyle','none','HandleVisibility','off','CapSize',0);
    end
    if numel(x)>1
        uistack(l,'top');
    end
    
    for cci=1:length(colorconds)
        if isstruct(colorconds)
            is_ok = colorconds(cci).is_ok(colorxpypxy{1},colorxpypxy{2},colorxpypxy{3},colorxpypxy{4});
            thiscolor = colorconds(cci).color;
        elseif isnan(colorconds)
            is_ok = true(size(fill_ok));
            thiscolor = [0 0 0];
        else
            error('expected colorconds to be struct or nan');
        end
        if mean(fill_ok)<1
            scatter(jitteredx(is_ok & ~fill_ok),jitteredy(is_ok & ~fill_ok),markersize,thiscolor,'MarkerFaceColor',[1 1 1],'LineWidth',0.07*markersize,'HandleVisibility','off');
        end           
        scatter(jitteredx(is_ok & fill_ok),jitteredy(is_ok & fill_ok),markersize,thiscolor,'MarkerFaceColor','flat','LineWidth',0.07*markersize,'DisplayName',lgtxt{cci});
    end
    
    axis square;
    
    xline(x0,'k','HandleVisibility','off');
    yline(y0,'k','HandleVisibility','off');
    
    % need to do it again?
    xlim(xlims);
    ylim(ylims);

    if numel(x)>1
        text(max(xlim)-0.95*range(xlim),max(ylim)-0.95*range(ylim),['n=' num2str(numel(x))],'FontSize',fontsize,'HorizontalAlignment','left','VerticalAlignment','bottom');
        [c,p] = corr(x,y,'type','Spearman');
        text(0.1*range(xlim)+min(xlim),0.9*range(ylim)+min(ylim),['rho = ' num2str(c)],'FontSize',fontsize);
        text(0.1*range(xlim)+min(xlim),0.86*range(ylim)+min(ylim),['p = ', num2str(p)],'FontSize',fontsize);
        stats.corr(1,:) = [c p];
        [c,p] = corr(x,y,'type','Pearson');
        text(0.1*range(xlim)+min(xlim),0.75*range(ylim)+min(ylim),['r = ' num2str(c)],'FontSize',fontsize);
        text(0.1*range(xlim)+min(xlim),0.71*range(ylim)+min(ylim),['p = ', num2str(p)],'FontSize',fontsize);
        stats.corr(2,:) = [c p];
    end
    
    xlabel(xname,'FontSize',fontsize,'FontWeight','bold');
    ylabel(yname,'FontSize',fontsize,'FontWeight','bold');
    if isstruct(colorconds)
        if mean(fill_ok)<1
            scatter(NaN,NaN,markersize,0.5*[1 1 1],'MarkerFaceColor',0.5*[1 1 1],'LineWidth',0.07*markersize,'DisplayName',fillname);
            scatter(NaN,NaN,markersize,0.5*[1 1 1],'MarkerFaceColor',[1 1 1],'LineWidth',0.07*markersize,'DisplayName',['not ' fillname]);
%             plot(max(xlims),max(ylims),'Marker','o','MarkerSize',markersize,'MarkerEdgeColor',0.5*[1 1 1],'MarkerFaceColor',0.5*[1 1 1],'LineWidth',0.07*markersize,'DisplayName');
%             plot(max(xlims),max(ylims),'Marker','o','MarkerSize',markersize,'MarkerEdgeColor',0.5*[1 1 1],'MarkerFaceColor',[1 1 1],'LineWidth',0.07*markersize);
%             lgtxt = [lgtxt,{fillname},{['not ' fillname]}];
        end
        if ~strcmp(lglocation,'none')
            legend('Location',lglocation);
        end
    end
    
    hold off;
    
    if ~nohist
        scatpos = get(gca,'Position');
        %% need to implement histogram
        histylim = 0;
        hisths = {};
        for tmpi=1:2
            switch tmpi
                case 1
    %                 hisths{tmpi} = subplot(5,5,1:4);
                    hisths{tmpi} = nexttile(1,[1 4]);
                    xy = x; xyp = xp;
                    thislim = xlims;
                case 2
    %                 hisths{tmpi} = subplot(5,5,(10:5:25));
                    hisths{tmpi} = nexttile(10,[4 1]);
                    xy = y; xyp = yp;
                    thislim = ylims;
                    hbw = range(ylims)/nbin;
            end
            nbin = 20;
            hbw = range(thislim)/nbin;
            histogram(xy,'FaceColor',0.7*[1 1 1],'BinWidth',hbw);
            hold on;
            histogram(xy(xyp < sigthresh),'FaceColor',basecolors{tmpi},'BinWidth',hbw);

            xlim(thislim);
            histylim = max([histylim ylim]);
        end
        for tmpi=1:2
            negi = 2*tmpi-3;
            switch tmpi
                case 1
    %                 hisths{tmpi} = subplot(5,5,1:4);
                    hisths{tmpi} = nexttile(1,[1 4]);
                    xy = x; xyp = xp; xy0 = x0;
                case 2
    %                 hisths{tmpi} = subplot(5,5,(10:5:25));
                    hisths{tmpi} = nexttile(10,[4 1]);
                    xy = y; xyp = yp; xy0 = y0;
            end
            ylim([0 histylim]);
            for tmpj=1:2
                neg = 2*tmpj-3;
                text([negi*0.05+0.5+neg*0.43]*range(xlim)+min(xlim),0.85*range(ylim)+min(ylim),[char(string(round(100*mean(neg*(xy-xy0) > 0 & xyp<sigthresh),1))),'%'],'FontSize',fontsize,'Rotation',(tmpi-1)*[-90]);
                text([negi*0.05+0.5+neg*0.43]*range(xlim)+min(xlim),0.7*range(ylim)+min(ylim),['p = ', char(string(round(myBinomTest(sum(neg*(xy-xy0) > 0 & xyp<sigthresh),length(xy),sigthresh/2),4)))],'FontSize',fontsize,'Rotation',(tmpi-1)*[-90]);
            end

            all_mean = mean(xy);
            all_median_p = signrank(xy);
            l = line([all_mean all_mean],ylim,'LineStyle',':','Color',0.7*[1 1 1],'LineWidth',3); uistack(l,'top');
            text(all_mean+0.01*diff(xlim),max(ylim)-0.05*diff(ylim),['p = ' num2str(all_median_p)],'FontSize',fontsize,'Rotation',(tmpi-1)*[-90],'VerticalAlignment','top','HorizontalAlignment',ternary(tmpi==1,'left','right'));
            
            sig_mean = nan; sig_median_p = sig_mean;
            if any(xyp < sigthresh)
                sig_mean = mean(xy(xyp < sigthresh));
                sig_median_p = signrank(xy(xyp < sigthresh));
            end
            l = line([sig_mean sig_mean],ylim,'LineStyle',':','Color',basecolors{tmpi},'LineWidth',3); uistack(l,'bottom');
            text(sig_mean+0.01*diff(xlim),max(ylim)-0.1*diff(ylim),['p = ' num2str(sig_median_p)],'FontSize',fontsize,'Rotation',(tmpi-1)*[-90],'VerticalAlignment','top','HorizontalAlignment',ternary(tmpi==1,'left','right'));

            l = line([0 0],ylim,'LineStyle','-','Color',[0 0 0],'LineWidth',1/2); uistack(l,'bottom');

            xticks([]);
            ylabel('#','FontSize',fontsize,'FontWeight','bold');
            tmp = get(gca,'Position'); tmp(tmpi+[0 2]) =  scatpos(tmpi+[0 2]);
%             set(gca,'Position',tmp);
        end
        view([90 -90]); % just for vertical hist
        set(gcf,'Position',markersize*[0 0 25 25])
    end
    if ~exist('h','var')
        h = [];
    end
end