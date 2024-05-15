function h = myheatmap(cc,cp,clim,sigtextfontsize,ngrad,firststarp,xl,yl)
    if ~exist('clim','var') || isempty(clim)
            clim = [floor(min(cc(~isinf(cc)),[],'all')),ceil(max(cc(~isinf(cc)),[],'all'))];
    end
    if ~exist('sigtextfontsize','var') || isempty(sigtextfontsize)
            sigtextfontsize = 10;
    end
    if ~exist('ngrad','var') || isempty(ngrad)
            ngrad = 256;
    end
    if ~exist('firststarp','var') || isempty(firststarp)
            firststarp = 0.05;
    end
    if ~exist('cp','var') || isempty(cp)
        cp = NaN(size(cc));
    end
    blue_to_white = colorgradient([0 0 1],[1 1 1],ngrad);
    white_to_red = colorgradient([1 1 1], [1 0 0],ngrad);
    cm = [];
    if ~any(isnan(clim))
        if min(clim)<0
            cm = vertcat(cm,blue_to_white(1:end-1,:));
        end
        if inbounds(0,clim)
            cm = vertcat(cm,[1 1 1]);
        end
        if max(clim)>0
            cm = vertcat(cm,white_to_red(2:end,:));
        end
    end
    if 1
        imcc = cc;
        imcc(isinf(cc) & cc<0 & cp<firststarp) = clim(1); imcc(isinf(cc) & cc>0 & cp<firststarp) = clim(2);
        h = imagesc(imcc,'AlphaData',~isnan(imcc) & ~isinf(imcc));
    else
        h = imagesc(cc,'AlphaData',~isnan(cc));
    end
    if any(~isnan(cp),'all')
        if all(isnan(cp))
            for cci=1:size(cc,1)
                for ccj=1:size(cc,2)
                    text(1/(2*size(cc,2))+(ccj-1)/size(cc,2),1-(1/(2*size(cc,1))+(cci-1)/size(cc,1)),num2str(cc(cci,ccj)),'Units','Normalized','HorizontalAlignment','center','FontSize',sigtextfontsize)
                end
            end
        else
            assert(all(size(cc)==size(cp)));
            for cpi=1:size(cp,1)
                for cpj=1:size(cp,2)
                    if ~(isnan(cc(cpi,cpj)) || isinf(cc(cpi,cpj)))
                        %text(1/(2*size(cp,2))+(cpj-1)/size(cp,2),1-(1/(2*size(cp,1))+(cpi-1)/size(cp,1)),significance_text(cp(cpi,cpj),firststarp),'Units','Normalized','HorizontalAlignment','center','FontSize',sigtextfontsize)
                        text(cpj,cpi,sprintf('p=%.4f',cp(cpi,cpj)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',sigtextfontsize);
                    end
                end
            end
        end
    end
    colormap(gca,cm);
    colorbar;
    axis equal;
    axis tight;
    
    if ~any(isnan(clim))
        caxis(clim);
    end
    set(gca,'color',0.5*[1 1 1]);
    xticks(1:size(cc,2));
    yticks(1:size(cc,1));
    if exist('xl','var')
        xticklabels(xl);
    end
    if exist('yl','var')
        yticklabels(yl);
    end
end