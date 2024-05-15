%% psth
scriptname = 'selective_subpop_psths';

comparisonstatmode = 'diff';
dotile = 0;
if 1
    psth_bothoffi_relevant_bs = {};
    psth_bothoffi_relevant_ps = {};
    psth_bothoffi_relevant_ses = {};
    for halfi_=1:3
            [psth_bothoffi_relevant_bs{offi,halfi_}, psth_bothoffi_relevant_ps{offi,halfi_}, psth_bothoffi_relevant_ses{offi,halfi_}] = get_relevant_bs_ps(glm_analysis_periods,offi,halfi_,monkeyabbr,xname,allu_ok);
    end
end

flipmode = 1;
wi = 1;
valmodelis = [2 3];
attr2cutoffs = [0.03];
attrmodeli = 1;
valmode = 1;

attrefits = vertcat(glm_analysis_periods{1}(1).efits(:,wi,attrmodeli,3),glm_analysis_periods{2}(1).efits(:,wi,attrmodeli,3));
valefits = vertcat(glm_analysis_periods{1}(1).efits(:,wi,valmodelis(1)+(valmode-1)*3,3),glm_analysis_periods{2}(1).efits(:,wi,valmodelis(2)+(valmode-1)*3,3));
attr2 = cellfun(@(x)field_or_nil(x,'r2',NaN),attrefits);
valr2 = cellfun(@(x)field_or_nil(x,'r2',NaN),valefits);

shufattr2 = squeeze(vertcat(glm_analysis_periods{1}(1).shufr2(:,wi,attrmodeli,3,:),glm_analysis_periods{2}(1).shufr2(:,wi,attrmodeli,3,:)));
shufvalr2 = squeeze(vertcat(glm_analysis_periods{1}(1).shufr2(:,wi,valmodelis(1)+(valmode-1)*3,3,:),glm_analysis_periods{2}(1).shufr2(:,wi,valmodelis(2)+(valmode-1)*3,3,:)));

issigattrresponsive = false(size(attrefits));
valb(:,valmode) = NaN(size(valefits));
valp(:,valmode) = valb(:,valmode);
for ii=1:length(attrefits)
    tmp = field_or_nil(field_or_nil(attrefits{ii},'stats'),'p',[]);
    if ~isempty(tmp)
        issigattrresponsive(ii) = any(tmp(1:(end-3))<sigthresh);
    end

    tmp = field_or_nil(valefits{ii},'b',[]);
    if ~isempty(tmp)
        valb(ii,valmode) = tmp(1);
    end
    
    tmp = field_or_nil(field_or_nil(valefits{ii},'stats'),'p',[]);
    if ~isempty(tmp)
        valp(ii,valmode) = tmp(1);
    end
end

attr2shufcorr = attr2-mean(shufattr2,2);
valr2shufcorr = valr2-mean(shufvalr2,2);
valcodidx(:,valmode) = clamp(valr2shufcorr./attr2shufcorr,0,1);

x = ps_analysis_periods{1}.epochs(offi).t;

basettl = ['cells'];

thesearevalueinfo = cellfun(@(x)iscell(x) && length(x)==3 && strcmp(x{1},'value') && ismember(7,x{2}{1}.keys),ps_analysis_periods{1}.epochs(offi).attributes);
thesearevalue = cellfun(@(x)iscell(x) && length(x)==3 && strcmp(x{1},'value'),ps_analysis_periods{1}.epochs(offi).attributes);

for doval = 0:2    
    if doval
        nvalbins = 2*(doval-1)+5;
        thesearevaluenvalbins = cellfun(@(x)iscell(x) && length(x)==3 && strcmp(x{1},'value') && x{3} == nvalbins,ps_analysis_periods{1}.epochs(offi).attributes);

        if sum(thesearevaluenvalbins)~=2
            error('expected exactly 2 attribute sets for value analysis with this nvalbins');
        end
        
        valbins = [1:nvalbins]';

        if nvalbins == 5
            valcols = {[0 0 1] [0.5 0 1] [1 0 1] [1 0 0.5] [1 0 0]};
        elseif nvalbins == 7
            valcols = {[0 0 1] [0.33 0 1] [0.67 0 1] [1 0 1] [1 0 0.67] [1 0 0.33] [1 0 0]};
        else
            error('undefined nvalbins');
        end

        conds_ = struct('name','value');
        plotparams = struct(...
            'xname',{'value'},...
            'levelnames',  {arrayfun(@(x)num2str(x),1:nvalbins,'UniformOutput',false)},...
            'colors', {valcols}...
        );
            
    else
        plotparams = struct(...
            'xname',{'R','O','S'},...
            'levelnames',  {{'High','Low'},{'Early','Late'},{'Uncertain',' Certain'}},...
            'colors', {{[1 0 0] [1 0.5 0.5]} {[0 0 1] [0.5 0.5 1]} {[0 0.7 0],[0.6 0.9 0.6]}}...
        );

        halfi_=3;
        thesexname = {'R','O','S'};
        colors = {[1 0 0],[0 0 1],[0 1 0],[0 0 0]};
        
        relevant_ps = psth_bothoffi_relevant_ps{offi,halfi_};
        relevant_bs = psth_bothoffi_relevant_bs{offi,halfi_};
        [uncconds_,uncps,uncbs] = generate_uncconds(thesexname,sigthresh,xname,xnamelabels,this_allu_ok,relevant_ps,relevant_bs,colors);
        conds_ = {};
        for txi=1:length(thesexname)
            condname = xnamelabels(strcmp(xname,thesexname{txi})); condi=strcmp({uncconds_.name},condname); conds_ = [conds_ uncconds_(condi)];
        end

        conds_ = [conds_{:}];
    
        condunits = arrayfun(@(x)num2str(find(x.thisalluok)'),uncconds_,'UniformOutput',false);
    end
    
    for condi=1:length(conds_)
        for ppi=1:length(plotparams)
            if strcmp(plotparams(ppi).xname,'value')
                tmpattr = 'value';
                tmplevels = num2cell(valbins)';
            else
                tmplevels = {1,2};
                switch plotparams(ppi).xname
                    case 'R'
                        tmpattr = {'ev'};
                    case 'O'
                        tmpattr = {'t_out_s'};
                    case 'S'
                        tmpattr = {'entropy'};
                end
            end
            tmpattr = repmat({tmpattr},1,length(tmplevels));
            these_trialconds_ = struct(...
                'name',plotparams(ppi).levelnames, ...
                'attributes',tmpattr, ...
                'levels',tmplevels, ...
                'color',plotparams(ppi).colors);
    
            condname = conds_(condi).name;
            condname = regexprep(condname,'-only','');
            
            if doval
                plotname = plotparams(ppi).xname;
                flipregsel = [];
            else
                plotname = xnamelabels{strcmp(xname,plotparams(ppi).xname)};
                flipregsel = strcmp(xnamelabels,condname);
            end

            if ~(strcmp(condname,plotname) || strcmp(condname,'E[R]') || doval)
                continue;
            else
                if doval==1
                    figname_ = ['figure_6a'];
                elseif doval==2
                    figname_ = ['supp_figure_10c'];
                elseif strcmp(condname,plotname) 
                    figname_ = ['figure_3b_' plotname];
                elseif strcmp(condname,'E[R]')
                    figname_ = ['figure_5b_' plotname];
                else
                    error('figname undefined');
                end
            end

            if strcmp(condname,plotname) && ~doval
                halfi_=[1 2]; halfiplot=[2 1];
            else
                halfi_=3; halfiplot=3;
            end

            dosigsel = ~strcmp(condname,plotname);

            %% calculate stats
            psns_ = [];
            comparisonstatwin_ = [];
            comparisonstat_ = [];
            for hi=1:length(halfi_)
                if doval
                    switch flipmode
                        case 1
                            flipbs = ones(size(valb(:,valmode)));
                            flipmodetxt = 'valnoflip';
                            flipps = NaN(size(flipbs));
                    end
                else
                    psth_relevant_bs = psth_bothoffi_relevant_bs{offi,halfi_(hi)};
                    psth_relevant_ps = psth_bothoffi_relevant_ps{offi,halfi_(hi)};
                    flipbs = psth_relevant_bs(:,flipregsel,offi);
                    flipps = psth_relevant_ps(:,flipregsel,offi);
                end
                if strcmp(condname,'T[R]')
                    flipbs = -flipbs; 
                end
                if strcmp(condname,'Unc[R]') && ~strcmp(condname,plotname) 
                    flipbs(225:end) = -flipbs(225:end);
                end

                tmp = max(size(ps_analysis_periods{1}.epochs(offi).ys,1),size(ps_analysis_periods{2}.epochs(offi).ys,1));
                level = {};
                ys = {};
                for mi_=1:2
                    level{mi_} = ps_analysis_periods{mi_}.epochs(offi).level;
                    ys{mi_} = ps_analysis_periods{mi_}.epochs(offi).ys;
                    if size(ps_analysis_periods{mi_}.epochs(offi).ys,1)<tmp
                        tmpdiff = tmp-size(ps_analysis_periods{mi_}.epochs(offi).level,1);
                        
                        tmp_level = size(ps_analysis_periods{mi_}.epochs(offi).level);
                        tmp_level(1) = tmpdiff;
                        level{mi_} = vertcat(level{mi_},NaN(tmp_level));
                        tmp_ys = size(ps_analysis_periods{mi_}.epochs(offi).ys);
                        tmp_ys(1) = tmpdiff;
                        ys{mi_} = vertcat(ys{mi_},NaN(tmp_ys));
                    end
                    if doval 
                        if mi_==1 
                            tmp = thesearevalue & ~(thesearevaluenvalbins & thesearevalueinfo);
                        else
                            tmp = thesearevalue & ~(thesearevaluenvalbins & ~thesearevalueinfo);
                        end
                    else
                        tmp = find(thesearevalue);
                        tmp(end) = []; 
                    end
                    level{mi_}(:,:,tmp,:) = [];
                end
                t = ps_analysis_periods{1}.epochs(offi).t;
                attributes = ps_analysis_periods{1}.epochs(offi).attributes;
                tmp = cellfun(@(x)iscell(x),attributes);
                attributes(tmp) = []; attributes{end+1} = 'value';
                level = cat(2,level{1},level{2});
                ys = cat(3,ys{1},ys{2});
                noktrall = cat(2,sum(ps_analysis_periods{1}.oktrall),sum(ps_analysis_periods{2}.oktrall));
                
                tmp = length(these_trialconds_);
                psns = NaN(length(allu_ok),length(t),tmp);
                for tci=1:size(psns,3)
                    tcattr = these_trialconds_(tci).attributes;
                    tclevels = these_trialconds_(tci).levels;
                    attri = ismember(attributes,tcattr);
                    for ui=1:length(allu_ok)
                        if this_allu_ok(ui)
                            switch halfiplot(hi)
                                case 1
                                    tmpstart = 2;
                                    tmpby = 2;
                                case 2
                                    tmpstart = 1;
                                    tmpby = 2;
                                case 3
                                    tmpstart = 1;
                                    tmpby = 1;
                            end
                            tmptr = tmpstart:tmpby:noktrall(ui);
                            tmpys = squeeze(ys(tmptr,:,ui));
                            tmplevel = squeeze(level(tmptr,ui,attri,halfiplot(hi)));
                            tmplevel(tmplevel(:,strcmp(attributes(attri),'entropy'))>1,strcmp(attributes(attri),'entropy')) = 2; 
                            tmp = tmpys(all(tmplevel==tclevels,2),:);
                            psns(ui,:,tci) = nanmean(tmp);
                        end
                    end
                end
                psns = (psns - a_unit_norm_param(:,1))./a_unit_norm_param(:,2);

                if strcmp(plotparams(ppi).xname,'S') && ~strcmp(condname,plotname)
                    tmp = flip(psns((224+1):end,:,:),3);
                    psns = vertcat(psns(1:224,:,:),tmp);
                end
                if strcmp(plotname,'T[R]')
                    psns = flip(psns,3);
                end
                comparisonstatwin = NaN(size(this_allu_ok));
                comparisonstat = NaN(length(this_allu_ok),length(ps_analysis_periods{1}.epochs(offi).t));
                switch comparisonstatmode
                    case 'ROC'
                        error('have not edited to reflect new flipping code');
                        for ui=1:length(this_allu_ok)
                            comparisonstatwin(ui) = rocarea3(mean(ys{ui,halfiplot(hi),1}(:,inbounds(x,anawind)),2),mean(ys{ui,halfiplot(hi),2}(:,inbounds(x,anawind)),2));
                            comparisonstat(ui,:) = rocarea3(efilternan(ys{ui,halfiplot(hi),1},'gauss',gauss_smoothing_sigma),efilternan(ys{ui,halfiplot(hi),2},'gauss',gauss_smoothing_sigma));
                        end
                    case 'diff'
                        lowpsn = psns(:,:,1);
                        highpsn = psns(:,:,2);
                        comparisonstat = highpsn-lowpsn;
                        comparisonstatwin = mean(comparisonstat(:,inbounds(x,anawind)),2);
                end
                % flip and NaN out not sig
                nannotsig = ones(size(flipps)); 
                if dosigsel
                    nannotsig(flipps>=sigthresh) = NaN;
                end
                psns = psns.*nannotsig;
                tmp = psns; psns(flipbs<0,:,1) = tmp(flipbs<0,:,2); psns(flipbs<0,:,2) = tmp(flipbs<0,:,1);
                comparisonstatwin = comparisonstatwin.*nannotsig.*sign(flipbs);
                comparisonstat = comparisonstat.*nannotsig.*sign(flipbs);

                % concatenate halves
                psns_ = cat(4,psns_,psns);
                comparisonstatwin_ = cat(4,comparisonstatwin_,comparisonstatwin);
                comparisonstat_ = cat(4,comparisonstat_,comparisonstat);
            end

            psns_ = nanmean(psns_,4);
            comparisonstatwin_ = nanmean(comparisonstatwin_,4);
            comparisonstat_ = nanmean(comparisonstat_,4);
            if doval
                %% value psth
                if doval==1
                    ti = 1;
                else
                    ti = 4;
                end
                dosubtract = 0;
                basewind = [-250 -1];
                tiledefs = struct('name',{'alluok','attr-responsive','attr-responsive & valcodidx>=0.6','E[R]-signaling'},'condok',{this_allu_ok,this_allu_ok & attr2shufcorr>attr2cutoff & issigattrresponsive,this_allu_ok & attr2shufcorr>attr2cutoff & issigattrresponsive & valcodidx(:,valmode) >=0.6,this_allu_ok & ertr_conds(6).thisalluok});
                smoothpsn = efilternan(squeeze(psns),'gauss',gauss_smoothing_sigma);
                smoothpsn = smoothpsn(:,inbounds(x,[-500 1000]),:);
                baseresp = mean(mean(psns(:,inbounds(x,basewind),:),2),3);
                if dosubtract
                    smoothpsn = smoothpsn - baseresp;
                end
                figure; set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.5 0.8])
                smooththensepsn = squeeze(nanstd(smoothpsn(tiledefs(ti).condok,:,:))/sqrt(sum(tiledefs(ti).condok)));
                smooththenmeanpsn = squeeze(nanmean(smoothpsn(tiledefs(ti).condok,:,:),1));    
                for i=1:size(smooththenmeanpsn,2)
                    shadedErrorBar(x(inbounds(x,[-500 1000])),smooththenmeanpsn(:,i),smooththensepsn(:,i),'lineProps',{'color',plotparams.colors{i}});
                end
                axis square;
                title(tiledefs(ti).name);
                xlabel('Time from Offer 1 onset (s)');
                xticklabels(arrayfun(@(x)num2str(x/1000),xticks,'UniformOutput',false));
                ylabel('Normalized activity');
                text(max(xlim)-0.05*range(xlim),max(ylim)-0.05*range(ylim),['n=' num2str(sum(tiledefs(ti).condok))],'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','top');
                tmp = ylim;
                rectheightp = 0.03;
                tmp(1) = tmp(1)-rectheightp*range(tmp);
                ylim(tmp);
                rectangle('Position',[min(anawind) min(ylim) range(anawind) rectheightp*range(tmp)],'FaceColor',[0 0 0]);
                rectangle('Position',[min(basewind) min(ylim) range(basewind) rectheightp*range(tmp)],'FaceColor',[0 0 0]);
                xline(0,'k');

                figname = [figname_ '_psth'];
                savepath=fullfile(neuralsavedir,figname);
                saveas(gcf,[savepath '.png']);
                saveas(gcf,[savepath '.fig']);
                close(gcf);

                %% dots
                dosubtract = doval-1;
                figure; set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.9 0.9]);
                binwinmeans = squeeze(nanmean(psns(:,inbounds(t,these_params.offer.anawind),:),2));
                binbasemeans = squeeze(nanmean(psns(:,inbounds(t,basewind),:),2));
                subtract_this = ternary(dosubtract,binbasemeans,0);
                z = binwinmeans-subtract_this;
                zmean = mean(z(tiledefs(ti).condok,:));
                zerr = std(z(tiledefs(ti).condok,:))/sqrt(sum(tiledefs(ti).condok));
                y_ = reshape(z(tiledefs(ti).condok,:),[],1);
                [~,tmp] = meshgrid(1:sum(tiledefs(ti).condok),1:nvalbins);
                x_ = reshape(tmp',[],1);
                [cc,cp] = corr(x_,y_,'type','Spearman');
                nexttile;
                errorbar(zmean,zerr,'CapSize',0,'Color','k');
                hold on;
                linex(median(1:nvalbins),'Color',0.5*[1 1 1]);
                liney(0,'Color',0.5*[1 1 1]);
                for zi=1:length(zmean)
                    plot(zi,zmean(zi),'Marker','o','MarkerFaceColor',valcols{zi},'MarkerEdgeColor','none','MarkerSize',10);
                end
                xlim([0.5 nvalbins+0.5]);
                text(max(xlim)-0.95*range(xlim),max(ylim)-0.05*range(ylim),['rho=' num2str(cc)],'FontSize',10,'HorizontalAlignment','left','VerticalAlignment','top');
                text(max(xlim)-0.95*range(xlim),max(ylim)-0.1*range(ylim),['p=' num2str(cp)],'FontSize',10,'HorizontalAlignment','left','VerticalAlignment','top');
                axis square;
                title(tiledefs(ti).name);
                ylabel(['Normalized activity' ternary(dosubtract,' (baseline subtracted)','')]);
                xlabel('Integrated value quantile');

                figname = [figname_ '_dots'];
                savepath=fullfile(neuralsavedir,figname);
                saveas(gcf,[savepath '.png']);
                saveas(gcf,[savepath '.fig']);
                close(gcf);

            else
                %% normalized activity
                ttl = [conds_(condi).name ' ' basettl];
                ttl = [monktxt ', ' ttl];
                ttl = [ttl ', normalized activity'];
                ttl = [ttl ', Offer ' num2str(offi)];
                ttl = [ttl ', ' plotname];
                figure('Visible','on');
    
                plot_pref_vs_nonpref(plotparams(ppi).levelnames,this_allu_ok,psns_,x,plotparams(ppi).colors,get_psname(ps_analysis_periods{1}.epochs(offi).ps),gauss_smoothing_sigma,0,offer_plot_t);
    
                title(ttl);
                ttl = [ttl '_' ternary(dosigsel,'dosigsel','alltaskresponsive')];
                if ~dotile
                    figname = [figname_ '_activity'];
                    savepath=fullfile(neuralsavedir,figname);
                    saveas(gcf,[savepath '.fig']);
                    saveas(gcf,[savepath '.png']);
                    close(gcf);
                end
    
                %% comparison stat
                ttl = [conds_(condi).name ' ' basettl];
                ttl = [monktxt ', ' ttl];
                ttl = [ttl ', ' comparisonstatmode];
                ttl = [ttl ', Offer ' num2str(offi)]; 
                ttl = [ttl ', ' plotname];
                figure('Visible','on');
    
                plot_pref_comparison(this_allu_ok,x,plotparams(ppi).colors{1},comparisonstat_,comparisonstatwin_,strjoin(plotparams(ppi).levelnames,' vs '),anawind,comparisonstatmode,get_psname(ps_analysis_periods{1}.epochs(offi).ps),gauss_smoothing_sigma,offer_plot_t);
    
                title(ttl);
                ttl = [ttl '_' ternary(dosigsel,'dosigsel','alltaskresponsive')];
                if ~dotile
                    figname = [figname_ '_comparison'];
                    savepath=fullfile(neuralsavedir,figname);
                    saveas(gcf,[savepath '.fig']);
                    saveas(gcf,[savepath '.png']);
                    close(gcf);
                end
            end
        end
    end
end


%% helper functions
function psname = get_psname(ps)
    if contains(ps,'off')
        psname = regexprep(ps,'off','Offer ');
    else
        psname = [upper(ps(1)) ps(2:end)];
    end
end