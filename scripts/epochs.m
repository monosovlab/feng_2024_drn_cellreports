%%
scriptname = 'epochs';

do_save = 1;
epochfigparams = struct(...
    'name',{'fix on','offer 1','offer 2','choice made','cue','reveal','reward','stimulus off'},...
    'ps',{{'efix'},{'eoff1'},{'eoff2'},{'ecstart'},{'cue'},{'rev'},{'eout'},{'estimoff'}});
plotmodes = {'rpe'}; efps = 5:6;
regtitles = {'Expected reward','Reward delay'};
efit_attr_names = {'R','O'};
cicolors = {[1 0.5 0],[1 0.75 0.5];[0.5 0 1],[0.75 0.5 1]};
cistyles = {'-','--'};
cinames = {{'RPE +','RPE 0','RPE -'}};

for do255025only = 0:1
    if do255025only
        trialconds = struct(...
            'name',{'Info RPE +','Info RPE 0','Info RPE -';'Noinfo RPE +','Noinfo RPE 0','Noinfo RPE -'}, ...
            'attributes',{{'rpe','info','entropy'},{'rpe','info','entropy'},{'rpe','info','entropy'};{'rpe','info','entropy'},{'rpe','info','entropy'},{'rpe','info','entropy'}}, ...
            'levels',{[3 2 3],[2 2 3],[1 2 3];[3 1 3],[2 1 3],[1 1 3]}, ...
            'color',{cicolors{1,1},cicolors{1,2},cicolors{1,1};cicolors{2,1},cicolors{2,2},cicolors{2,1}},...
            'style',{cistyles{1},cistyles{1},cistyles{2};cistyles{1},cistyles{1},cistyles{2}});
        figname = 'supp_figure_8c';
    else
        trialconds = struct(...
            'name',{'Info RPE +','Info RPE 0','Info RPE -';'Noinfo RPE +','Noinfo RPE 0','Noinfo RPE -'}, ...
            'attributes',{{'rpe','info'},{'rpe','info'},{'rpe','info'};{'rpe','info'},{'rpe','info'},{'rpe','info'}}, ...
            'levels',{[3 2],[2 2],[1 2];[3 1],[2 1],[1 1]}, ...
            'color',{cicolors{1,1},cicolors{1,2},cicolors{1,1};cicolors{2,1},cicolors{2,2},cicolors{2,1}},...
            'style',{cistyles{1},cistyles{1},cistyles{2};cistyles{1},cistyles{1},cistyles{2}});
        figname = 'figure_7b';
    end
    allymax = 0;
    trange = {rpe_plot_t};
    flipxname = ''; % if empty, don't flip
    if isempty(flipxname)
        flipbs = ones(size(allu_ok));
    else
        flipbs = sign(relevant_bs(:,strcmp(xname,flipxname),offi));
    end

    condgroups = {[6]};
    modei=1;

    for cgi=1:length(condgroups)
        if isempty(condgroups{cgi})
            conds_ = struct('name','test','color',[0 0 0],'thisalluok',allu_ok & formalGLMoutcome_veryweakpred_ratio);
        else
            conds_ = conds(condgroups{cgi});
        end
        anastats = {};

        basettl = [monktxt ', cue rev'];
        figure('Visible','on'); nrow = length(conds_); ncol = 2;
        tiledlayout(nrow,ncol);
        ylims = zeros(length(conds_),length(efps),2);
        h = {};
        for condi=1:length(conds_)
            thisalluok = conds_(condi).thisalluok;
            for efpsi=1:length(efps)
                cellx = mod(efpsi-1,2)+1;
                efpi = efps(efpsi);
                h{condi,efpsi}=nexttile;
                if sum(thisalluok)==0
                    continue;
                end
                switch plotmodes{modei}
                    case 'rpe'
                        api = efpi-2;
                        
                        tmp = max(size(ps_analysis_periods{1}.epochs(api).ys,1),size(ps_analysis_periods{2}.epochs(api).ys,1));
                        level = {};
                        ys = {};
                        for mi_=1:2
                            level{mi_} = ps_analysis_periods{mi_}.epochs(api).level;
                            ys{mi_} = ps_analysis_periods{mi_}.epochs(api).ys;
                            if size(ps_analysis_periods{mi_}.epochs(api).ys,1)<tmp
                                tmpdiff = tmp-size(ps_analysis_periods{mi_}.epochs(api).level,1);
                                
                                tmp_level = size(ps_analysis_periods{mi_}.epochs(api).level);
                                tmp_level(1) = tmpdiff;
                                level{mi_} = vertcat(level{mi_},NaN(tmp_level));
                                tmp_ys = size(ps_analysis_periods{mi_}.epochs(api).ys);
                                tmp_ys(1) = tmpdiff;
                                ys{mi_} = vertcat(ys{mi_},NaN(tmp_ys));
                            end
                        end
                        t = ps_analysis_periods{1}.epochs(api).t;
                        attributes = ps_analysis_periods{1}.epochs(api).attributes;
                        level = cat(2,level{1},level{2});
                        ys = cat(3,ys{1},ys{2});
                        noktrall = cat(2,sum(ps_analysis_periods{1}.oktrall),sum(ps_analysis_periods{2}.oktrall));

                        if strcmp(epochfigparams(efpi).ps{1},'cue')
                            these_trialconds = trialconds(1,:);
                        else
                            these_trialconds = trialconds(2,:);
                        end
                        
                        psns = NaN(length(allu_ok),length(t),length(these_trialconds));
                        for tci=1:length(these_trialconds)
                            for ui=1:length(allu_ok)
                                tmpys = squeeze(ys(1:noktrall(ui),:,ui));
                                tmplevel = squeeze(level(1:noktrall(ui),ui,ismember(attributes,these_trialconds(tci).attributes),3));
                                tmp = tmpys(all(tmplevel==these_trialconds(tci).levels,2),:);
                                psns(ui,:,tci) = nanmean(tmp);
                            end
                        end
                        psns = (psns - a_unit_norm_param(:,1))./a_unit_norm_param(:,2);
                end
                fullps = psns.*flipbs;
                smoothfullps = efilternan(fullps,'gauss',gauss_smoothing_sigma);
                meanfullps = nanmean(smoothfullps(thisalluok,:,:),1); % added nanmean/nanstd, see above note
                sefullps = nanstd(smoothfullps(thisalluok,:,:),[],1)/sqrt(sum(thisalluok));

                if sum(thisalluok)==1
                    sefullps = zeros(size(sefullps));
                end

                for fpi=1:size(fullps,3)
                    shadedErrorBar(t,meanfullps(1,:,fpi),sefullps(1,:,fpi),'lineProps',{'color',these_trialconds(fpi).color,'linewidth',3,'linestyle',these_trialconds(fpi).style});
                    xline(0,'HandleVisibility','off');
                    hold on;
                end
                axis square;

                if 1 && ismember(efpi,5:6)
                    if modei == 1
                        switch epochfigparams(efpi).name
                            case 'cue'
                                anawind_ = [250 750];
                            case 'reveal'
                                anawind_ = [250 1050];
                        end
                        anastat1 = fullps(:,inbounds(t,anawind_),1)-fullps(:,inbounds(t,anawind_),2);
                        anastat2 = fullps(:,inbounds(t,anawind_),2)-fullps(:,inbounds(t,anawind_),3);
                        anastat = nanmean(cat(3,anastat1,anastat2),3);
                        anastats{efpsi,1} = anastat1;
                        anastats{efpsi,2} = anastat2;
                        anastats{efpsi,3} = anastat;
                    end
                    anastatp = signrank(nanmean(anastat(thisalluok,:),2));
                    text(anawind_(1),min(ylim)+0.1*range(ylim),['p=' num2str(anastatp)],'FontSize',10);
                    rectangle('Position',[min(anawind_) min(ylim) range(anawind_) 0.05*range(ylim)],'FaceColor',[0 0 0]);
                end
                yline(0,'k','HandleVisibility','off');
                ylims(condi,efpsi,:) = ylim;
            end
        end
        %%
        for condi=1:length(conds_)
            thisalluok = conds_(condi).thisalluok;
            if sum(thisalluok)>0    
                for efpsi=1:length(efps)
                    efpi = efps(efpsi);
                    set(gcf,'CurrentAxes',h{condi,efpsi});
                    xlim(trange{modei});
                    if allymax
                        ylim([min(ylims(:,:,1),[],'all') max(ylims(:,:,2),[],'all')]);
                    else
                        ylim([min(ylims(condi,:,1),[],'all') max(ylims(condi,:,2),[],'all')]);
                    end

                    tmp = h{condi,efpsi}.Children(2).Position;
                    h{condi,efpsi}.Children(2).Position = [tmp(1) min(ylim)+0.1*range(ylim) tmp(3)];
                    tmp = h{condi,efpsi}.Children(1).Position;
                    h{condi,efpsi}.Children(1).Position = [tmp(1) min(ylim) tmp(3) 0.05*range(ylim)];

                    xlabel(['Time from ' epochfigparams(efpi).name ' onset (ms)']);
                    if 1
                        ylabel('Normalized activity');
                        if efpsi==1
                            text(max(xlim)-0.05*range(xlim),max(ylim)-0.05*range(ylim),{conds_(condi).name,['(n=' num2str(sum(thisalluok)) ')']},'FontSize',12,'HorizontalAlignment','right','VerticalAlignment','top');
                        end
                        if modei==1
                            title(epochfigparams(efpi).name);
                        end
                    end
                    if condi==1
                        title(ternary(efpsi==1,'Info','Noinfo'));
                        legend(cellfun(@(x)string([ternary(efpsi==1,'Info ','Noinfo ') x]),cinames{modei}),'Location','northwest','FontSize',5);
                    end
                end
            end
        end
        ttl = basettl;
        ttl = [ttl ', Offer ' num2str(offi) ' ' strjoin({conds_.name})];
        ttl = [ttl ', ' plotmodes{modei}];
        savepath=fullfile(neuralsavedir, figname);
        if do_save
            saveas(gcf,[savepath '.fig']);
            saveas(gcf,[savepath '.png']);
            close(gcf);
        end
    end
end