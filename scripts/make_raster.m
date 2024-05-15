savedir = neuralsavedir;
for api=[1 3 4]
    switch api
        case {1}
            plot_t = offer_plot_t;
            %% select neurons
            us = [10]; uymax = {[0 6.5]}; ucolor = [0 0 0];
            %% task conditions
            conds_ = struct(...
                'name',{'E[R]','T[R]','Unc[R]'},...
                'attribute',{{'ev'},{'t_out_s'},{'entropy'}},...
                'names',{{'Large E[R]','Small E[R]'},{'Short Delay','Long Delay'},{'Uncertain','Certain'}},...
                'condi',{{2,1},{1,2},{2,1}},...
                'color',{{[1 0 0],[1 0.5 0.5]},{[0 0 1],[0.5 0.5 1]},{[0 0.7 0],[0.6 0.9 0.6]}}...
            );
        case {3,4}
            plot_t = rpe_plot_t;
            us = [314]; uymax = {[0 25]}; ucolor = [0 0 0];
            if api==3
                conds_ = struct(...
                    'name',{'infoRPE'},...
                    'attribute',{{'rpe','info'}},...
                    'names',{{'Info RPE +','Info RPE 0','Info RPE -'}},...
                    'condi',{{6,5,4}},...
                    'color',{{[1 0.5 0],[1 0.75 0.5],[1 0.5 0]}},...
                    'linestyle',{{'-','-','--'}}...
                );
            else
                conds_ = struct(...
                    'name',{'noinfoRPE'},...
                    'attribute',{{'rpe','info'}},...
                    'names',{{'Noinfo RPE +','Noinfo RPE 0','Noinfo RPE -'}},...
                    'condi',{{3,2,1}},...
                    'color',{{[0.5 0 1],[0.75 0.5 1],[0.5 0 1]}},...
                    'linestyle',{{'-','-','--'}}...
                );
            end
    end     
    
    %%
    rasIntv=4;
    Hist_template = ps_analysis_periods{1}.epochs(api).t;
    
    for ui=1:length(us)
        u = us(ui);
        if u<=224
            mi_=1;
            au = u;
        else
            mi_=2;
            au = u - 224;
        end
    
        alloktr = ps_analysis_periods{mi_}.oktrall(:,au);
    
        basettl=['unit ' num2str(us(ui)) ' raster and sdf'];
        %%
        for ci=1:length(conds_)

            if api==1
                figname = ['figure_3a_' conds_(ci).name];
            else
                figname = ['figure_7a_' conds_(ci).name];
            end

            ttl = [basettl ' ' strjoin(conds_(ci).names,' vs ')];
            yoffset=0; trn = 0;
            figure('Visible','on');
            tiledlayout(4,1,'TileSpacing','none','Padding','none');
            for li=flip(1:length(conds_(ci).condi))
                if isfield(conds_,'linestyle')
                    linestyle = conds_(ci).linestyle{li};
                else
                    linestyle = '-';
                end
                nexttile(1,[1 1]);
                
                time_start = rasterdat{api}(ci).time_start{li};
                sptimes =  rasterdat{api}(ci).sptimes{li};
    
                [Hist_raster, ax] = Raster_PSTH_plot(sptimes, time_start, Hist_template,gca,conds_(ci).color{li},yoffset,rasIntv,0,35,2);
                yoffset=yoffset+rasIntv*length(time_start);
                
                nexttile(2,[3 1]);
                hold on;
                plot(Hist_template(1:end-1),1000*efilternan(mean(Hist_raster),'gauss',gauss_smoothing_sigma),'Color',conds_(ci).color{li},'LineWidth',3,'LineStyle',linestyle);
                hold on;
                trn = trn+length(time_start);
            end
            nexttile(1,[1 1]);
            hold on;
            ylim([1 rasIntv*trn]+rasIntv/2*[-1 1]);
            xlim(plot_t);
            xticks([]);
            yticks([]);
            title(ttl);
    
            nexttile(2,[3 1]);
            hold on;
            xline(0,'k');
            hold on;

            if iscell(uymax)
                uymax_ = uymax{ui};
            else
                uymax_ = uymax(ui);
            end
            if ~isempty(uymax_)
                if length(uymax_)>1
                    ylim(uymax_);
                else
                    ylim([0 uymax_]);
                end
            end
            xlim(plot_t);
            ylabel('Firing rate (Hz)');
            xlabel(['Time from ' ps_analysis_periods{1}.epochs(api).ps ' onset (ms)']);
            for li=flip(1:length(conds_(ci).condi))
                text(max(xlim)-0.05*range(xlim),max(ylim)-(0.05*(li-1)+0.05)*range(ylim),conds_(ci).names{li},'Color',conds_(ci).color{li},'FontSize',20,'HorizontalAlignment','right','VerticalAlignment','top');
            end

            thisucolor = [0 0 0];

            text(max(xlim)-0.95*range(xlim),max(ylim)-0.05*range(ylim),['Unit ' num2str(u)],'Color',thisucolor,'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','top');

            savepath=fullfile(sigthreshsavedir,figname);
            saveas(gcf,[savepath '.fig']);
            saveas(gcf,[savepath '.png']);
            close(gcf);
        end
    end
end
