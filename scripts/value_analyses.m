%% value coding index
scriptname = 'value_analyses';

do_plot = 1;
valcodidx = NaN(length(allu_ok),4);
valb = valcodidx;
valp = valcodidx;

hypotheses = {'single attribute','random mixtures','partial integration','full integration'};
valmodelnames = {'attr','with_Info','without_Info'};
wi = 1;
valmodelis = [2 3];
attrmodeli = 1;
nshuf = 100;
%%
for valmode=[1 3 4]
    switch valmode
        case 1
            vmtxt = 'Linear value, all behavior from recording';
        case 3
            vmtxt = 'Linear value, each unit''s session';
        case 4
            vmtxt = 'Prelec-2 value, all behavior from recording';
    end

    attrefits_ = vertcat(glm_analysis_periods{1}(1).efits(:,wi,attrmodeli,3),glm_analysis_periods{2}(1).efits(:,wi,attrmodeli,3));
    valefits_ = vertcat(glm_analysis_periods{1}(1).efits(:,wi,valmodelis(1)+(valmode-1)*3,3),glm_analysis_periods{2}(1).efits(:,wi,valmodelis(2)+(valmode-1)*3,3));

    attrefits = attrefits_;
    valefits = valefits_;

    if valmode==1
        value_sims;
        hypis = 0:length(hypotheses);
    else
        hypis = 0;
    end
%%
    for hypi=hypis
        if hypi
            hyptxt = hypotheses{hypi};
            attr2 = simattrr2{hypi};
            valr2 = simvalr2{hypi};
            shufattr2 = simshufattrr2{hypi};
            shufvalr2 = simshufvalr2{hypi};
            issigattrresponsive = simissigattrresponsive{hypi};
            valb = simvalb{hypi};
        else
            hyptxt = 'real';
            attrefits = attrefits_;
            valefits = valefits_;
            attr2 = cellfun(@(x)field_or_nil(x,'r2'),attrefits);
            valr2 = cellfun(@(x)field_or_nil(x,'r2'),valefits);
    
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
        end

        attr2shufcorr = attr2-mean(shufattr2,2);
        valr2shufcorr = valr2-mean(shufvalr2,2);
        valcodidx(:,valmode) = clamp(valr2shufcorr./attr2shufcorr,0,1);

        thisuok = this_allu_ok & attr2shufcorr>attr2cutoff & issigattrresponsive==1;
    
        thisvalr2shufcorr = valr2shufcorr;
        thisvalcodidx = valcodidx(:,valmode);
        thisvalb = valb(:,valmode);

        thisbinwidth = 0.2;

        if do_plot
            figure; tcl = tiledlayout(1,2); set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.6 1]);
            for plotmodei = 1:2
                nexttile();
                switch plotmodei
                    case 1
                        plotmodetext = 'hist';
                        histogram(thisvalcodidx(thisuok),'BinWidth',thisbinwidth,'FaceColor',structure_color,'Normalization','probability');
                        xlim([0 1]);

                        xlabel('Value coding index');
                        ylabel('% of neurons with strong attribute effects');
                        set(gca,'box','off');
                        axis square;
                    case 2
                        plotmodetext = 'scatter';
                        scatter(attr2shufcorr(~thisuok & this_allu_ok,end),thisvalr2shufcorr(~thisuok & this_allu_ok,end),[],0.7*[1 1 1],'.');
                        hold on;
                        scatter(attr2shufcorr(thisuok,end),thisvalr2shufcorr(thisuok,end),20,structure_color);
                        leftbotlim = min(-attr2cutoff,min([attr2shufcorr(thisuok,end) thisvalr2shufcorr(thisuok,end)],[],'all'));
                        xlim([leftbotlim,1]);
                        ylim([leftbotlim,1]);
                        xline(0,'k');
                        yline(0,'k');

                        xlabel('Attribute model r_{corr}^2');
                        ylabel('Value model r_{corr}^2');

                        uistack(line(xlim,ylim,'Color',0.7*[1 1 1]),'bottom');
                        linex(attr2cutoff,'Color',0*[1 1 1],'LineStyle','--');
                        axis square;
                end
                if 1
                    title(escape({['(n=' num2str(sum(thisuok)) '/' num2str(sum(this_allu_ok)) ')'],'Note: May look slightly different from published figures','given the non-deterministic nature of the shuffle procedure','used to compute the value coding index'}));
                end
            end
        end

        if do_plot
            if hypi==0
                if valmode==1
                    figname = ['figure_6bc_supp_figure_5_supp_figure_6c_realdata_' vmtxt];
                else
                    figname = ['supp_figure_6c_realdata_' vmtxt];
                end
            else
                figname = ['supp_figure_5_' hypotheses{hypi} '_' vmtxt];
            end
            savepath=fullfile(neuralsavedir, figname);
            saveas(gcf,[savepath '.png']);
            saveas(gcf,[savepath '.fig']);

            close(gcf);
        end
        
        if hypi==0
            %% choice predictive index
            thisoffi = 1;
            choice_predictive_index = NaN(size(valefits,1),2,2);
            for ui = 1:length(valefits)
                if ~allu_ok(ui)
                    continue;
                end
                if ui<(224+1)
                    mi_=1;
                    unit_i = ui;
                else
                    mi_=2;
                    unit_i = ui-224;
                end
                x = attrefits{ui}.stats.resid; y = 2*(thisoffi-1.5)*unityresids{mi_}{unit_i}; oktr = abs(y)>=0.05;
                x = x(oktr); y = y(oktr);
            
                [choice_predictive_index(ui,1,1), choice_predictive_index(ui,2,1)] = corr(x,y,'type','Pearson');
                [choice_predictive_index(ui,1,2), choice_predictive_index(ui,2,2)] = corr(x,y,'type','Spearman');
                
            end
            choice_predictive_index = choice_predictive_index.*sign(valb(:,valmode));
            
            %%
            figure; set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.5 0.8]);
            tcl = tiledlayout(2,2);
            
            for cti=1:2
                switch cti
                    case 1
                        corrtype = 'pearson';
                    case 2
                        corrtype = 'spearman';
                end
                for tauoki=1
                    switch tauoki
                        case 1
                            thisuok = this_allu_ok;
                            oktxt = 'task-responsive';
                    end
                    nexttile();
                    x = choice_predictive_index(:,1,cti); xp = choice_predictive_index(:,2,cti);
                    
                    xline(0,'k');
                    hold on;
            
                    binwidth = 0.05;
                    tmp = histogram(x(thisuok),'FaceColor',0.5*[1 1 1],'BinWidth',binwidth); 
                    xline(mean(x(thisuok)),'color',0.5*[1 1 1],'LineStyle','--');
                    srp = signrank(x(thisuok));
            
                    histogram(x(thisuok & xp<sigthresh),'FaceColor',0*[1 1 1],'BinWidth',binwidth);
                    xline(mean(x(thisuok & xp<sigthresh)),'color',0*[1 1 1],'LineStyle','--');
                    sigsrp = signrank(x(thisuok & xp<sigthresh));
                    ns = sum((x*[-1 1] > 0) & thisuok & xp<sigthresh);
                    binop = myBinomTest(ns,sum(thisuok),sigthresh/2);
            
                    xlim([-1.3 1.3]*max(abs(xlim)));
                    
                    text(max(xlim)-0.05*range(xlim),max(ylim)-0.05*range(ylim),['signed-rank p =' num2str(srp)],'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','top');
                    text(max(xlim)-0.05*range(xlim),max(ylim)-0.15*range(ylim),['signed-rank p =' num2str(sigsrp)],'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','top');
                    text(max(xlim)-0.95*range(xlim),max(ylim)-0.90*range(ylim),[char(string(round(100*ns(1)/sum(thisuok),1))),'%'],'FontSize',10,'HorizontalAlignment','left','VerticalAlignment','bottom');
                    text(max(xlim)-0.95*range(xlim),max(ylim)-0.95*range(ylim),['p = ', char(string(round(binop(1),2)))],'FontSize',10,'HorizontalAlignment','left','VerticalAlignment','bottom');
                    text(max(xlim)-0.05*range(xlim),max(ylim)-0.90*range(ylim),[char(string(round(100*ns(2)/sum(thisuok),1))),'%'],'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','bottom');
                    text(max(xlim)-0.05*range(xlim),max(ylim)-0.95*range(ylim),['p = ', char(string(round(binop(2),2)))],'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','bottom');
                    
                    ylabel('# of neurons');
                    xlabel('Choice Predictive Index');
                    title([oktxt ', ' corrtype ' choice predictive index (n=' num2str(sum(thisuok)) ')']);
                end
            end
            
            figname = ['supp_figure_6ab_' vmtxt];
            savepath=fullfile(neuralsavedir, figname);
            saveas(gcf,[savepath '.png']);
            saveas(gcf,[savepath '.fig']);
            close(gcf);

            if valmode==1
                %% valcoding vs rpe
                valcodingcutoff = 0.6;
                
                thisuok = this_allu_ok & attr2shufcorr>attr2cutoff & issigattrresponsive;
                isvalcoding = valcodidx >= valcodingcutoff;
                thisanastat = mean(anastatfull(:,1:2,end),2);
                thisanastatp = combine_pvalues(anastatpfull(:,1:2,end),2);
                
                figure; tcl = tiledlayout(1,2);

                nexttile;
                h = gca;
                A = [sum(isvalcoding(thisuok)),sum(isvalcoding(thisuok) & thisanastatp(thisuok)<sigthresh & sign(valb(thisuok))==sign(thisanastat(thisuok)-0.5)),sum(thisanastatp(thisuok)<sigthresh)];

                [H,S] = venn_(A([1 3]),A(2));
                h = h.Children([end end-1]);

                axis square;
                axis equal;
                axis off;
               
                %Now label each zone 
                for zi = 1:3
                    switch zi
                        case 1
                            tmpx = mean([min(h(1).XData) min(h(2).XData)]);
                            tmpn = A(1)-A(2);
                            binop = myBinomTest(tmpn,sum(thisuok),A(1)/sum(thisuok));
                        case 2
                            tmpx = mean([max(h(1).XData) min(h(2).XData)]);
                            tmpn = A(2);
                        case 3
                            tmpx = mean([max(h(1).XData) max(h(2).XData)]);
                            tmpn = A(3)-A(2);
                    end
                    tmptxt = {[num2str(tmpn) '/' num2str(sum(thisuok))] [num2str(round(100*tmpn/sum(thisuok))) '%']};
                    text(tmpx,S.ZoneCentroid(1,2),tmptxt,'HorizontalAlignment','center','VerticalAlignment','middle');
                end
                text(prctile(h(1).XData,30),max(ylim),{'Strong' 'value coding'},'Color',0.8*h(1).FaceColor,'HorizontalAlignment','center','VerticalAlignment','top');
                text(prctile(h(2).XData,70),max(ylim),'RPE','Color',0.8*h(2).FaceColor,'HorizontalAlignment','center','VerticalAlignment','top');
                
                nexttile;
                x = sign(valb(thisuok)).*valcodidx(thisuok); 
                y = thisanastat(thisuok);
                myScatterWithHistogram('Signed Value Coding Index','RPE Index',x,y,[],[],[],[],[],{[-1.02,1.02],[0,1]},0,0.5,[],sigthresh,1,[],[],'none');
                axis square;

                title(tcl,{'Note: May look slightly different from published figures','given the non-deterministic nature of the shuffle procedure','used to compute the value coding index'});

                figname = ['supp_figure_10ab'];
                savepath=fullfile(neuralsavedir,figname);
                saveas(gcf,[savepath '.png']);
                saveas(gcf,[savepath '.fig']);
                close(gcf);
            end
        end
    end
end