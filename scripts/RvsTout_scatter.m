scriptname = 'RvsTout_scatter';

if mi==0
    modeis=1:3;
else
    modeis=1:2;
end

for modei=modeis
    switch modei
        case 1
            selxnames = {'R','O'};
            basecolors = {[1 0 0],[0 0 1]};
            if mi==0
                figname = 'figure_4_left_and_center_hists_figure_5a_left_scatter';
            else
                figname = ['supp_figure_2g_left_and_center_hists_monkey_' monkeyabbr(mi)];
            end
        case 2
            selxnames = {'R','S'};
            basecolors = {[1 0 0],[0 0.7 0]};
            if mi==0
                figname = 'figure_4_left_and_right_hists_figure_5a_right_scatter';
            else
                figname = ['supp_figure_2g_left_and_right_hists_monkey_' monkeyabbr(mi)];
            end
        case 3
            selxnames = {'R','IxS'};
            basecolors = {[1 0 0],[0 0.7 0]};
            figname = 'supp_figure_4c_scatter';
    end
    
    selxnamelabels = get_xnamelabels(selxnames);
    tmp=sum(vertcat(basecolors{:}));
    newColors = [0 0 0; basecolors{1}; basecolors{2}; tmp/max(tmp)];
    
    basettl=[strjoin(selxnamelabels,' vs ') ' Offer ' num2str(offi) ' coding by GLM'];

    %% scatter and hists
    ttl = [monktxt ', ' basettl];
    ttl = [ttl ', with hists'];
    doflip=1;
    xylims = [-0.45 0.45];
    figure('Visible','on');
    [xp,yp] = plot_offer_scatter(befits(1,:),relevant_bs(:,:,offi),relevant_ps(:,:,offi),[],relevant_bs(:,:,offi),relevant_ps(:,:,offi),[],xname,selxnames,selxnamelabels,this_allu_ok,monkeyabbr,1,2,doflip,'GLM',sigthresh,xylims,[],basecolors);
    title(ttl);
    savepath=fullfile(neuralsavedir, figname);
    saveas(gcf,[savepath '.fig']);
    saveas(gcf,[savepath '.png']);
    close(gcf);
end