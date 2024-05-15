function [uncconds,uncps,uncbs] = generate_uncconds(thesexname,sigthresh,xname,xnamelabels,this_allu_ok,relevant_ps,relevant_bs,colors)
    uncps = [];
    uncbs = [];
    sigfor = NaN(size(thesexname));
    uncconds = struct('name','Task-responsive','thisalluok',this_allu_ok,'color',[0 0 0],'thesexname',{thesexname},'sigfor',sigfor); % all and marginals first (not none)
    for txi=1:length(thesexname)
        uncps = [uncps relevant_ps(:,strcmp(xname,thesexname{txi}))];
        uncbs = [uncbs relevant_bs(:,strcmp(xname,thesexname{txi}))];
        sigfor = NaN(size(thesexname)); sigfor(txi) = 1;
        newunccond = struct('name',xnamelabels{strcmp(xname,thesexname{txi})},'thisalluok',this_allu_ok & uncps(:,end)<sigthresh,'color',colors{txi},'thesexname',{thesexname},'sigfor',sigfor);
        uncconds = [uncconds newunccond];
        sigfor = NaN(size(thesexname)); sigfor(txi) = 0;
        newunccond = struct('name',['Non-' xnamelabels{strcmp(xname,thesexname{txi})}],'thisalluok',this_allu_ok & ~(uncps(:,end)<sigthresh),'color',colors{txi},'thesexname',{thesexname},'sigfor',sigfor);
        uncconds = [uncconds newunccond];
    end
    
    % joint (not all)
    [~,~,ib] = intersect(thesexname,xname,'stable');
    thesexnamelabels = xnamelabels(ib);
    mutexconds = 2^length(thesexname);
    for condi=1:mutexconds
        bi = dec2bin(condi-1,length(thesexname));
        bi = bi=='1';
        sigfor = bi;
        newunccond = struct('name',get_newunccondname(bi,thesexnamelabels),'thisalluok',this_allu_ok & all((uncps<sigthresh) == bi,2),'color',sum([0 0 0; vertcat(colors{bi})],1),'thesexname',{thesexname},'sigfor',sigfor);
        uncconds = [uncconds newunccond];
    end
end

function newunccondname = get_newunccondname(bi,thesexnamelabels)
    condxnamelabels = thesexnamelabels(bi);
    if sum(bi)==0
        newunccondname = 'None';
    elseif sum(bi)==1
        newunccondname = [condxnamelabels{1} '-only'];
    elseif all(bi)
        newunccondname = 'All';
    else
        newunccondname = strjoin(condxnamelabels,' and ');
    end
end