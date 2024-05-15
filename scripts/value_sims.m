% value_sims.m

%% estimate session by session fluctuation
% use sessionbs from extra_behavior
usbs = {};
usses = {};
for mi_=1:2
    [uniquethissessionbs, tmp] = unique(sessionbs{mi_},'rows','stable');
    uniquethissessionses = sessionses{mi_}(tmp,:);
    usbs{mi_} = uniquethissessionbs';
    usses{mi_} = uniquethissessionses';
end

%% generate sim data
% pseudocode
simattrr2 = repmat({NaN(length(allu_ok),1)},size(hypotheses));
simvalr2 = repmat({NaN(length(allu_ok),1)},size(hypotheses));
simvalb = repmat({NaN(length(allu_ok),1)},size(hypotheses));
simattrb = repmat({NaN(length(allu_ok),4)},size(hypotheses));
simattrp = repmat({NaN(length(allu_ok),4)},size(hypotheses));
simissigattrresponsive = repmat({NaN(length(allu_ok),1)},size(hypotheses));
simshufattrr2 = repmat({NaN(length(allu_ok),nshuf)},size(hypotheses));
simshufvalr2 = repmat({NaN(length(allu_ok),nshuf)},size(hypotheses));
rng(4123,'twister');
for ui = 1:length(allu_ok)
    if ~this_allu_ok(ui)
        continue;
    end
    fprintf(' value simulations, unit %d\n',ui);
    shufinds = nan(attrefits{ui}.n,nshuf);
    for shufi=1:nshuf
        shufinds(:,shufi) = randperm(attrefits{ui}.n);
    end
    for hypi=1:length(hypotheses)
        if ui<225
            mi_ = 1;
            nuisance_xregs = {'Lft','Rgt','2'};
        else
            mi_ = 2;
            nuisance_xregs = {'I','Lft','Rgt','2'};
        end
        % monkey behav stuff
        monkey_xreg = befits{1,mi_}.xreg;
        monkey_xname = befits{1,mi_}.xname;
        monkey_xreg_is_value = ~ismember(monkey_xname,nuisance_xregs);
        monkey_b = befits{1,mi_}.b(monkey_xreg_is_value); % bigfit to behavior, use only value regs
        monkey_sv = befits{1,mi_}.x(:,monkey_xreg_is_value)*monkey_b;
        monk_persession_b_sd = sqrt(clamp(std(usbs{mi_}(monkey_xreg_is_value,:),[],2).^2 - mean(usses{mi_}(monkey_xreg_is_value,:),2).^2,0,inf)); % session by session variability
        addnoise = ternary(strcmp(hypotheses{hypi},'full integration'),monk_persession_b_sd,0);
        multnoise = 0;
        monkey_b_var = (monkey_b + addnoise.*randn(size(monkey_b))) .* (1 + multnoise.*randn(size(monkey_b)));
        
        % some neuron stuff
        neuron_xreg = attrefits{ui}.xreg;
        neuron_xname = attrefits{ui}.xname;
        neuron_xreg_is_value = ~ismember(neuron_xname,nuisance_xregs);
        xcur = attrefits{ui}.x(:,neuron_xreg_is_value);
        val = xcur * monkey_b_var;
        val_to_normact = valefits{ui}.b(1); % valfit,attrfit from glm_analysis_periods, as used in value_analyses
        att_in_normact = xcur * attrefits{ui}.b(neuron_xreg_is_value);
        att_to_normact = std(att_in_normact) ./ std(val); % CHECK THIS
        
        % estimate neural noise based on deviation from
        % attribute model fit (vmodi 1)
        % to original data (run 1)
        noise_sd = std(attrefits{ui}.stats.resid);
        noise = noise_sd*randn(size(attrefits{ui}.y));   

        switch hypotheses{hypi}
            case 'single attribute'
                atti = 1+mod(ui-1,numel(monkey_b));
                attb = max(abs(monkey_b(monkey_xreg_is_value)));
                hypb = zeros(size(monkey_b));
                hypb(atti) = attb;
                signal_to_normact = att_to_normact;
            case 'random mixture'
                hypb = monkey_b_var(randperm(numel(monkey_b_var))) .* sign(rand(size(monkey_b_var))-0.5);
                signal_to_normact = att_to_normact;
            case 'partial integration'
                n_half_weights = ceil(numel(monkey_b_var)*0.5);
                weights_to_randomize = [ones(n_half_weights,1) ; zeros(numel(monkey_b_var) - n_half_weights,1)] == 1;
                weights_to_randomize = weights_to_randomize(randperm(numel(weights_to_randomize)));
                
                % sign-flip and shuffle randomized
                % weights
                b_mixed_weights = monkey_b_var;
                b_rand_weights = b_mixed_weights(weights_to_randomize);
                b_rand_weights = b_rand_weights(randperm(numel(b_rand_weights))) .* sign(rand(size(b_rand_weights))-0.5);
                b_mixed_weights(weights_to_randomize) = b_rand_weights;
                hypb = b_mixed_weights;
                signal_to_normact = att_to_normact;
            case 'full integration'
                hypb = monkey_b_var;
                signal_to_normact = val_to_normact;
        end
        coded_var = xcur*hypb*signal_to_normact;
        cury = coded_var + noise;

        %% glm
        % fit simdata to attr and val models
        % 1 unshuf, 10 shuf
        % just need to save r2 for all and valb,issigattrresponsive for unshuf
        for shufi=0:nshuf
            if shufi
                inds = shufinds(:,shufi);
            else
                inds = 1:size(shufinds,1);
            end
            attrefit_ = calcr2(eglm_fit(neuron_xreg,cury(inds),'normal'));
            valefit_ = calcr2(eglm_fit(valefits{ui}.xreg,cury(inds),'normal'));

            if shufi
                simshufattrr2{hypi}(ui,shufi) = attrefit_.r2;
                simshufvalr2{hypi}(ui,shufi) = valefit_.r2;
            else
                simattrr2{hypi}(ui) = attrefit_.r2;
                simvalr2{hypi}(ui) = valefit_.r2;
                simissigattrresponsive{hypi}(ui) = any(attrefit_.stats.p(neuron_xreg_is_value));
                [~,~,ii] = intersect({'R','O','S','IxS'},attrefit_.xname,'stable');
                simattrb{hypi}(ui,:) = attrefit_.b(ii);
                simattrp{hypi}(ui,:) = attrefit_.stats.p(ii);
                simvalb{hypi}(ui) = valefit_.b(1);
            end
        end
    end
end
