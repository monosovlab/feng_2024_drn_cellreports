function fit = eglm_fit(xreg, yreg, fitting_style, varargin)
% fit = eglm_fit(xreg, yreg, fitting_style[,...])
%
% Fit a GLM to data, using code originally meant for binary choices between
% offers that have multiple features that can contribute to the decision.
% The main inputs are data structures specifying the x (regressors) and y
% (variabled to be predicted), and the type of GLM to use.
%
% Note that this does NOT automatically add a regressor that is constant.
% If you want a constant regressor, you need to add it to the list of
% regressors in xreg. This is because in our usual type of models we
% already have a constant regressor built in (e.g. main effect of offer
% presentation order)
%
% inputs:
%  xreg - a Nregressors-length cell array of specs for each regressor. Each
%         entry xi should be a data structure with the following field:
%
%         xreg{xi}.terms: a cell array where each entry is either one of
%          the individual underlying terms whose product will be the 
%          regressor, or a 'flag' variable indicating that the following
%          terms should be treated as flag variables. The terms should all
%          be either (Ntrials x 1) vectors, or (Ntrials x 2) matrices where
%          the first and second columns represent the properties of the
%          first and second offers on that trial.
%
%         Optionally, it can have the following fields:
%
%         xreg{xi}.normalization: specifies what normalization style to
%          use. Options are:
%
%           'default' - apply z-scoring normaliation to all terms except the
%            flag variables.
%
%           'only interaction terms' - same as default, but do not 
%            normalize the first term in the list of terms. So all effects
%            will be essentially in units of 'effect of change in first 
%            listed variable, PER-STANDARDIZED-CHANGE in the other listed 
%            variables'.
%
%           'none' - do not normalize any of the terms.
%
%           {'oktrials', Z} - a vector Z (logical or numeric) of indexes of 
%            trials. Same effect as default, but when computing the mean 
%            and SD for the z-scoring normalization of each term, uses all 
%            trials specified by Z, NOT the 'oktrials' vector that is 
%            normally used to restrict the analysis.
%
%            This option is meant for cases where e.g. we are doing
%            separate regressions for each of 20 separate behavioral
%            sessions, but in each of those fits we want to normalize the
%            terms of each regressor based on their overall distribution in
%            ALL 20 of the sessions (to ensure that they are perfectly
%            comparable across sessions).
%
%         xreg{xi}.manually_specified_beta: forces the regression weight to
%          be a specific value. That beta weight's SE and p-value will be
%          set to NaN. Use this when you have a set of models in the same
%          general framework, but some of them operate by forcing specific
%          parameters to be fixed weights rather than freely fit weights.
% 
%  yreg - the variable to be predicted. Either a (Ntrials x 1) column
%         vector of its values on each trial, or a data structure with 
%         that stored in its yreg.value field.
%
%  fitting_style - specifies the type of GLM. Can be either a 2-element
%         cell array of strings in the format
%         {distribution_name,link_function_name}, or else a string
%         specifying one of the following common fitting styles for our
%         data:
%          'binary choices' = {'binomial','logit'}
%          'poisson' = {'poisson','identity'}
%          'normal' = {'normal','identity'}
%
% Optional inputs, specified as a string indicating the name of the input
%  parameter, followed by its value:
%
%  'name',name - give the fit a name. This will be stored in fit.name. May
%         be useful for reference e.g. when you have a large
%         collection of fits.
%
%  'oktrials',oktrials - a (Ntrials x 1) logical vector or a (k x 1) 
%         numeric vector (where k <= Ntrials) specifying a subset of 
%         trials. The analysis will be restricted to those trials.
%
%  'calculate predictions',mode - specify whether to calculate predicted
%         y values for each trial, log likelihoods of the data, and related
%         values. mode must be logical (true or false).
%
%  'rescale betas' - rescales the beta weights and their SEs so that
%         instead of being in the raw units of the predicted variable
%         ('s transformation by the link function), they are in units
%         of 'effect of changing this regressor by +1, AS A FRACTION OF 
%         the effect that would occur by changing the regressor named 
%         xreg_name by an amount equal to xreg_change'. For instance, when
%         predicting behavior in a 'pay per view' experiment, you may want
%         to estimate the subject's willingness to pay to obtain each 
%         regressor, so you may want to scale the beta weights to be in 
%         units of 'effect of +1 mL offered juice' or 'effect of +$1
%         offered money'. Specify this as a struct with the fields:
%
%          'name' - name of the new units (e.g. 'Effect on log odds of
%           choice relative to +1 mL juice')
%          'xreg_name' - name of the regressor to use for rescaling (e.g.
%           'E[Reward]')
%          'xreg_change' - number indicating the change in the regressor to
%           use for setting the amount of rescaling (e.g. if the regressor
%           is in un-normalized units of uL juice, then you could use 
%           xreg_change = 1000 to scale beta weights to be relative to the
%           effect of 1000 x +1 uL = +1 mL of juice)
%         
%         The name of the 
%
% output:
%  fit - a data structure with the fitting results. Fields include:
%
%  fit.warning - information about whether a warning was detected while
%   calling glmfit to fit the GLM. Will print out a message if a warning is
%   detected, in which case you should be wary of trusting the results!
%
%  fit.name - name of the fit, if specified by the user
%
%  fit.param - basic parameters of the fit, including:
%
%  fit.param.distribution - GLM distribution function
%  fit.param.link - GLM link function
%
%  fit.param.oktrials - which trials in the data were used to fit the model
%  fit.param.normalization{xi} - parameters use to normalize the xi-th
%   regressor, including whether it is normalized, and the mean & sd that
%   were used to do the z-scoring normalization
%
%  fit.param.betascale - parameters used to do scaling of the beta weights
%   (if scaling was done). 
%  fit.param.betascale.name - always holds the name of the units of the
%   beta weights (e.g. "Effect on logs odds of choosing offer 2")
%
%  fit.xreg, fit.yreg - data specifying x (regressors) 
%   and y (variable to be predicted)
%
%  fit.x, fit.y - regressors (x) and variable to be predicted (y)
%  fit.b, fit.stats - outputs of glmfit
%
%  fit.glmfun_inputs - data structure holding the raw inputs given to
%   glmfit / glmval, which are the same as fit.x and fit.b in typical
%   cases, but may be different if you set certain regressors to have
%   betas that are manually-specified instead of fitted by the GLM.
%
%  The GLMs predictions (actually, in the current implementation, 
%  postdictions) include:
%
%  fit.ypred - predicted value of y on each trial
%  fit.loglik_tr - log likelihood of y on each trial
%  fit.loglik - total log likelihood of the entire dataset
%  fit.mean_loglik_per_trial - log likelihood divided by # of data points
%  fit.p_correct_prediction - for binomial data, the fraction of trials
%   where the model predicts y 'correctly' (i.e. with probability > 0.5).
%
%  fit.calibration - data structure holding results of 'calibration', i.e.
%   splitting the data into bins with different predicted values of y, and
%   testing whether the mean value of y in each group is close to its
%   predicted value. Key fields are:
%  fit.calibration.nbin - number of bins
%  fit.calibration.ypredmean - each bin's mean predicted value 
%  fit.calibration.ymean - each bin's mean observed value
%  fit.calibration.ysd - SD of each bin's observed value
%  fit.calibration.yse - SE of each bin's observed value
%
% written by ESBM


fit = struct();
fit.name = '';

% whether to calculate predictions for each trial, 
% and the log likelihood of the data given the fitted model
fit.param.calculate_predictions = true;

% whether to use a subset of the original data
fit.param.oktrials = 'default';

% whether to rescale the fitted beta weights so that they represent the
% effect of a specific amount of change in a specific xregressor's
% raw data values (e.g. scaling from "log odds of biasing choice" to
% "equivalent effect of adding +1 mL of juice reward")
fit.param.betascale.name = [];
fit.param.betascale.xreg_name = [];
fit.param.betascale.xreg_change = [];
fit.param.betascale.rescaling_factor = [];

% read in variable arguments.
vi = 1;
while vi <= numel(varargin)
    switch varargin{vi}
        case 'name'
            assert(numel(varargin) >= vi+1 && ischar(varargin{vi+1}),'variable argument "%s" must be followed by a string',varargin{vi});
            fit.name = varargin{vi+1};
            vi = vi + 1;
        case 'calculate predictons'
            assert(numel(varargin) >= vi+1 && islogical(varargin{vi+1}),'variable argument "%s" must be followed by a logical',varargin{vi});
            fit.param.calculate_predictions = varargin{vi+1};
            vi = vi + 1;
        case 'oktrials'
            assert(numel(varargin) >= vi+1 && (isnumeric(varargin{vi+1}) || islogical(varargin{vi+1})),'variable argument "%s" must be followed by a numeric or logical',varargin{vi});
            fit.param.oktrials = varargin{vi+1};
            vi = vi + 1;
        case 'rescale betas'
            assert(numel(varargin) >= vi+1 && isstruct(varargin{vi+1}) && all(isfield(varargin{vi+1},{'name','xreg_name','xreg_change'})),'variable argument "%s" must be followed by a data structure with fields "name", "xreg_name", "xreg_change"',varargin{vi});
            fit.param.betascale = varargin{vi+1};
            vi = vi + 1;
        otherwise
            error('unknown variable argument "%s"',varargin{vi});
    end
    vi = vi + 1;
end

% convert raw regressor inputs into actual regressors
fit.xreg = xreg;
fit.yreg = yreg;
[fit.x,fit.y,fit.param.normalization,fit.glmfun_inputs] = eglm_make_regression_inputs(fit.xreg,fit.yreg,fit.param.oktrials);
fit.b_is_manually_specified = ~isnan(fit.glmfun_inputs.manually_specified_beta);
fit.nx = size(fit.x,2);
fit.n = numel(fit.y);

% get names of x and y variables
fit.xname = cellfun(@(z) z.name,fit.xreg,'uniform',0);
if isstruct(fit.yreg)
    fit.yname = fit.yreg.name;
else
    fit.yname = 'predicted variable';
end

% fitting style specifies the probability distribution the data is assumed
% to follow, and the link function. It is either specified as:
% {distribution,link function}, or shorthand for a commonly used such pair.
if iscell(fitting_style)
    assert(numel(fitting_style) == 2 && ischar(fitting_style{1}) && ischar(fitting_style{2}),'if fitting_style is a cell array, must have two strings indicating the distribution and the link function');
else
    switch fitting_style
        case 'binary choices'
            fitting_style = {'binomial','logit'};
        case 'poisson'
            fitting_style = {'poisson','identity'};
        case 'normal'
            fitting_style = {'normal','identity'};
        otherwise
            error('unknown fitting style');
    end
end
fit.param.distribution = fitting_style{1};
fit.param.link = fitting_style{2};

% track whether the fit gives a warning, and if so, what type
fit.warning.detected = false;
fit.warning.msg = [];
fit.warning.id = [];
[lastmsg_before_fit,lastid_before_fit] = lastwarn;
lastwarn('');

% fit the model
% [fit.glmfun_inputs.b,~,fit.stats] = glmfit(fit.glmfun_inputs.x,fit.y,fit.param.distribution,'link',fit.param.link,'constant','off','offset',fit.glmfun_inputs.offset);

% call glmfit with old syntax to avoid massive slowdown due to glmfit using
% try/error/catch to check if 1st variable argument is in a list of strings
% to determine if it is being called with the new syntax
% argument order is: x, y, distribution, link function, estdisp, offset, weights, constant, rankwarn, options, b0
[fit.glmfun_inputs.b,~,fit.stats] = glmfit(fit.glmfun_inputs.x,fit.y,fit.param.distribution,fit.param.link,'off',fit.glmfun_inputs.offset,[],'off');

% adjust fit so that our variables reflect the ones that were manually
% specified to be fixed values
fit = eglm_adjust_fit_for_manually_specified_beta(fit);


% check if warnings were given
[fit.warning.msg,fit.warning.id] = lastwarn;
fit.warning.detected = ~isempty(fit.warning.msg) || ~isempty(fit.warning.id);
if fit.warning.detected
    fprintf(' eglm_fit detected at least one warning when calling glmfit to fit the GLM! Be wary of trusting the fitting results!\n');
    fprintf('  last warning msg: %s\n  last warning  id: %s\n',fit.warning.msg,fit.warning.id);
end
% if no warning was detected, restore prior state of Matlab's warning flags
if ~fit.warning.detected
    lastwarn(lastmsg_before_fit,lastid_before_fit);
end



% if required, calculate model's prediction about each trial
% and calculate the log likelihood
if fit.param.calculate_predictions
    % calculate model's prediction about each data point
    fit.ypred = glmval(fit.glmfun_inputs.b,fit.glmfun_inputs.x,fit.param.link,'constant','off','offset',fit.glmfun_inputs.offset);

    % calculate log likelihood of each individual trial
    % based on the model's fitted parameters and its assumption about
    % the distribution the data follow
    fit.loglik_tr = nans(size(fit.ypred));
    switch fit.param.distribution
        case 'binomial'
            fit.loglik_tr = log((fit.y==1).*fit.ypred + (fit.y==0).*(1-fit.ypred));
        case 'poisson'
            fit.loglik_tr = log(poisspdf(fit.y,fit.ypred));
        case 'normal'
            fit.loglik_tr = log(normpdf(fit.y,fit.ypred,fit.stats.sfit));
        otherwise
            error('have not implemented calculation of log likelihood for GLM distribution "%s"',fit.param.distribution);
    end

    % total log likelihood of entire dataset
    fit.loglik = sum(fit.loglik_tr);

    % mean log likelihood per trial
    fit.mean_loglik_per_trial = fit.loglik ./ fit.n;
    
    % if binomial data, calculate p(correct prediction)
    if strcmp(fit.param.distribution,'binomial')
        fit.p_correct_prediction = mean(sign(fit.ypred-0.5) == sign(fit.y-0.5));
    else
        fit.p_correct_prediction = nans(fit.n,1);
    end
end

% if required, rescale beta weights and their SEs (but no other stats)
% to be in units based on a fixed change in an original raw variable
% (e.g. "relative to the effect of adding +1 mL of juice")
if ~ischar(fit.param.betascale.name)
    % no scaling. Set the name of the beta weight units based on the name
    % of the y variable and its transformation by the link function
    curyname = fit.yname;
    if strcmp(fit.param.distribution,'binomial')
        if strcmp(fit.param.link,'logit')
            curyname = ['log odds of ' curyname];
        else
            curyname = ['p(' curyname ')'];
        end
    end
    
    fit.param.betascale.name = ['Effect on ' curyname];
else
    % scaling needs to be done
    xi = find(strcmp(fit.xname,fit.param.betascale.xreg_name));
    assert(numel(xi)==1,'did not find unique regressor with name "%s" to do rescaling of beta weights',fit.param.betascale.name);
    
    % get effect of adding the specified change to the regressor
    normparam = fit.param.normalization{xi};
    nterms = numel(normparam.sd);
    assert(nterms==1,'cannot rescale beta weights using regressor "%s", can only rescale beta weights based on a regressor that was computed using a single term (not the product of multiple terms)',fit.param.betascale.name);
    
    % get SD used to normalize the regressor
    nsd = normparam.sd;
    % if was not normalized, we treat it as if normalization SD = 1
    % (i.e., no normalization by SD)
    if isnan(normparam.sd)
        nsd = 1;
    end
    
    % fitted beta (in units of change in predicted variable per addition 
    %  of +1 to the normalized variable)
    b_normalized = fit.b(xi);
    % fitted beta (in units of the change in predicted variable per
    %  addition of +1 to the ORIGINAL variable)
    b_orig = b_normalized ./ nsd;
    
    % rescaling factor (i.e. factor to multiply all betas so that are in
    %  units of the equivalent change in predicted variable per addition of
    %  +xreg_change to the ORIGINAL variable)
    fit.param.betascale.rescaling_factor = 1 ./ (b_orig .* fit.param.betascale.xreg_change);
    
    % apply the rescaling
    fit.b = fit.b .* fit.param.betascale.rescaling_factor;
    fit.stats.beta = fit.stats.beta .* fit.param.betascale.rescaling_factor;
    fit.stats.se = fit.stats.se .* fit.param.betascale.rescaling_factor;
end