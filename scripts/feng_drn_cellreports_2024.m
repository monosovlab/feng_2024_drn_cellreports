% feng_drn_cellreports_2024.m

% instructions for running this script:
% 1. download .mat data files from:
%     https://wustl.box.com/s/8i2tbnvt16wlu8fsp2xymrn175ryyubd
% 2. put the downloaded .mat data files in this location: 
%     feng_2024_drn_cellreports/data
% 3. open matlab and set your current directory to the scripts folder:
%     feng_2024_drn_cellreports/scripts
% 4. run this script


do_load = 1;
do_run = 1;

if do_load

    basepath = pwd();
    scripts_suffix = [filesep() 'scripts'];
    if endsWith(basepath,scripts_suffix)
        basepath = basepath(1:(end-numel(scripts_suffix)));
        fprintf('detected that current directory is the "scripts" subdirectory\n setting base path to "%s"\n',basepath);
    else
        fprintf('setting base path to current directory "%s"\n',basepath);
    end
    expected_base_path_suffixes = {'feng_2024_drn_cellreports','feng_2024_drn_cellreports-main'};
    assert(any(cellfun(@(z) endsWith(basepath,z),expected_base_path_suffixes)),'base path was "%s", was expected to end with %s or %s',basepath,expected_base_path_suffixes{1},expected_base_path_suffixes{2});

    fprintf('setting up paths and parameters\n');

    addpath(basepath);
    addpath(fullfile(basepath,'data'));
    addpath(fullfile(basepath,'scripts'));
    addpath(fullfile(basepath,'scripts','MATLAB'));
    addpath(fullfile(basepath,'scripts','util'));
    savebasepath = basepath;


    %% set params
    figbasedir = fullfile(savebasepath,'figures');
    fpbasedir = figbasedir;
    sigthreshsavedir = figbasedir;
    neuralsavedir = figbasedir;

    structure = 'DRN';
    structure_color = [1 0.5 0.5];
    monkeyabbr = ['Z','P'];
    these_params = struct();
    offer_plot_t = [-500 1000];
    anticip_plot_t = [-1500 500];
    rpe_plot_t = [-500 1500];
    gauss_smoothing_sigma = 50;
    attr2cutoff = 0.03;
    task = 'IMTv11';
    
    load('these_params.mat');
    
    sigthresh = these_params.sigthresh; 
    anawind = these_params.offer.anawind;
    
    % make directories to store the figures
    dirs_to_create = {figbasedir,fpbasedir,sigthreshsavedir};
    for dtc = 1:numel(dirs_to_create)
        if ~(exist(dirs_to_create{dtc},'dir') == 7)
            fprintf('did not find directory "%s", creating it\n',dirs_to_create{dtc});
            [success,msg,msgid] = mkdir(dirs_to_create{dtc});
            fprintf(' success: %d\n message: %s\n messageid: %s\n',success,msg,msgid);
        end
    end
    
    %% load data
    
    f = {'allu_ok.mat'
    'altbefits.mat'
    'befits.mat'
    'combinedAtlasAligned.mat'
    'rasterdat.mat'
    'extrabehavdat.mat'
    'keyprops.mat'
    'monkeydef.mat'
    'out_rew_ul.mat'
    'revoutplotconds.mat'
    'rpebehavdat.mat'
    'chose_info.mat'
    'chose_off1.mat'
    'glm_analysis_periods.mat'
    'ps_analysis_periods.mat'
    'combined_a_unit_norm_param_IMTv11_DRN_zo_ze.mat'};
    
    for fi=1:length(f)
        cur_filename = fullfile(basepath,'data',f{fi});
        fprintf('loading data file %2d/%2d : %s\n',fi,length(f),cur_filename);
        load(cur_filename);
        fprintf(' unitbs %d unitys %d unitypreds %d unityresids %d\n', ...
            exist('unitbs','var'), exist('unitys','var'), ...
            exist('unitypreds','var'), exist('unityresids','var'));
    end
    % re-create needed anonymous functions for altbefits from string format
    % (since matlab had trouble saving the function handles directly)
    for i = 1:numel(altbefits)
        if ~isempty(altbefits{i}.probfun_string)
            probfun = eval(altbefits{i}.probfun_string);
        else
            probfun = [];
        end
        if ~isempty(altbefits{i}.magfun_string)
            magfun = eval(altbefits{i}.magfun_string);
        else
            magfun = [];
        end
        altbefits{i}.nonlinutilfun = eval(altbefits{i}.nonlinutilfun_string);
    end
    
    fn = fieldnames(extrabehavdat);
    for fni=1:length(fn)
        eval([fn{fni} '=extrabehavdat.' fn{fni} ';']);
    end
    
    fn = fieldnames(rpebehavdat);
    for fni=1:length(fn)
        eval([fn{fni} '=rpebehavdat.' fn{fni} ';']);
    end

end

if do_run
    %% perform analyses
    % please run these in order, otherwise you may experience issues
    
    % figure 1c; supp figure 4b
    fprintf('running script: behav_from_recording\n');
    behav_from_recording 

    % supp figure 2a; supp figure 2b; supp figure 2e; supp figure 2f; supp table 4; supp table 5; supp table 6
    fprintf('running script: extra_behavior\n');
    extra_behavior

    % calls other subscripts to make anything not listed here, see comments in script for specific figures/tables
    fprintf('running script: run_offer_analysis\n');
    run_offer_analysis

    % figure 3a, figure 7a
    fprintf('running script: make_raster\n');
    make_raster

    % supp figure 2c; supp figure 2d; supp figure 7b;
    fprintf('running script: rpe_and_behavior\n');
    rpe_and_behavior

    % supp figure 7a;
    fprintf('running script: analog_analysis\n');
    analog_analysis
end
