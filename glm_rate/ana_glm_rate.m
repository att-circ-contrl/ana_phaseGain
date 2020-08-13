% wrapper function to fit GLM to outcomes


% Copyright 2020, Benjamin Voloh
% Distributed under a GNU GENERAL PUBLIC LICENSE

[rootpath,datapath] = set_ana_paths(0);

stspath = [datapath '/MAT/STS_ALIGN_fft_distal_toi2'];

%% settings
doRun = 0;
    nparallel = 1;
    overwrite = 0;
    onlySaveSignificant = 0;

% saving stuff
targetDir = 'glm_rate_all_R';
savepath = [datapath '/' targetDir];
if ~exist(savepath); mkdir(savepath); end

%init parallelization
if nparallel > 1 && isempty(gcp('nocreate'))
    disp('===============================================')
    disp('starting parallel pool...')
    pp = parcluster();
    pp.NumWorkers = nparallel;
    parpool('local',nparallel); 
    disp(pp)
end


%% loop over all cells and fit GLM
if doRun
    d = dir([stspath '/*rewStart_sts.mat']);
    parfor id=1:numel(d)
        tic 
        name = d(id).name;
        disp('================================================================')
        disp([num2str(id) ': ' name])
        
        sname = [savepath '/' name(1:end-7) 'glm.mat'];
        if ~overwrite && exist(sname,'file')
            disp('skipping...')
            continue
        end
        
        in = load([stspath '/' name])
        
        % unpack
        sts_align = in.sts_align;
        tr_align = in.tr_align;
        trlcond_align = in.trlcond_align;
        hd_align = in.hd_align;
        ev_align = in.ev_align;

        % fit glm
        fcfg = [];
        fcfg.toilim = [0.1 0.7];
        fcfg.minTrials = 25;
        fcfg.minRate = 1;
        fcfg.family = 'poisson';
        fcfg.useMinLambda = 1;
        fcfg.nfold = 10;
        fcfg.nrand = 5;
        fcfg.verbose = 1;
        
        out = fit_glm_rate_R(fcfg,sts_align,tr_align,trlcond_align);
            
        out.name = name;
        out.iso = sts_align.iso;
        out.area = sts_align.area;
        
        if isempty(out)==1; continue; end %model fit failed

        % save
        if onlySaveSignificant && out.p_all > 0.05
            warning('not significant, not saving...')
            continue
        end
        
        sname = [savepath '/' name(1:end-7) 'glm.mat'];
        S = [];
        S.out = out;
        parsave(sname,S)

        toc 
        foo=1;
    end
end

%% concatenate all output into a single structure
glmname = [savepath '/glm_all.mat'];

% concatenate other stuff
d = dir([savepath '/*glm.mat']);

glm_all = [];
for id=1:numel(d)
    dotdotdot(id,0.1,numel(d));

    name = d(id).name;

    in=load([savepath '/' name]);
    glm_all = cat(1,glm_all,in.out);

end

% resave
fprintf('saving: %s\n',glmname)
save(glmname,'glm_all')
