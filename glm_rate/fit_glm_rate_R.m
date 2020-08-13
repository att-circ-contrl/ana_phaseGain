function out = fit_glm_rate_R(f_cfg,sts_align,tr_align,trlcond_align)
% out = fit_glm_rate_R(CFG,sts_align,tr_align,trlcond_align)

%% settings

f_cfg = checkfield(f_cfg,'toilim','needit');
f_cfg = checkfield(f_cfg,'minTrials',25);
f_cfg = checkfield(f_cfg,'minBlocks',2);
f_cfg = checkfield(f_cfg,'minRate',1);

f_cfg = checkfield(f_cfg,'family','poisson');
f_cfg = checkfield(f_cfg,'useMinLambda',1);
f_cfg = checkfield(f_cfg,'nfold',10);
f_cfg = checkfield(f_cfg,'nrand',1000);
f_cfg = checkfield(f_cfg,'verbose',0);


% get overall rate
[TC] = condition_selection(trlcond_align,'spk');

cfg = [];
cfg.latency = f_cfg.toilim;
cfg.trials = TC.trials_all;
cfg.keeptrials = 'yes';
cfg.trackcallinfo = 'no';
rate_overall = ft_spike_rate(cfg,sts_align);
            
out = []; % start output if we exit early

%% prepare data

gcfg = [];
gcfg.onlyLearnedTrials = 1;
gcfg.ignoreReversals = 1;
ti = getRegressors(tr_align,trlcond_align,gcfg);
trl = ti.trl;
if numel(trl) < f_cfg.minTrials; return; end

% build output
cfg = [];
cfg.latency = f_cfg.toilim;
cfg.trials = trl;
cfg.keeptrials = 'yes';
cfg.outputunit = 'spikecount'; %rate
rate = ft_spike_rate(cfg,sts_align);

y = rate.trial;

% model specification
X = [ti.outcomes_past,ti.outcomes];
idx = 1:size(ti.outcomes_past,2);
idx = [idx-idx(end)-1, 0];
cnames = cellfun(@(x) strrep(['Outcome' num2str(x)],'-','_'),num2cell(idx),'un',0);

% delete nans
bad = any(isnan(X),2);
X(bad,:) = [];
y(bad) = [];
trl(bad) = [];

seltrl = false(size(ti.seltrl));
seltrl(trl) = 1;

% checks
if rate_overall.avg < f_cfg.minRate
    warning('rate is too low')
    return
end

if numel(trl) < f_cfg.minTrials
    warning('too few trials')
    return
end

nblock = sum(ti.learnedBlock_moveavg);
if nblock < f_cfg.minBlock
    warning('too few blocks')
end

%% fit the GLM
dat = [];
dat.X = X;
dat.y = y;
dat.family = f_cfg.family;
dat.useMinLambda = f_cfg.useMinLambda;
dat.nfold = f_cfg.nfold;
dat.nrand = f_cfg.nrand;
dat.verbose = f_cfg.verbose; %1;

rFuncPath = '/Users/ben/Desktop/phase_code/_orig/_glm_rate';
rFunc = 'cvglmnetFromMatlab';
out = sendToR(dat,rFuncPath,rFunc);

%% output

% ignore if we failed
if isempty(out)
    return
end

% stats
out.p1se = 1 - sum(out.B1se~=0) ./ out.nrand;
out.pmin = 1 - sum(out.Bmin~=0) ./ out.nrand;
out.p = 1 - sum(out.Ball~=0) ./ out.nrand;
out.p_all_min = 1-sum(any(out.Bmin(:,2:end)~=0,2)) ./ out.nrand;
out.p_all_1se = 1-sum(any(out.B1se(:,2:end)~=0,2)) ./ out.nrand;
out.p_all = 1-sum(any(out.Ball(:,2:end)~=0,2)) ./ out.nrand;

% data
out.cnames = cnames;
out.taskInfo = ti;
out.trl = trl;
out.X = X;
out.y = y;
out.seltrl = seltrl;
out.warning = '';

% used for later checks
out.isBadlyConditioned = 0;
out.isencoding = any(out.p<0.05);
out.rate_overall = rate_overall.avg;
out.learnedBlocks = ti.learnedBlock_moveavg;
out.ntrl = numel(trl);

