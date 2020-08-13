function out = getRegressors(tr_align,trlcond_align,cfg)
% out = getRegressors(tr_align,trlcond_align,cfg)
%
% cfg:
%   accWin: defualt [10 0]
%   accThresh: defualt 0.8
%   onlyLearnedTrials: 0(default), 1(acc thresh), 2 (EM algo)
%   ignoreReversals: 1/0 (default=1), ignore trials where reversals occured
%   forceEM: true/false(default)
%   outcomepad: default, [-5 -4 -3 -2 -1]
%   outcomeTypes: default, [0 6]

% Copyright 2020, Benjamin Voloh
% Distributed under a GNU GENERAL PUBLIC LICENSE

%% settings
cfg = checkfield(cfg,'accWin',[10 0]);
cfg = checkfield(cfg,'accThresh',0.8);
cfg = checkfield(cfg,'onlyLearnedTrials',1);
cfg = checkfield(cfg,'ignoreReversals',1);
cfg = checkfield(cfg,'forceEM',0);
cfg = checkfield(cfg,'outcomepad',[-5 -4 -3 -2 -1]);
cfg = checkfield(cfg,'outcomeTypes',[0 6]);


%% prepare data

% select trials
if ~isfield(trlcond_align,'is_notLostIso')
    trlcond_align.is_notLostIso = true(size(trlcond_align.trialinblock));
end
trl = find( ismember(tr_align.trialOutcome,cfg.outcomeTypes) & trlcond_align.is_notLostIso );

% block stuff
block_orig = zeros(size(tr_align.trialOutcome));
d = [ 0; find(diff(tr_align.blocksLocal)~=0) ]; % ignore first block, no reversal
if d(end)~=numel(tr_align.trialOutcome); d = [d;numel(tr_align.trialOutcome)]; end

for is=1:numel(d)-1
    block_orig(d(is)+1:d(is+1)) = is;
end

% get learning trials
trl_learning = find( ismember(tr_align.trialOutcome,[0 6]) & trlcond_align.is_notLostIso );
tmp = get_learnedTrials(tr_align,trl_learning,trl,'perf',cfg.accWin,cfg.accThresh);
learnedBlock_moveavg = tmp.learnedBlocks;
learnedTrials_moveavg = tmp.learnedTrials;

if cfg.onlyLearnedTrials==2 || cfg.forceEM
    tmp = get_learnedTrials(tr_align,trl_learning,trl,'em',cfg.accWin,cfg.accThresh);
    learnedBlock_em = tmp.learnedBlocks;
    learnedTrials_em = tmp.learnedTrials;
    BLinfo = tmp.BLinfo;
else
    learnedBlock_em = [];
    learnedTrials_em = [];
    BLinfo = [];
end
    

% should we only use learned trials in our selection?
if cfg.onlyLearnedTrials==1 % thresh
    good = ismember(trl,find(learnedTrials_moveavg));
elseif cfg.onlyLearnedTrials==2
    good = ismember(trl,find(learnedTrials_em));
else
    good = true(size(tr_align.trialOutcome));
end

%should we ignore reversal trials?
if cfg.ignoreReversals && numel(d)>1
    idx = ismember(trl,d(2:end));
    good(idx) = 0;
end

% final trial selection
trl(~good) = [];
seltrl = false(size(tr_align.trialOutcome));
seltrl(trl) = 1;

% other
outcomes = tr_align.isChoice_HighReward(seltrl)==1;
choice = tr_align.Choice_Color(seltrl)-1;
choiceLocation = tr_align.Choice_Location(seltrl)-1;
choiceShapeAction = tr_align.Choice_ShapeAction(seltrl)-1;
correctColor = tr_align.blocksLocal(seltrl)-1;
sacDir = tr_align.SaccadeDir(seltrl) - 1;

staySwitch_future = [choice(1:end-1) ~= choice(2:end); nan];
staySwitch_past = [nan; choice(2:end) ~= choice(1:end-1)];
winStaySwitch = staySwitch_future;
winStaySwitch(outcomes==0) = nan;
loseStaySwitch = staySwitch_future;
loseStaySwitch(outcomes==1) = nan;
        
perf = movemean(outcomes,cfg.accWin,0,1);
if ~isempty(learnedTrials_em); learnedTrials_em(~seltrl) = []; end
learnedTrials_moveavg(~seltrl) = [];
block = block_orig(seltrl);

c = cumsum(~seltrl);
nInterveningTrials = [0;diff(c(seltrl))];
nInterveningTrials = normalizerange(nInterveningTrials,[0 1]);

% past outcomes
outcomepad = cfg.outcomepad;

outcomes_past = nan(numel(outcomes),numel(outcomepad));
for it=1:numel(outcomes)

    ind = outcomepad + it;
    ind2 = outcomepad - outcomepad(1) + 1;
    bad = ind<1 | ind>numel(outcomes);
    ind(bad) = [];
    ind2(bad) = [];

    tmp = nan(1,numel(outcomepad));
    tmp(ind2) = outcomes(ind);

    outcomes_past(it,1:numel(tmp)) = tmp;
end

%% output
out = [];
out.outcomes = outcomes;
out.choice = choice;
out.choiceLocation = choiceLocation;
out.choiceShapeAction = choiceShapeAction;
out.choiceShapeAction = choiceShapeAction;
out.correctColor = correctColor;
out.staySwitch_future = staySwitch_future;
out.staySwitch_past = staySwitch_past;
out.outcomes_past = outcomes_past;
out.winStaySwitch = winStaySwitch;
out.loseStaySwitch = loseStaySwitch;
out.sacDir = sacDir;
out.goodTrl = ismember(tr_align.trialOutcome(seltrl),[0 6]);
out.nInterveningTrials = nInterveningTrials;

out.performance = perf;
out.BLinfo = BLinfo;
out.learnedTrials_em = learnedTrials_em;
out.learnedBlock_em = learnedBlock_em;
out.learnedTrials_moveavg = learnedTrials_moveavg;
out.learnedBlock_moveavg = learnedBlock_moveavg;
out.blockidx = [d(1:end-1)+1,d(2:end)];
out.block_orig = block_orig;
out.block = block;

out.cfg = cfg;
out.trl = trl;
out.seltrl = seltrl;
out.tr_align = tr_align;
out.trlcond_align = trlcond_align;
  