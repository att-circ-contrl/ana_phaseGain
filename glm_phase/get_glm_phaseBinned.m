
%% paths
%{
glmpath: path to glm_all, and where output will be saved
stspath: sts_align data
%}

[rootpath,datapath,pypath,pyenvpath] = set_decoding_paths(0);

glmpath = [datapath '/_decoding/glm_rate_outcomeHistory_compare_noReversal_R_all2'];
glmpath = [glmpath '/glm_goodBadPhase_revision'];

stspath = [datapath '/MAT/STS_ALIGN_fft_distal_toi2'];

% save path
suffix = '_bin6_glm_phaseBinned_test';
savepath = [glmpath '/_indiv' suffix];
if ~exist(savepath); mkdir(savepath); end


%% settings
nparallel = 1; %12;

% load, select glm
load([glmpath '/glm_all_bigBoot.mat'])

%% loop through all cells

tic
if nparallel>1 
    % make pool
    if isempty(gcp('nocreate'))
        pp = parcluster();
        pp.NumWorkers = nparallel;
        parpool('local',nparallel); 
    end

    % run
    spmd
        get_glm_phaseBinned_wrapper(glm_all,savepath,stspath)
    end
else
    get_glm_phaseBinned_wrapper(glm_all,savepath,stspath)
end

fprintf('TOTAL TIME: %g',toc)

%% main func
function get_glm_phaseBinned_wrapper(glm_all,savepath,stspath)


% settings
doOverwrite = 1; %overwrite if exits

randType = 1; % randomizePhase: 0==no rand, 1==permutation, 2==within trial
    nrand = 50;

nphaseBin = 6;

foi_ana = [2 60];
toi = [0.1 0.7];
toi_ppc = [0.1 1]; % this is like "_distal_toi"

% already done?
idx = 1:numel(glm_all);

if ~doOverwrite
    tmp = {glm_all(:).name};
    tmp = cellfun(@(x) [savepath '/' x(1:end-7) 'regPhase.mat'],tmp,'un',0);
    
    flag = cellfun(@exist,tmp) == 0;
    idx = find(flag);
    
    warning('skipping %g datasets...', numel(glm_all)-numel(idx))
end

% order of computation
w = cellfun(@numel,{glm_all.trl});
w = w(idx);
workerIndices = splitjobs('lpt',idx,numlabs,w);
IDX = workerIndices{labindex};


iid=0;
for id=IDX
    iid=iid+1;
    tic
    name = glm_all(id).name; % <dataset>_rewStart_sts.mat
    fprintf('%g/%g: %g: %s\n',iid,numel(IDX),id,name)

    load([stspath '/' name])

    % ----------------------------------------------------------------------
    % BUILD REGRESSORS
    
    % select trials
    trl = glm_all(id).trl;
    ti = glm_all(id).taskInfo;
    seltrl = ismember(ti.trl,trl);

    % regressors
    X = [ti.outcomes_past(:,end),ti.outcomes];
    X2 = [ones(size(X,1),1), X];
    X2(~seltrl,:) = [];

    cnames = {'int','out_1','out_0'};
   
    % ----------------------------------------------------------------------
    % fit glm to phase binned spikes

    % useful if you want to define specific spikes to use eg burst spikes
    spikeSelectionInfo = [];
    %spikeSelectionInfo.selspk = true(size(sts_align.time{1}));
    %spikeSelectionInfo.minTrl = 10;
    %spikeSelectionInfo.minSpike = 30;
    
    % get observed
    [freq,B,~,W,R,A,C,D,R2,sampleInfo] = fit_glm_phaseBinned(sts_align,nphaseBin,X2,toi,trl,0,foi_ana,[],spikeSelectionInfo);

    % get randomized
    if randType~=0
        Brand = nan([size(B),nrand]);
        Crand = nan([size(C),nrand]);
        Arand = nan([size(A),nrand]);
        Rrand = nan([size(R),nrand]);
        Drand = nan([size(D),nrand]);
        R2rand = nan([size(R2),1]);
        
        for ir=1:nrand
            [~,b,~,W,r,a,c,d,r2,~] = fit_glm_phaseBinned(sts_align,nphaseBin,X2,toi,trl,0,foi_ana,sampleInfo,spikeSelectionInfo);

            Brand(:,:,:,:,ir) = b;
            Crand(:,:,:,:,ir) = c;
            Arand(:,:,:,:,ir) = a;
            Rrand(:,:,:,:,ir) = r;
            Drand(:,:,:,:,ir) = d;
            R2rand(:,:,:,:,ir) = r2;
        end
    else
        Brand = nan([size(B),1]);
        Crand = nan([size(C),1]);
        Arand = nan([size(A),1]);
        Rrand = nan([size(R),1]);
        Drand = nan([size(D),1]);
        R2rand = nan([size(R2),1]);
    end

    % ----------------------------------------------------------------------
    % Other info
    
    % get ppc, peak locking
    pcfg = [];
    pcfg.trials = trl;
    pcfg.toi = toi_ppc;
    [~,ipk,ppcdat] = get_ppc_peaks(sts_align,pcfg);

    ppcdat = rmfield(ppcdat,{'ppcdat','raldat'});

    % shift so last phase bin is mean of first and last
    sz = ndims(R)+1;
    phaseCentre = R(:,:,1:nphaseBin,:);
    phaseCentre = cat(sz,phaseCentre,circshift(phaseCentre,-1,3));
    phaseCentre = circ_mean(phaseCentre,[],sz);

    sz = ndims(Rrand)+1;
    phaseCentreRand = Rrand(:,:,1:nphaseBin,:,:);
    phaseCentreRand = cat(sz,phaseCentreRand,circshift(phaseCentreRand,-1,3));
    phaseCentreRand = circ_mean(phaseCentreRand,[],sz);

    % some final stuff
    [~,ch,~,~] = lfp2spk(name,1);
    lab = sts_align.lfplabel;
    lab = strrep(lab,'CSC','');
    lab = strrep(lab,'LFP','');
    lab = strrep(lab,'CH','');
    lab = cellfun(@str2num,lab);
    sameChannel = lab==str2num(ch);

    ff = isnan(sts_align.fourierspctrm{1});
    nNan = sum(sum(ff,1),3) ./ (size(ff,1)*size(ff,3));
    badChannel = nNan==1;

    % ----------------------------------------------------------------------
    % store and save
    g = glm_all(id);

    out = [];
    out.name = name;
    out.trl = trl;
    out.iswarn = ~cellfun(@isempty,W);
    out.cnames = cnames;
    out.isencoding = g.p_all<0.05;
    try
        out.funcLabel = g.funcLabel;
        out.clustID = g.clustID;
    end
    
    out.Bphase = B;
    out.warn = W;
    out.phaseRange = R;
    out.meanAngle = A;
    out.deviance = D;
    out.R2 = R2;
    out.count = C;
    out.phaseRangeCentre = phaseCentre;

    out.meanAngleRand = Arand;
    out.countRand = Crand;
    out.BphaseRand = Brand;
    out.phaseRangeRand = Rrand;
    out.phaseRangeCentreRand = phaseCentreRand;
    out.devianceRand = Drand;
    out.R2rand = R2rand;

    out.lfplabel = sts_align.lfplabel;
    out.area = g.area;
    out.area_local = repmat(g.area,size(sts_align.area_lfp));
    out.area_lfp = sts_align.area_lfp;
    out.sameChannel = sameChannel;
    out.badChannel = badChannel;
    out.nNan = nNan;
    out.iso = g.iso;

    out.freq = freq;
    out.ppcdat = ppcdat;
    out.ipk = ipk;
    out.pk_freq = ppcdat.pk_freq;
  
    % save
    sname = [savepath '/' name(1:end-7) 'regPhase.mat'];
    save(sname,'out')
    
    % clean
    clear out sts_align in tmp
    clear B SE W R A C freq2 D R2 SMP R2adj Brand Crand Arand Rrand Drand R2rand R2adjrand

    toc
    
    foo=1;
end


% save
%disp('saving...')
%save(sname,'reg_phase_all','glm_all','clustOut','-v7.3')


end
