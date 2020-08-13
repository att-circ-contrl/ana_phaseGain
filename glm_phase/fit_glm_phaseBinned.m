function [freq,B,SE,W,phaseLim,A,C,D,R2,sampleInfo] = fit_glm_phaseBinned(sts_align,nphaseBin,X,toi,trl,randomizePhase,foi,oldPhaseInfo,spikeSelectionInfo)
%[freq,B,SE,W,phaseLim,A,C,D,R2,sampleInfo] = fit_glm_phaseBinned(sts_align,nphaseBin,X,toi,trl,randomizePhase,foi,oldPhaseInfo,spikeSelectionInfo)
%
% fits GLM on spikes binnned by phase.
%
% nphaseBin: number of phase bins (will outpit nPhaseBins+1)
% X: regressors. Must have a column of ones for intercept
% toi: time of interest to calculate spike count on each trial
% trl: trials to use
% randomizePhase: 0==none (observed), 1==permutation, 2==rand phase within trial
% foi: [minFreq maxFreq]
% oldPhaseInfo: sampleInfo output of previous run. useful if phase distribution does not change, so do not have to recalculate
% spikeSelectionInfo: structure with these fields
%         selspk: selection of spikes to use
%         minSpike: minimum number of spikes
%         minTrl: minumum number of trials

% Copyright 2020, Benjamin Voloh
% Distributed under a GNU GENERAL PUBLIC LICENSE

%turn off fitting warnigs
isWarn = 'off';
warning(isWarn,'stats:glmfit:IterationLimit')
warning(isWarn,'stats:glmfit:BadScaling')
warning(isWarn,'stats:glmfit:IllConditioned')

% checks
if nargin < 8
    oldPhaseInfo = [];
end

if nargin < 9
    spikeSelectionInfo = [];
end

% useful
channels = sts_align.lfplabel;
nch = numel(channels);
freq = sts_align.freq;

selfreq = freq >= foi(1) & freq <= foi(2);
freq(~selfreq) = [];
nfreq = numel(freq);

[nsmp,npred] = size(X);

% init
B = nan(nch,size(X,2),nphaseBin+1,nfreq);
SE = nan(nch,size(X,2),nphaseBin+1,nfreq);
W = cell(nch,1,nphaseBin+1,nfreq);
phaseLim = nan(nch,1,nphaseBin+1,nfreq);
A = nan(nch,1,1,nfreq);
C = nan(nch,1,nphaseBin+1,nfreq);
D = nan(nch,1,nphaseBin+1,nfreq);
R2 = nan(nch,1,nphaseBin+1,nfreq);
SMP = cell(nch,nfreq);

% loop over all channels, freqs
for ich=1:nch
    % build response in the foi
    for ifreq=1:nfreq
        thisFreq = sts_align.freq==freq(ifreq);
        FS = sts_align.fourierspctrm{1}; %assumes only one spike channel
        FS = squeeze(FS(:,ich,thisFreq));

        % select spikes on pref vs anti-pref phases
        tr = sts_align.trial{1};
        t = sts_align.time{1};

        selt = t >= toi(1) & t<= toi(2);
        seltr = ismember(tr,trl);
        selall = selt & seltr;

        % use pre-selected spikes
        if ~isempty(spikeSelectionInfo)
            selspk = spikeSelectionInfo.selspk;
            selall = selall & selspk{ich};
            
            utr=unique(tr(selall));
            nspk = sum(selall);
            
            if nspk < spikeSelectionInfo.minSpike || numel(utr) < spikeSelectionInfo.minTrl
                %warning('after spike selection, not enough data')
                continue
            end
        end
        
        % angle
        a = angle(FS);
        phmu = circ_mean(a(~isnan(a) & selt));
        a = a - phmu;
        a = wrapToPi(a);

        selall = selall & ~isnan(a);
        
        a(~selall) = [];
        tr(~selall) = [];
        t(~selall) = [];
        
        if randomizePhase==0 %observed
            [pl,smpInR,~,~] = get_phaseBins_equalSamples(a,nphaseBin);
        elseif randomizePhase==1 %permutation
            if isempty(oldPhaseInfo)
                [pl,smpInR,~,~] = get_phaseBins_equalSamples(a,nphaseBin);
            else
                smpInR = oldPhaseInfo.smpInR{ich,ifreq};
                pl = oldPhaseInfo.phaseLim(ich,:,:,ifreq);
            end
            
            smpInR = smpInR(randperm(numel(smpInR)));
        elseif randomizePhase==2 %rand within trial
            utr = unique(tr);
            
            for ii=1:numel(utr)
                tmp = rand * 2*pi;
                tmpsel = tr==utr(ii);
                a(tmpsel) = wrapToPi(a(tmpsel) + tmp);
            end
            
            [pl,smpInR,~,~] = get_phaseBins_equalSpikes2(a,nphaseBin);
        end

        SMP{ich,ifreq} = smpInR;
        
        % loop over good, bad phases
        b = [];
        se = [];
        w = {};
        c = [];
        d = [];
        r2 = [];
        for ip=1:nphaseBin+1
            % get response on phase
            if ip <= nphaseBin
                sel = smpInR==ip;
            else
                sel = ~isnan(a);
            end
            y = get_rate(t(sel),tr(sel),trl,toi,1);

            % fit model
            lastwarn('')
            [btmp,dev,stats] = glmfit(X,y,'poisson','Constant','off'); %#ok
            b(:,ip) = btmp;
            se(:,ip) = [stats.se];
            w{1,ip} = lastwarn;
            c(ip) = sum(sel);
            d(ip) = dev;
            
            % null poisson model
            pred = nanmean(y) * ones(size(y));
            devn = sum( 2*(y.*(log((y+(y==0))./pred))-(y-pred)) ); % from glmfit func for poisson
            tmpd = 1 - dev ./ devn;
            
            r2(ip) = tmpd;
            
            foo=1;
        end

        % store
        B(ich,:,:,ifreq) = b; % beta
        SE(ich,:,:,ifreq) = se; % standard error
        W(ich,:,:,ifreq) = w; %warnings
        phaseLim(ich,:,:,ifreq) = pl; %phase bin limits
        A(ich,:,:,ifreq) = phmu; %mean phase
        C(ich,:,:,ifreq) = c; % spike count per bin
        D(ich,:,:,ifreq) = d; % deviance
        R2(ich,:,:,ifreq) = r2; % deviance squared (R2 analogue for poisson)                
    end        
end

% final
sampleInfo = [];
sampleInfo.phaseLim = phaseLim;
sampleInfo.smpInR = SMP;