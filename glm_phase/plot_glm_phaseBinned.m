[rootpath,datapath,pypath,pyenvpath] = set_decoding_paths(0);

stspath = [datapath '/MAT/STS_ALIGN_fft_distal_toi2'];
glmpath = [datapath '/glm_rate_all_R'];

%% settings

dataType = 'beta'; %beta,r2

thisMonk = 'all'; %all, ke, ha

foi = [10 25];
toi = [0.1 0.7];

% useful
theseAreas = {{'46','8a','8'},{'ACC'},{'CD','VS'}};
theseAreasStr = {'lpfc','acc','str'};

statstr = {'outcome','rpe','integrate'}; % defined by clusters

% plotting stuff
suffix = 'bin6_rand'; % last
figname = ['figure_' suffix '_' thisMonk '_' dataType];

figdir = [glmpath '/' figname];
try, if ~exist(figdir); mkdir(figdir); end, end


% data flags
% for PPC results: encoding=0, lock=0, data=1
% for phase gain results: encoding=1, lock=1, data=1
getGLMPhase = 1;
    assertEncoding = 1;
    assertLock = 1;
    assertData = 1;
    theseData = {};

getPhaseGain = 1;
getPhaseGainFreq = 0;

savePhaseGainResults = 0; % usefull to input to other scripts that compare acorss phase gain computations
    extraSuffix = ''; %maybe if you do extra processing
    saveResultsPath = [glmpath '/phaseGainResults_' suffix extraSuffix '.mat'];


% plotting flags
masterPlotFlag = 1; % dont plot anything
saveFig = 0;

getPropIndivSign = 0;
plotSketchPhaseResults = 0;
plotPhaseGainFrequencyResolved = 0;
plotPhaseGainLockNonlock = 0;
plotPhaseGainAverage = 0;
plotEncodingPhase = 0;

plotLockingDistribution = 0;
plotInterarealPPC = 1;
plotCodeNoncode = 1;

% plotting variables
phaseGainYlim = []; % [-0.5000    1.1000];
violinPercentile = [5 95];



%% useful

if getGLMPhase
    gname = 'glm_all_bigBoot.mat';
    load([glmpath '/' gname])
    
    glm_orig = glm_all;
    
    indivPath = [glmpath '/_indiv_' suffix];
    
    cfg = [];
    cfg.assertEncoding = assertEncoding;
    cfg.assertLock = assertLock;
    cfg.assertData = assertData;
    cfg.foi = foi;
    cfg.theseData = theseData;
    glm_phase = cat_glm_phase(indivPath,glm_all,cfg);
    
    % downsample glm_all
    sel = ismember({glm_all.name},{glm_phase.name});
    glm_all(~sel) = [];
end



%% ==================================================================
% ==================================================================
%                        DATA PREPARATION
% ==================================================================
% ==================================================================
% to easily perform computations, unpack everything

% useful
narea = numel(theseAreas);
nclust = 3;
[~,npred,nphase,nfreq] = size(glm_phase(1).Bphase);
freq = glm_phase(1).freq;

% info about channels
area_lfp = cat(2,glm_phase(:).area_lfp);
lfpLabel = cat(1,glm_phase(:).lfplabel)';

% data that must be resized
% - some comes from glm_all, just in case it changes (eg if you re-run
% clustering)
sz = cellfun(@size,{glm_phase(:).area_lfp},'un',0);
szOrig = sz;
func = cellfun(@(x,s) repmat({x},s),{glm_all.funcLabel},sz,'un',0);
func = [func{:}];
isencoding = cellfun(@(x,s) repmat(x,s),{glm_phase.isencoding},sz,'un',0);
isencoding = [isencoding{:}];
area_spk = cellfun(@(x,s) repmat(x,s),{glm_phase.area},sz,'un',0);
area_spk = [area_spk{:}];

names = cellfun(@(x,s) repmat({x},s),{glm_phase(:).name},sz,'un',0);
names = [names{:}];

% PPC info
tmp = cellfun(@(x) x.ppc,{glm_phase.ppcdat},'un',0);
ppc = cat(1,tmp{:});
selfreq = ismember(glm_phase(1).ppcdat.freq,freq);
ppc(:,~selfreq) = [];
ppc = get_ppcToEffectsize(ppc);

pk_freq = cat(1,glm_phase(:).pk_freq);
tmp = cellfun(@isempty,pk_freq);
pk_freq(tmp) = {nan};

islock = cellfun(@(x) any(x>=foi(1) & x<=foi(2)),pk_freq)'; 

zthresh = 6;
mu = nanmean(ppc,1);
sd = nanstd(ppc,1);
z = (ppc-mu)./sd;

outlierPPC = any(abs(z) > zthresh,2)';


% phase-GLM outputs
Bphase = cat(1,glm_phase(:).Bphase);
meanAngle = cat(1,glm_phase(:).meanAngle);
phaseRange = cat(1,glm_phase(:).phaseRange);
phaseCentre = cat(1,glm_phase(:).phaseRangeCentre);
R2 = cat(1,glm_phase(:).R2);

BphaseRand = cat(1,glm_phase.BphaseRand);
meanAngleRand = cat(1,glm_phase.meanAngleRand);
phaseCentreRand = cat(1,glm_phase(:).phaseRangeCentreRand);
R2rand = cat(1,glm_phase(:).R2rand);


% ignore badly fit models
if 1
    iswarn = cat(1,glm_phase(:).iswarn);
    iswarn = iswarn(:,:,1:size(Bphase,3),:);
    iswarn = any(iswarn,3);
    iswarn = repmat(iswarn,1,size(Bphase,2),size(Bphase,3),1);

    if strcmp(dataType,'beta')
        thresh = 20;
        iswarn(abs(Bphase)>thresh) = 1;

        iswarn_rand = repmat(iswarn,[1 1 1 1 size(BphaseRand,5)]);
        iswarn_rand(abs(BphaseRand)>thresh) = 1;
    else
        iswarn_rand = repmat(iswarn,[1 1 1 1 size(BphaseRand,5)]);
    end

    Bphase(iswarn) = nan;
    BphaseRand(iswarn_rand) = nan;
    
    R2(iswarn(:,1,:,:)) = nan;
    R2rand(iswarn_rand(:,1,:,:,:)) = nan;
end


% get current and previous outcome
Bphase = Bphase(:,2:3,:,:);
BphaseRand = BphaseRand(:,2:3,:,:,:);

% get encoding stat
if strcmp(dataType,'beta')
    [Bstat,Bstat_orig] = get_encoding_metric(Bphase,0);
    [BstatRand,~] = get_encoding_metric(BphaseRand,0);
elseif strcmp(dataType,'r2')
    [Bstat,Bstat_orig] = get_encoding_metric(R2,1);
    [BstatRand,~] = get_encoding_metric(R2rand,1);
else
    error('what data type?')
end

% area by numebr
iarea_spk = nan(size(area_spk));
iarea_lfp = nan(size(area_lfp));
for ii=1:numel(theseAreas)
   sel1 = ismember(area_spk,theseAreas{ii}); 
   sel2 = ismember(area_lfp,theseAreas{ii}); 
   iarea_spk(sel1) = ii;
   iarea_lfp(sel2) = ii;
end

  

%% ==================================================================
% ==================================================================
%                        PHASE GAIN EXTRACTION
% ==================================================================
% ==================================================================


 % =================================================================
if getPhaseGain    
    disp('>> getting phase gain')
    
    gainAmp = [];
    prefEncodingPhase = [];
    gainAmpRand = [];
    prefEncodingPhaseRand = [];
    
    prefFiringPhase = [];
        
    for id=1:numel(names)
        dotdotdot(id,0.1,numel(names))
        p = ppc(id,:);
        ff = func{id};
        
        % select frquency with maximum sig/prominent PPC
        selfreq = ismembertol(freq,pk_freq{id},0.01);
        selfreq = selfreq & (freq >= foi(1) & freq <= foi(2));
        if ~any(selfreq) % didnt lock, find max PPC in this foi
            selfreq = freq >= foi(1) & freq <= foi(2);
        end
        selfreq = find(selfreq);
        
        [~,tmp] = max(p(selfreq));
        ifreq = selfreq(tmp);
        
        % mean firing phase at this freq
        prefFiringPhase(id,1) = meanAngle(id,:,:,ifreq);

        % phase gain
        for ic=1:3
            if ~strcmp(func{id},statstr{ic})
                continue
            end
            
            % observed
            doPlot = 0;
            x = wrapTo2Pi(squeeze(phaseCentre(id,:,:,ifreq)));
            y = squeeze(Bstat(id,ic,1:end-1,ifreq));
            
            [A1, b, r] = cosinefit2(x, y, [], doPlot);
            prefEncodingPhase(id,1) = b;
            gainAmp(id,1) = get_phaseGain(A1,nanmean(y));
            
            nrand = size(BstatRand,5);
            for ir=1:nrand
                x2 = wrapTo2Pi(squeeze(phaseCentreRand(id,:,:,ifreq,ir)));
                y2 = squeeze(BstatRand(id,ic,1:end-1,ifreq,ir));
                [A2, b, r] = cosinefit2(x2, y2, [], doPlot);
                prefEncodingPhaseRand(id,1,ir) = b;
                gainAmpRand(id,1,ir) = get_phaseGain(A2,nanmean(y2));
            end
            foo=1;
        end
    end
    
    EPFG = gainAmp - nanmedian(gainAmpRand,3);  
else
    % dont overwrite if we've already calculated
    if ~exist('gainAmp')
        gainAmp = nan;
        gainAmpRand = nan;
        EPFG = nan;
    end
end

 % =================================================================
% get phase gain per frequency
[npair,~,~,nfreq] = size(Bstat);

if getPhaseGainFreq
    disp('>> getting phase gain, freq resolved')

    for id=1:npair
        dotdotdot(id,ceil(npair*0.1),npair)
        for ic=1:3
            if ~strcmp(statstr{ic},func{id})
                continue
            end

            for ifreq=1:nfreq
                x = wrapTo2Pi(squeeze(phaseCentre(id,:,:,ifreq)));
                y = squeeze(Bstat(id,ic,1:end-1,ifreq));

                [A1, b, r] = cosinefit2(x, y, [], doPlot);
                gainAmpFreq(id,ifreq,1) = get_phaseGain(A1,nanmean(y));

                nrand = size(BstatRand,5);
                for ir=1:nrand
                    x2 = wrapTo2Pi(squeeze(phaseCentreRand(id,:,:,ifreq,ir)));
                    y2 = squeeze(BstatRand(id,ic,1:end-1,ifreq,ir));
                    [A2, b, r] = cosinefit2(x2, y2, [], doPlot);
                    gainAmpFreqRand(id,ifreq,1,ir) = get_phaseGain(A2,nanmean(y2));
                end
            end
        end
    end
    
    EPFG_freq = gainAmpFreq - nanmedian(gainAmpFreqRand,4);

else
    % dont overwrite if we've already calculated
    if ~exist('gainAmpFreq')
        gainAmpFreq = nan;
        gainAmpFreqRand = nan;
        EPFG_freq = nan;
    end
end

foo=1;


%% prepare

% this data should be ignored in all analyses
badData = isnan(EPFG);
badData = badData';

foo=1;

% save data for comparison later
if savePhaseGainResults
    dat=[];
    dat.EPFG = EPFG;
    dat.EPFG_freq = EPFG_freq;
    dat.gainAmp = gainAmp;
    dat.gainAmpRand = gainAmpRand;
    dat.func = func;
    dat.area_spk = area_spk;
    dat.area_lfp = area_lfp;
    dat.names = names;
    dat.lfpLabel = lfpLabel;
    dat.freq = freq;

    disp('saving phase gain...')
    disp(saveResultsPath)
    save(saveResultsPath,'-struct','dat')
    clear dat
end

%% ==================================================================
% ==================================================================
%                               PLOTS
% ==================================================================
% ==================================================================

if ~masterPlotFlag
    warning('not plotting...')
    return
end


%% prop individually significant units
if getPropIndivSign
    
    N = [];
    labs = {};
    for ia1=1:3
        for ia2=1:3
            selarea = ismember(area_spk,theseAreas{ia1}) & ismember(area_lfp,theseAreas{ia2});
            astr = sprintf('%s-%s',theseAreasStr{ia1},theseAreasStr{ia2});
            
            sel = ~badData & selarea;
            c = gainAmp(sel);
            cr = squeeze(gainAmpRand(sel,:,:));

            nrand = size(cr,2);
            p = sum(cr>c,2) ./ nrand;
            n = sum(p<0.05);
            
            N(ia1,ia2,1) = n; % # significant
            N(ia1,ia2,2) = numel(p)-n; % # total
            lab{ia1,ia2} = astr;
        end
    end

    % stats
    mu = N(:,:,1) ./ sum(N,3);
    ex = sum(N,3) * nanmean(mu(:));
    obs = N(:,:,1);
    
    [x,p] = x2test(obs,[],[],ex);
    
    prop = (obs ./ sum(N,3))';
    
    [xx,yy]=size(obs);
    xx=1:xx;
    yy=1:yy;
    
    ntot = sum(N(:));
    nEncode = sum(sum(N(:,:,1)));
    
    % plot
    figure
    imagesc(xx,yy,prop)
    
    s = sprintf('prop sign encoding\nX2=%.3g, p=%.3g\ntot encode=%g/%g, %g',x,p,nEncode,ntot,nEncode/ntot);
    title(s)
    xlabel('spk area')
    ylabel('lfp area')
    
    set(gca,'xtick',xx,'xticklabel',theseAreasStr,'ytick',yy,'yticklabel',theseAreasStr);
    axis square
    colorbar
    
    if saveFig
        sname = [figdir '/prop_sign_encoding.pdf'];
        save2pdf(sname,gcf)
    end
    
end


%% PLOT encoding metric as a funcvtion of freq
% Figure 4G
if plotPhaseGainFrequencyResolved
    disp('>> phase gain as a function of frequency')

       
    rstr = '';

    % =================================================================
    % plots
    avgtype = 'median';
        
    figure('name','freq resolved phase gain')
    
    %get data
    c = gainAmpFreq;
    cr = nanmedian(gainAmpFreqRand,4);
    
    d = c - cr;
    
    % cluster stats
    cfg = [];
    cfg.stat = 'signrank';
    cfg.permutetime = 0;
    [sigcluster,pval] = multcompcorr_cluster(c,cr,cfg);
    p = 1 - sigcluster;
    
    % observed
    [mu,se] = avganderror(d,avgtype,1,1,200);
    mx = max(mu+se') * 1.1;
    mup = ones(size(freq)) * mx;
    mup(p>0.05) = nan;
    
    % plot     
    shadedErrorBar(freq,mu,se,{'k-'},0);
    hold all
    plot(freq,mup,'k-','linewidth',2)

    axis square
    xlabel('freq')
    ylabel('median EPFG')
    set(gca,'ylim',[0 mx*1.1],'xscale','log','xlim',[freq(1)-0.1 freq(end)+1],'xtick',10:10:freq(end));

        
    set_bigfig(gcf,[0.2*nc,0.8])
    
    % svae
    if saveFig
        sname = [figdir '/phaseGain_freq' rstr];
        save2pdf(sname,gcf)
    end
end



%% phase gain, lock vs non-lock
% Figure 4H
if plotPhaseGainLockNonlock
    disp('>> phase gain, lock vs nonlock')
    
    
    % start figure
    figure('name','phase gain: lock vs nonlock')
    nr = 1; nc = 1;
    
    cosAmp2 = EPFG;
    
    % plot average difference
    mu = [];
    se = [];
    n = [];
    str = {'nonlock','lock'};

    thisSel = nan(size(islock));
    thisSel(islock==1 & isencoding==1) = 1;
    thisSel(islock==0 & isencoding==1) = 0;
    thisSel(badData) = nan;
    
    for ii=1:2
        tmp = EPFG(thisSel==ii-1);
        [mu(ii),se(ii)] = avganderror(tmp,'median',1,1,200);
        n(ii) = sum(thisSel==ii-1);
    end
    
    a = EPFG(~isnan(thisSel));
    g = thisSel(~isnan(thisSel));
    [pk,A,~] = kruskalwallis(a,g,'off');
    x2 = A{2,5};
    
    p = [];

    for ip=1:2
        p(ip) = signrank(EPFG(thisSel==ip-1),0);
    end
            
    % plot
    ns = 1;
    subplot(nr,nc,ns)
    
    tmpsel = ~isnan(thisSel);
    tmpd = EPFG(tmpsel);
    tmpc = thisSel(tmpsel);
    
    if 1
        % cull outlier for visualization purposes
        pp = prctile(tmpd,violinPercentile);
        bad = tmpd <= pp(1) | tmpd >= pp(2);
        tmpd(bad) = [];
        tmpc(bad) = [];
    end
    
    hv = violinplot(tmpd,tmpc,...
            'violincolor',ones(1,3)*0.6,'violinAlpha',0.1,'edgecolor',ones(1,3)*0.8,...
            'DataColor',ones(1,3)*0.9,'DataAlpha',1,'width',0.4);
        
    for ih=1:numel(hv)
        set(hv(ih).ScatterPlot,'SizeData',20)
        delete(hv(ih).MedianPlot)
        delete(hv(ih).WhiskerPlot)
        delete(hv(ih).BoxPlot)
    end
    errorbar(1:numel(mu),mu,se,'k.','markersize',20)

        
    set(gca,'xticklabel',str)
    axis square
    
    s = sprintf('phase gain, lock vs nonlock\nn=%s, KWp=%.3g, Z2=%.3g,\np=%s\nmu=%s+%s',...
        mat2str(n),pk,x2,mat2str(p,3),mat2str(mu,2),mat2str(se,2));
    title(s)
    ylabel('median phase gain')
    set(gca,'fontsize',14,'xlim',[0.5 2.5])
    
    if ~isempty(phaseGainYlim)
        set(gca,'ylim',phaseGainYlim)
    end
    
    plotcueline('y',0)

    % final
    set_bigfig(gcf,[0.3 0.6])
    
    if saveFig
       sname = [figdir '/phaseGain_lock_vs_nonlock'];
       save2pdf([sname '.pdf'],gcf)
    end
end

foo=1;

  
%% PHASE GAIN, split by area OR func
% Figure 4D-F
if plotPhaseGainAverage
    disp('>> average phase gain, split by func or area')
    
    [~,ifunc] = ismember(func,statstr);
            
    % start fig
    thisStr = {'FUNC','SPK-AREA','LFP-AREA'};
    
    figure('name','phase gain')
    nr = 1; nc = 3;
    hax = [];
    for ii=1:3
        if ii==1
            thisData = ifunc;
            xticklab = statstr;
        elseif ii==2
            thisData = iarea_spk;
            xticklab = theseAreasStr;
        elseif ii==3
            thisData = iarea_lfp;
            xticklab = theseAreasStr;
        end
        

        % get data
        mu=[]; se=[];
        mur = []; ser = []; p =[]; Z = [];
        dat = {};
        for ic=1:3
            sel = thisData==ic & ~badData;
            w = EPFG(sel);

            avgtype = 'median';
            if sum(sel)>5
                [mu(ic),se(ic)] = avganderror(w,avgtype,1,1,200);
                [p(ic),~,stat] = signrank(w,0);
                try
                    Z(ic) = stat.zval;
                catch
                    Z(ic) = stat.signedrank;
                end
                foo=1;
            else
                mu(ic) = nan;
                se(ic) = nan;
                p(ic) = 1;
                Z(ic) = nan;
            end
            dat{ic} = w;
        end
    
        cmb = combnk(1:3,2);
        pp=[];
        for ic=1:size(cmb,1)
            d1 = dat{cmb(ic,1)};
            d2 = dat{cmb(ic,2)};
            d = [d1;d2];
            g = [zeros(size(d1));ones(size(d2))];
            [tmpp,a,~] = kruskalwallis(d,g,'off');
            pp(ic,:) = [cmb(ic,:), tmpp,a{2,5}];
        end

        
        [pk,A,~] = kruskalwallis(cosAmp2,thisData,'off');
        x2 = A{2,5};
    
        ns = ii;
        subplot(nr,nc,ns)
        
        tmpd = cat(1,dat{:});
        tmpc = cellfun(@(x,y) ones(size(x))*y,dat,num2cell(1:numel(dat)),'un',0);
        tmpc = cat(1,tmpc{:});
        tmpc = statstr(tmpc);
                
        if 1
            % cull outlier for visualization purposes
            pr = prctile(tmpd,violinPercentile);
            bad = tmpd <= pr(1) | tmpd >= pr(2);
            tmpd(bad) = [];
            tmpc(bad) = [];
        end
        
        mup = ones(size(mu)) * max(tmpd)*1.05;
        mup(p>0.05) = nan;
    
        
        hv = violinplot(tmpd,tmpc,...
            'violincolor',ones(1,3)*0.6,'violinAlpha',0.1,'edgecolor',ones(1,3)*0.8,...
            'DataColor',ones(1,3)*0.9,'DataAlpha',1,'width',0.4,'GroupOrder',statstr);
        

        for ih=1:numel(hv)
            set(hv(ih).ScatterPlot,'SizeData',20)
            delete(hv(ih).MedianPlot)
            delete(hv(ih).WhiskerPlot)
            delete(hv(ih).BoxPlot)
        end
        errorbar(1:numel(mu),mu,se,'k.','markersize',20)
        hold all
        plot(1:numel(mu),mup,'k*')
        plotcueline('y',0)
        
        set(gca,'xticklabel',xticklab)
        
        s = sprintf('phase gain, split by %s\np=%s\nZ=%s\nmu+se=%s\nKWp=%.3g, X2=%.3g,\n%s',...
            thisStr{ii},mat2str(p,2),mat2str(Z,3),mat2str([mu;se]',3),pk,x2,mat2str(pp,3));
        title(s)
        ylabel('phase gain')
        axis square
        hax(ii) = gca;
    end
    setaxesparameter(hax,'ylim')
    ylim = get(hax(1),'ylim');
    ylim = ylim +[0 ylim(2)*0.1];
    setaxesparameter(hax,'ylim',ylim,'xlim',[0.5 3+0.5])

    % finish
    set_bigfig(gcf,[0.7 0.5])
    
    % save
    if saveFig
        sname = [figdir '/phaseGain'];
        save2pdf(sname)
    end
end


%% PREF PHASE vs ENCODING PHASE
% Figure 5 A,C
if plotEncodingPhase
    disp('>> encoding phases')

    nphase = size(Bphase,3)-1;
    
    nbin = 8;
    [~,ifunc] = ismember(func,statstr);
    ifunc = ifunc';
    
    thetaTick = linspace(-pi,pi,4+1);
    thetaTick(end) = [];
    
    %========================================================================
    %========================================================================
    % PLOT
    
    figure('name','encoding phases')
    nr = 2; nc = 4;
    hax = [];
    for ic=1:3
        sel = ismember(func,statstr(ic)) & ~badData;
        
        % ----------------------------------------------------------------------
        % pref firing phase
        tmp = prefFiringPhase(sel);
        mu = circ_mean(tmp);
        [lowB, upB, mu2, cST,c] = circstat_confidenceOnMeanPhase_03(tmp);

        [t,r] = rose3(tmp,nbin,1);
        [po,m] = circ_otest(tmp);
        
        ns = ic;
        subplot(nr,nc,ns)
        polarplot(t,r)
        hp = addMeanDirection(tmp,r,[],1);
        s = sprintf('pref phase, %s\nn=%g,mu=%.3g+%.3g\nn=%g,HAp=%.3g, M=%.3g',...
            statstr{ic},sum(sel),mu,cST,sum(sel),po,m);
        title(s)
        set(gca,'ThetaAxisUnits','radians','thetalim',[-pi pi],'thetatick',thetaTick)
        
        hax(ic,1) = gca;
        
        % ----------------------------------------------------------------------
        % pref (relative) encoding phase
        tmp = prefEncodingPhase(sel);
                
        [t,r] = rose3(tmp,nbin,1);
        r = r ./ sum(r);
        
        tmpr = prefEncodingPhaseRand(sel);
        
        [tr,rr] = rose3(tmpr,nbin,1);
        rr = rr ./ sum(rr);
        
        mu = circ_mean(tmp);
        [lowB, upB, mu2, cST,c] = circstat_confidenceOnMeanPhase_03(tmp);

        [po,mo] = circ_otest(tmp);
        p = circ_medtest(tmp,0);
        
        ns = ic+nc;
        subplot(nr,nc,ns)
        hp = polarplot(t,r);
        hp2 = addMeanDirection(tmp,r,[],1);


        s = sprintf('prop of MAX %s\nmu=%.3g+%.3g\nHAp=%.2g,m=%.2g\nMEDp=%.2g',...
            statstr{ic},mu,cST,po,mo,p);
        title(s)
        set(gca,'ThetaAxisUnits','radians','thetalim',[-pi pi],'thetatick',thetaTick)
        
        hax(ic,2) = gca;
    end
    setaxesparameter(hax(:,2),'rlim');
    
    % ----------------------------------------------------------------------
    % stats on pref phase
    [pa,A] = circ_wwtest(prefFiringPhase,ifunc);
    paf = A{2,5};
    
    mu = []; se = [];
    dat = {};
    for ic=1:3
        sel = ismember(func,statstr{ic});
        tmp = prefFiringPhase(sel);
        [~, ~, mu(ic), se(ic),~] = circstat_confidenceOnMeanPhase_03(tmp);
        dat{ic} = tmp;
    end

    ns = nc;
    subplot(nr,nc,ns)
    barwitherr(se,mu)
    hold all
    set(gca,'xticklabel',statstr)
    ylabel('mean angle')
    
    mx = max(mu+se)*1.01;
    p_all = [];
    cmb = combnk(1:3,2);
    STAT = cmb;
    for ii=1:size(cmb,1)
        sel = ismember(ifunc,cmb(ii,:));
        tmp = prefFiringPhase(sel);
        itmp = ifunc(sel);

        [p,A] = circ_wwtest(tmp,itmp);
        p_all(ii) = p;
        STAT(ii,3:4) = [p,A{2,5}];

    end
    
    ylim = get(gca,'ylim');
    set(gca,'ylim',ylim + [0 1])
    s = sprintf('mean angle\nWWp=%.3g, F=%.3g\n%s',pa,paf,mat2str(STAT,2));
    title(s)

     % ----------------------------------------------------------------------
    % stats on relative encoding phase
    bad = isnan(prefEncodingPhase);
    [pa,A] = circ_wwtest(prefEncodingPhase(~bad),ifunc(~bad));
    paf = A{2,5};

    mu = []; se = [];
    dat = {};
    for ic=1:3
        sel = ismember(func,statstr{ic}) & ~isnan(prefEncodingPhase)';
        tmp = prefEncodingPhase(sel);
        [~, ~, mu(ic), se(ic),~] = circstat_confidenceOnMeanPhase_03(tmp);
        dat{ic} = tmp;
    end

    ns = nc*2;
    subplot(nr,nc,ns)
    barwitherr(se,mu)
    hold all
    set(gca,'xticklabel',statstr)
    ylabel('mean angle')
    
    mx = max(mu+se)*1.01;
    p_all = [];
    cmb = combnk(1:3,2);
    STAT = cmb;
    for ii=1:size(cmb,1)
        sel = ismember(ifunc,cmb(ii,:))  & ~isnan(prefEncodingPhase);
        tmp = prefEncodingPhase(sel);
        itmp = ifunc(sel);

        [p,A] = circ_wwtest(tmp,itmp);
        p_all(ii) = p;
        STAT(ii,3:4) = [p,A{2,5}];
    end
    
    ylim = get(gca,'ylim');
    set(gca,'ylim',ylim + [0 1])
    s = sprintf('mean angle\nWWp=%.3g, F=%.3g\n%s',pa,paf,mat2str(STAT,3));
    title(s)

   
    % finish
    set_bigfig(gcf,[0.7,0.8])
    
    % save
    if saveFig
        sname = [figdir '/encodingPhases'];
        save2pdf([sname,'.pdf'],gcf)
    end
end

foo=1;


%%  pseudo-sketch of RELATIVE PHASE, FUNC
% Figure 5B,D
if plotSketchPhaseResults
    disp('>> phase results sketch')

    figure('name','phase sketch')
    nr = 1; nc = 2;
    set_bigfig(gcf,[0.6 0.6])

    for ip=1:2
        
        % plot this
        doPrefPhase = ip-1;
        if doPrefPhase
            thisAngle = prefEncodingPhase;
            angstr = 'encodingAngle';
        else
            thisAngle = prefFiringPhase;
            angstr = 'firingAngle';
        end

        cols = get_safe_colors(0,[1 2 5]);

        hp = [];
        h = [];
        mu = [];
        ci = [];
        
        % loop
        for ic=1:3
            ns = ip;
            subplot(nr,nc,ns)
            
            yoff = -ic; % plot backwards

            % plot cosine
            off = 2;
            t = -off*pi:0.01:off*pi;
            hp(1) = plot(t,cos(t)+yoff,'k-','linewidth',1);
            hold all

            % data
            sel = strcmp(func,statstr{ic}) & ~badData;

            ar = thisAngle(sel);

            nboot = 10000;
            ff = @circ_mean;
            prc = [2.5 97.5];
            [pr,mur] = phase_bootstrap_conf(nboot,ff,ar,prc);


            mu(ic,1) = mur;
            ci(ic,:) = pr;
            

            selx = t >= pr(1) & t <= pr(2);
            off = 0.1;
            yy = cos(t(selx));
            yy = [yy+off, fliplr(yy)-off] + yoff;
            xx = t(selx);
            xx = [xx,fliplr(xx)];

            htmp = patch(xx,yy,cols(ic,:));
            h(ic) = htmp;
            hold all
            tmpmu = [mur mur];
            plot(tmpmu,cos(tmpmu) + [-off off]+yoff,'color','k','linewidth',2)

            set(htmp,'edgecolor','none')

        end
        
        
         % finish
        set(gca,'tickdir','out')
        set(gca,'ylim',[yoff*1.5, 1])
        set(gca,'ytick',[])
        xtick = round(-2*pi:pi/2:2*pi,2);
        set(gca,'xtick',xtick)

        plotcueline('x',0)

        tmp = [[1:3]',mu,ci];
        s = sprintf('%s: func, mean, 95CI\n%s',angstr,mat2str(tmp,2));
        title(s)
        xlabel([angstr ' phase'])

        axis square
        legend(h,statstr,'location','southoutside')

    end
    
    % save fig
    if saveFig
        sname = [figdir '/sketch_phase_results'];
        save2pdf(sname,gcf)
    end
end



%% average PPC, coding vs noncoding
% Figure 3B
if plotCodeNoncode
    % stuff
    avgtype = 'mean';
    doLog = 1;
    
    % start figure
    figure
    nr = 1; nc = 3;
    
    
    % freq resolved
    ns = 1:2;
    subplot(nr,nc,ns)
    
    h = [];
    cols = get_safe_colors(0,1:2);
    mx = 0;
    
    % plot
    for is=1:2
        sel = isencoding==is-1 & ~outlierPPC;
        p = ppc(sel,:);
        [mu,se] = avganderror(p,avgtype,1,1,200);
        
        htmp = shadedErrorBar(freq,mu,se,{'color',cols(is,:)},0);
        h(is) = htmp.mainLine;
        hold all
        
        mx = max(mx,max(mu+se'));
    end
        
    s = sprintf('ppc between distal sites\nn=%g',numel(isencoding));
    title(s)
    axis square
    if doLog
        set(gca,'xscale','log','xlim',[freq(1)-0.1,freq(end)+1]);
    else
        set(gca,'xlim',[freq(1)-1,freq(end)+1])
    end
    set(gca,'xtick',10:10:freq(end))
    
    legend(h,{'noncoding','coding'},'location','east')
    
    
    % for max beta peak
    ppc_foi = [];

    for ii=1:numel(pk_freq)
        pp = ppc(ii,:);
        p = pk_freq{ii};

        selpk = p >= foi(2) & p <= foi(2);
        if any(selpk)
            selfreq = find(ismembertol(freq,p(selpk),0.001));
        else
            selfreq = find(freq >= foi(1) & freq <= foi(2));
        end

        [~,imx] = max(pp(selfreq));
        ifreq = selfreq(imx);
        ppc_foi(ii,1) = ppc(ii,ifreq);

        if ~(freq(ifreq) >= foi(1) && freq(ifreq) <= foi(2))
            foo=1;
        end
    end
    
    p1 = ppc_foi(isencoding==0 & ~outlierPPC,:);
    p2 = ppc_foi(isencoding==1 & ~outlierPPC,:);
    
    lab = {'noncode','code'};
    
    mu = [];
    se = [];
    for ii=1:2
        sel = isencoding==ii-1 & ~outlierPPC;
        [mu(ii),se(ii)] = avganderror(ppc_foi(sel),avgtype);
    end
    [~,p,~,stat] = ttest2(p1,p2);
    
    % plot
    ns = 3;
    subplot(nr,nc,ns)
    xx = 1:numel(mu);
    errorbar(xx,mu,se,'k.','markersize',15)
    set(gca,'xtick',xx,'xticklabel',lab,'xlim',[0.5, xx(end)+0.5])
    
    s = sprintf('PPC at %s peak\nTp=%.3g, T=%.3g',mat2str(foi),p,stat.tstat);
    title(s)
    ylabel([avgtype ' PPC'])
        
    % finish
    set_bigfig(gcf,[0.5,0.4])
    
    if saveFig
        logstr = {'','_log'};
        sname = [figdir '/ppc_codeVSnoncode' logstr{doLog+1}];
        save2pdf(sname,gcf)
    end
end

foo=1;


%% average PPC, inter-areal matrix
% Figure 3D-G
if plotInterarealPPC
    avgtype = 'mean';
    doOnlyEncoding = 1;
    
    hax = [];
    
    % get avergae beta power
    ppc_foi = [];

    for ii=1:numel(pk_freq)
        pp = ppc(ii,:);
        p = pk_freq{ii};

        selpk = p >= foi(2) & p <= foi(2);
        if any(selpk)
            selfreq = find(ismembertol(freq,p(selpk),0.001));
        else
            selfreq = find(freq >= foi(1) & freq <= foi(2));
        end

        [~,imx] = max(pp(selfreq));
        ifreq = selfreq(imx);
        ppc_foi(ii,1) = ppc(ii,ifreq);

        if ~(freq(ifreq) >= foi(1) && freq(ifreq) <= foi(2))
            foo=1;
        end
    end
    
    % all specific inter-areal pairs
    MU = [];
    SE = [];
    DAT = {};
    for ia1=1:numel(theseAreas)
        for ia2 = 1:numel(theseAreas)
            sel = ismember(area_spk,theseAreas{ia1}) &...
                ismember(area_lfp,theseAreas{ia2}) & ~outlierPPC;
            
            if doOnlyEncoding; sel = sel & isencoding==1; end
            
            tmp = ppc_foi(sel);
            
            [mu,se] = avganderror(tmp,avgtype,1,1,200);
            
            MU(ia2,ia1) = mu;
            SE(ia2,ia2) = se;
            DAT{ia2,ia1} = tmp;
        end
    end
    
    % plot
    figure
    nr = 1; nc = 4;
    
    % ========================================================
    % overall average
    ns = 1;
    subplot(nr,nc,ns)
    
    xx = 1:3;
    yy = 1:3;
    imagesc(xx,yy,MU)
    
    colorbar
    axis square
    xlabel('spike site')
    ylabel('LFP site')
    set(gca,'xtick',xx,'xticklabel',theseAreasStr,'ytick',yy,'yticklabel',theseAreasStr)
    
    hax(1) = gca;
    
    % ========================================================
    % intra-areal
    sel = iarea_spk == iarea_lfp & ~outlierPPC;
    if doOnlyEncoding; sel = sel & isencoding==1; end
    tmpp = ppc_foi(sel);
    tmpi = iarea_spk(sel);
    
    [pk,A,stats] = kruskalwallis(tmpp,tmpi,'off');
    x2 = A{2,5};
    C = multcompare(stats,'display','off');

    if strcmp(avgtype,'median')
        mu = grpstats(tmpp, tmpi, {@nanmedian});
        se = grpstats(tmpp, tmpi,{@(x) nanstd(bootstrp(200,@(x) nanmedian(x),x))});
        
        [pk,A,stats] = kruskalwallis(tmpp,tmpi,'off');
        C = multcompare(stats,'display','off');
    else
        [mu,se] = grpstats(tmpp, tmpi, {'mean','sem'});
        [pk,A,stats] = anovan(tmpp,{tmpi},'display','off');
        C = multcompare(stats,'display','off');
    end

    % plot
    ns = 2;
    subplot(nr,nc,ns)
    
    xx = 1:numel(mu);
    errorbar(xx,mu,se,'k.','markersize',10)
    
    
    pc = C(:,[1 2 6]);
    s = sprintf('INTRA-areal PPC\nKWp=%.3g, X2=%.3g\nmultcompare: %s',pk,x2,mat2str(pc,2));
    title(s)
    xlabel('spike site')
    ylabel([avgtype ' PPC effect size'])
    
    axis square
    ylim = get(gca,'ylim');
    ylim = [ylim(1) ylim(2) + diff(ylim)*0.2];
    set(gca,'xtick',xx,'xticklabel',theseAreasStr,'xlim',[0.5 3.5],'ylim',ylim)

    hax(2) = gca;

    % ========================================================
    % inter-areal
    
    sel = iarea_spk ~= iarea_lfp & ~outlierPPC;
    if doOnlyEncoding; sel = sel & isencoding==1; end
    tmpp = ppc_foi(sel);
    tmpi = iarea_spk(sel);
    

     if strcmp(avgtype,'median')
        mu = grpstats(tmpp, tmpi, {@nanmedian});
        se = grpstats(tmpp, tmpi,{@(x) nanstd(bootstrp(200,@(x) nanmedian(x),x))});
        
        [pk,A,stats] = kruskalwallis(tmpp,tmpi,'off');
        x2 = A{2,5};
        C = multcompare(stats,'display','off');
    else
        [mu,se] = grpstats(tmpp, tmpi, {'mean','sem'});
        [pk,A,stats] = anovan(tmpp,{tmpi},'display','off');
        C = multcompare(stats,'display','off');
        x2 = A{2,6};

     end
    
    % plot
    ns = 3;
    subplot(nr,nc,ns)
    
    xx = 1:numel(mu);
    errorbar(xx,mu,se,'k.','markersize',10)
    
    
    pc = C(:,[1 2 6]);
    s = sprintf('INTER-areal PPC\nKWp=%.3g, X2=%.3g\nmultcompare: %s',pk,x2,mat2str(pc,2));
    title(s)
    xlabel('spike site')
    ylabel([avgtype ' PPC effect size'])
    
    axis square
    ylim = get(gca,'ylim');
    ylim = [ylim(1) ylim(2) + diff(ylim)*0.2];
    set(gca,'xtick',xx,'xticklabel',theseAreasStr,'xlim',[0.5 3.5],'ylim',ylim)

    hax(3) = gca;

    % ========================================================
    % pairwise
    tmp = combnk(1:3,2);
    cmb = nan(size(tmp,1)*2,2);
    cmb(1:2:end,:) = tmp;
    cmb(2:2:end,:) = fliplr(tmp);
    
    [~,tmpi] = ismember([iarea_spk;iarea_lfp]',cmb,'rows');
    
    sel = tmpi~=0 & ~outlierPPC';
    if doOnlyEncoding; sel = sel & isencoding'==1; end
    tmpi(~sel) = [];
    tmpp = ppc_foi(sel);
    
    if strcmp(avgtype,'median')
        mu = grpstats(tmpp, tmpi, {@nanmedian});
        se = grpstats(tmpp, tmpi,{@(x) nanstd(bootstrp(200,@(x) nanmedian(x),x))});
    else
        [mu,se] = grpstats(tmpp, tmpi, {'mean','sem'});
    end
    
    % stats
    p = []; x2 = [];
    for ip=1:3
        id1 = 2*(ip-1) + 1;
        id2 = 2*ip;
        
        sel = ismember(tmpi,[id1 id2]);
        xx = tmpp(sel);
        g = tmpi(sel);
        
        if strcmp(avgtype,'median')
            [p(ip),A,~] = kruskalwallis(xx,g,'off');
            x2 = A{2,5};
        else
            [p(ip),A,~] = anovan(xx,{g},'display','off');
            x2(ip) = A{2,6};
        end
    end
    
    % plot
    ns = 4;
    subplot(nr,nc,ns)
    
    xx = 1:numel(mu);
    errorbar(xx,mu,se,'k.','markersize',10)

    pc = C(:,[1 2 6]);
    s = sprintf('PAIRWISE PPC\nKWp = %s, X2=%s',mat2str(p,2),mat2str(x2,2));
    title(s)
    %xlabel('spike site')
    ylabel([avgtype ' PPC effect size'])
    
    axis square
    
    ylim = get(gca,'ylim');
    ylim = [ylim(1) ylim(2) + diff(ylim)*0.2];
    xlab = cellfun(@(ii) sprintf('%s-%s',theseAreasStr{cmb(ii,1)},theseAreasStr{cmb(ii,2)}),num2cell(1:size(cmb,1)),'un',0);
    set(gca,'xtick',xx,'xticklabel',xlab,'xlim',[0.5 size(cmb,1)+0.5],'ylim',ylim)

    hax(4) = gca;
    
    % ========================================================
    % finish
    set_bigfig(gcf,[1 0.4])
    
    if saveFig
        if doOnlyEncoding; encodestr ='_onlyEncoding';
        else, encodestr = '_all';
        end
        
        sname = [figdir '/ppc_mean_areaConfusionMatrix_' avgtype encodestr];
        save2pdf(sname,gcf)
    end
end



%% distribution of all freq peaks
if plotLockingDistribution
    foi_plot = [freq(1) freq(end)];
    doLog = 1;
    
    theseAreasStr2 = [theseAreasStr, {'all'}];
    
    figure;
    nr = 1; nc = 4;
    
    % loop
    for ia=1:3+1
        sel = ~outlierPPC;
        areas = cellfun(@(x,y) repmat({x},size(y)),area_spk(sel)',pk_freq(sel),'un',0);
        areas = [areas{:}];
        if ia<=3 
            selarea = ismember(areas,theseAreas{ia});
            selarea2 = ismember(area_spk(sel),theseAreas{ia});
        else
            selarea = true(size(areas));
            selarea2 = true(size(area_spk(sel)));
        end
        
    
        selfreq = freq >= foi_plot(1) & freq <= foi_plot(2);
        pk = pk_freq(sel);
        pk = [pk{:}];
        
        bad = isnan(pk) | ~selarea;
        pk(bad) = [];


        if doLog
            bins = logspace(log10(foi_plot(1)), log10(foi_plot(end)),15);
        else
            bins = 1:2:foi_plot(end);
        end

        ns = ia;
        subplot(nr,nc,ns)
        hb = histogram(pk,bins,'Normalization','probability');
        set(hb,'facecolor',ones(1,3),'facealpha',1,'EdgeColor',ones(1,3)*0.5)
        hold all


        prop = numel(pk)/sum(selarea2)*100;
        
        xi = hb.BinEdges;
        xi = nanmean([xi(1:end-1); xi(2:end)]);
        p = hb.Values;
        p = smooth(xi,p,4,'lowess');

        plot(xi,p,'r','linewidth',2)
        xlabel('freq bin')
        ylabel('# neurons w sig peak')
        str = sprintf('prop spk-lfp pairs w/ sig peak\nn=%g/%g,prop=%.3g\n%s',numel(pk),sum(selarea2),prop,theseAreasStr2{ia});

        title(str)
        axis square
        set(gca,'fontsize',12,'xlim',[bins(1),bins(end)],'xscale','log')
        grid on
        
    end
    
    % finish
    set_bigfig(gcf,[0.7 0.3])
    if saveFig
        logstr = {'','_log'};
        sname = [figdir '/ppc_peakDistirbution' logstr{doLog+1}];
        save2pdf(sname,gcf)
    end
end


