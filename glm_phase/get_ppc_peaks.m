function varargout = get_ppc_peaks(sts_align,cfg)
% [pk,ipk] = get_ppc_peaks(sts_align,cfg)
% [pk,ipk,dat] = get_ppc_peaks(sts_align,cfg)
%
% calcualte PPC and rayleigh for each channnel in sts align, and extract
% significant, prominent PPC peaks
%
% cfg:
%   userange: defualt=0, use the range of values to get minpeakheight
%               (sets  minpeakprominence=0)
%   thresh: deault 0.25, used if userange=1
%   minpeakprominence: default=0.005
%   minpeakheight: default = 0.005;
%   minpeakdistance: default=3;
%   foi
%   toi
%   trials

% Copyright 2020, Benjamin Voloh
% Distributed under a GNU GENERAL PUBLIC LICENSE

%% checks/settings
cfg = checkfield(cfg,'trials','all');
cfg = checkfield(cfg,'userange',0);
cfg = checkfield(cfg,'thresh',0.25);
cfg = checkfield(cfg,'minpeakprominence',0.005);
cfg = checkfield(cfg,'minpeakheight',0.005);
cfg = checkfield(cfg,'minpeakdistance',3);
cfg = checkfield(cfg,'foi',[min(sts_align.freq), max(sts_align.freq)]);
cfg = checkfield(cfg,'toi',[min(sts_align.time{1}), max(sts_align.time{1})]);
cfg = checkfield(cfg,'foi_lock',nan);
cfg = checkfield(cfg,'doplot',0);
cfg = checkfield(cfg,'channel','all');

foi = cfg.foi;
toi = cfg.toi;

% freq selection
selfreq = sts_align.freq >= foi(1) & sts_align.freq <= foi(2);
thisFreq = sts_align.freq(selfreq);

%% get ppc, rayleigh stats

% PPC
scfg = [];
scfg.method = 'ppc2';
scfg.latency = toi;
scfg.foi = thisFreq;
scfg.trials = cfg.trials;
scfg.channel = cfg.channel;
ppc = ft_spiketriggeredspectrum_stat(scfg,sts_align);

% Rayeligh
scfg = [];
scfg.method = 'ral';
scfg.latency = toi;
scfg.foi = thisFreq;%
scfg.trials = cfg.trials;
scfg.channel = cfg.channel;
ray = ft_spiketriggeredspectrum_stat(scfg,sts_align);

% get freqpks
freq = ppc.freq;
pp = ppc.ppc2;
pp(pp<0) = 0;

if cfg.userange
    cfg.minpeakheight = cfg.thresh * range(pp);
    cfg.minpeakprominence = 0;
end
mph = cfg.minpeakheight;
mpd = cfg.minpeakdistance;
mpp = cfg.minpeakprominence;

% loop over alll channels
nch = size(pp,1);

pk_orig = cell(nch,1);
ipk_orig = cell(nch,1);
pk = cell(nch,1);
ipk = cell(nch,1);
pk_freq = cell(nch,1);
P_orig = cell(nch,1);
P = cell(nch,1);

for ich=1:nch
    [pk_orig{ich},ipk_orig{ich},~,P_orig{ich}] = findpeaks(pp(ich,:),'MinPeakHeight',mph,'MinPeakDistance',mpd,'MinPeakProminence',mpp);
    
    % return peaks wiyth significant phase locking
    selpk = ismember(freq,freq(ipk_orig{ich}));
    selsig = ray.ral(ich,selpk) < 0.05;
    pk{ich} = pk_orig{ich}(selsig);
    ipk{ich} = ipk_orig{ich}(selsig);
    P{ich} = P_orig{ich}(selsig);
    pk_freq{ich} = freq(ipk{ich});
end

%% output
varargout{1} = pk;
varargout{2} = ipk;
if nargout>2
    dat =[];
    dat.ipk_orig = ipk_orig;
    dat.pk_orig = pk_orig;
    dat.ipk = ipk;
    dat.pk = pk;
    dat.pk_freq = pk_freq;
    dat.ppc = ppc.ppc2;
    dat.ral = ray.ral;
    dat.freq = thisFreq;
    dat.prominence_orig = P_orig;
    dat.prominence = P;
    dat.cfg = cfg;
    dat.ppcdat = ppc;
    dat.raldat = ray;
    
    varargout{3} = dat;
end

%% plot?
if cfg.doplot
    figure
    
    subplot(1,2,1)
    plot(freq,pp)
    hold all
    plot(freq(ipk_orig),pk_orig,'go')
    plot(freq(ipk),pk,'rx')
    legend({'ppc','pk','sig-pk'},'location','eastoutside')
    xlabel('freq')
    ylabel(scfg.method)
    axis square
    
    subplot(1,2,2)
    plot(freq,ray.ral)
    plotcueline('y',0.05)
    xlabel('freq')
    ylabel('rayleigh P')
    axis square
    
    set_bigfig(gcf,[0.5 0.3])
end



