[rootpath,datapath,pypath,pyenvpath,lfppath] = set_decoding_paths(0);

savepath = [datapath '/_decoding/glm_rate_outcomeHistory_compare_noReversal_R_all2'];
savepath = '/Volumes/SSD_Q/RES_CODES3/_decoding/glm_rate_outcomeHistory_compare_noReversal_R_all2/glm_goodBadPhase_revision/';

figdir = [savepath '/_figures_glm_cluster'];
if ~exist(figdir); mkdir(figdir); end

if 0
   load([savepath '/glm_all_bigBoot.mat'])
end


%% settings
saveFig = 0;

plotMeanBetaByArea = 1;
plotClusterConfusion = 1;


%% useful stuff
cnames = glm_all(1).cnames(2:end);

theseAreas = {{'46','8','8a'},{'ACC'},{'CD','VS'}};
theseAreasStr = {'lpfc','acc','str'};

statstr = {'outcome','rpe','integrate'};

area = [glm_all(:).area];
[uarea,~,iarea] = unique(area);
n = accumarray(iarea,1);


%% stuff
selsig = [glm_all.p_all]<0.05;

%% PLOT mean beta weights, counts
if plotMeanBetaByArea
    % start fugure
    figure
    nr = 2; nc = numel(theseAreas);

    % plot beta weights
    for ia=1:numel(theseAreas)

        Bnorm = clustOut(ia).Bnorm;
        Bnorm = bsxfun(@rdivide,Bnorm,max(abs(Bnorm),[],2));
        ncell = size(Bnorm,1);

        % plot
        ns = ia;
        subplot(nr,nc,ns)
        cols = get_safe_colors(0,[1 2 5]);

        xx = -4:0;

        h=[];
        for ic=1:numel(statstr)
            selc = clustOut(ia).clustID == find(strcmp(clustOut(ia).clust2func,statstr{ic}));
            b = Bnorm(selc,:);
            [mu,se] = avganderror(b,'median',1,1,200);
            htmp = shadedErrorBar(xx,mu,se,{'-','color',cols(ic,:)},0);
            h(ic) = htmp.mainLine;

            if ismember(3,find(selc)) && ia==3
                foo=1;
            end
            hold all
        end

        set(gca,'xtick',xx,'fontsize',14,'ylim',[-1.2 1.2],'xlim',[xx(1)-0.5 xx(end)+0.5])
        plotcueline('y',0)
        xlabel('trials from current')
        ylabel('median norm Beta')

        legend(h,statstr,'location','southwest','fontsize',8)

        s = sprintf('%s,n=%g',cell2str(theseAreas{ia}),ncell);
        title(s)
        %grid on
        axis square
    end
    setaxesparameter('ylim')

    func = {glm_all.funcLabel};

    % counts
    N = [];
    for ia=1:numel(theseAreas)
        for ic=1:3
            sel = ismember(area,theseAreas{ia}) & strcmp(func,statstr{ic});

            N(ia,ic,1) = sum(sel & selsig);
            N(ia,ic,2) = sum(sel & ~selsig);

        end
    end

    P = N(:,:,1) ./ sum(N,3);

    % overall, diffferences between areas?
    [x,p] = x2test(N(:,:,1));
    fprintf('=========\n%s\nX2 test for difference in frequency:\nX2=%.3g, p=%.3g\n\n',mat2str(N(:,:,1)),x,p);

    % difference in cluster props?
    ob = mean(P) .* sum(N(:,:,1));
    ex = 1/3 .* sum(N(:,:,1));
    [x,p] = x2test(ob,[],[],ex);
    
    fprintf('=========\nX2 test for difference in frequency per cluster:\nX2=%.3g, p=%.3g\n\n',x,p);

    % mean proportions per clust?
    mu = nanmean(P);
    fprintf('====================================\n')
    for ip=1:numel(mu)
        fprintf('%s, MU(prop)=%.3g\n',statstr{ip},mu(ip));
    end
    
    % frewquenncies per monkey?
    fprintf('====================================\n')
    monks  = {'ke','ha'};
    for im=1:2
        NM = [];
        for ia=1:numel(theseAreas)
            for ic=1:3
                sel = ismember(area,theseAreas{ia}) & strcmp(func,statstr{ic}) & strncmp({glm_all.name},monks{im},2);

                NM(ia,ic,1) = sum(sel & selsig);
                NM(ia,ic,2) = sum(sel & ~selsig);

            end
        end
        
        s1 = sum(NM(:,:,1));
        s2 = sum(NM(:,:,2));
        st = s1+ s2;
        p = s1 ./ st;
        
        for ii=1:3
            fprintf('%s, %s: %g/%g, p=%.3g\n',monks{im},statstr{ii},s1(ii),st(ii),p(ii)); 
        end
        foo=1;
    end
    
    % plot porprotion of encoding vs non-encoding

    hp = [];
    for ia=1:3
        ns = ia + nc;
        subplot(nr,nc,ns)

        p = P(ia,:);
        h = pie(p,statstr);

        for ic=1:3
            ii=(ic-1)*2+1;
            set(h(ii),'facecolor',cols(ic,:))

            s = sprintf('%s,n=%g',get(h(ii+1),'string'),N(ia,ic,1));
            set(h(ii+1),'string',s)
        end

        % stats
        n = N(ia,:,1);
        ne = ones(size(n))*sum(n)/3;
        [x,p] = x2test(n,0,[],ne);

        s = sprintf('X=%.3g, X2p=%.3g',x,p);
        title(s)
    end

    % finish
    set_bigfig(gcf,[0.5,0.6])

    % save
    if saveFig
        sname = [figdir '/cluster_proportions.pdf'];
        save2pdf(sname,gcf)
    end
end

foo=1;
%% cluster confusion matrices
if plotClusterConfusion
        
    figure
    nr = 1; nc = 3;
    for ia=1:3
        ds = 1-clustOut(ia).DS;
        xx = 1:size(ds,1);
        yy = xx;
        
        % sort by function, 
        [~,itmp] = ismember(clustOut(ia).func,statstr);
        [~,is] = sort(itmp);
        ds2 = ds(:,is);
        ds2 = ds2(is,:);
        
        clustID = clustOut(ia).clustID(is);
        bord = find(diff(clustID)~=0)+1;

        % plot
        ns = ia;
        subplot(nr,nc,ns)
        
        imagesc(xx,yy,ds2)
        
        axis square
        
        s = sprintf('%s\nn=%g\n%s',clustOut(ia).areaStr{1},size(ds,1),cell2str(clustOut(ia).clust2func));
        title(s)
        xlabel('sorted neuron ID')
        plotcueline('x',bord,'k-','linewidth',2)
        plotcueline('y',bord,'k-','linewidth',2)
        hc = colorbar;
        set(hc.Title,'string','similarity')
    end
    
    % finish
    set_bigfig(gcf,[0.9 0.4])
    
    % save
    if saveFig
        sname = [figdir '/cluster_similarity.pdf'];
        save2pdf(sname,gcf)
    end
end

