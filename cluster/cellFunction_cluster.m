function out = cellFunction_cluster(Bnorm)
% out = cellFunction_cluster(Bnorm)
%
% Clusters samples on the basis of beta weights. Clusters are defind in
% terms of stability (how often they are clustered together when
% bootstrapping)

% Copyright 2020, Benjamin Voloh
% Distributed under a GNU GENERAL PUBLIC LICENSE

%% useful
[ncell,nreg] = size(Bnorm);

clust = 'kmeans';
clustalgo = 'Silhouette';
dist = 'cosine';

nboot = 100000;

fprintf('functional clustering: %s, eval=%s,dist=%s\n',clust,clustalgo,dist)

%% run it many times to get the optimum cluster

if 0 % if nclust is stable, can skip this step
    nrun = 1000;

    K = nan(nrun,1);
    eva_all = cell(nrun,1);
    fprintf('evaluating optimal K')
    for ib=1:nrun
        dotdotdot(ib,ceil(nrun*0.1),nrun)
        tmp = evalclusters(Bnorm,clust,clustalgo,'klist',1:10,'Distance',dist);
        K(ib) = tmp.OptimalK;
        eva_all{ib} = tmp;
    end

    nclust = mode(K);
    eva = eva_all{find(K==nclust,1)};
else
    nclust = 3;
    eva = [];
    eva_all = {};
    K = nan;
end

%% run it many times to determine neuron cluster membership
fprintf(' evaluating cluster membership')

dat = zeros(ncell);
clustN_boot = nan(nboot,nclust);
for ib=1:nboot
    dotdotdot(ib,ceil(nboot*0.1),nboot)
    idx = randi(ncell,1,ncell);
    
    tmp = Bnorm(idx,:);
    [tmpc,~] = kmeans(tmp,nclust,'Distance',dist);
    [~,~,ic] = unique(tmpc);
    clustN_boot(ib,:) = accumarray(ic,1);
    
    for ic=1:nclust
        sel = idx(tmpc==ic);

        for is1=1:numel(sel)
            for is2=1:numel(sel)
                if sel(is1)==sel(is2)
                    continue
                end
                dat(sel(is1),sel(is2)) = dat(sel(is1),sel(is2))+1;
            end
        end

    end

    foo=1;
end
    
% final clustering based on stability
D = dat ./ max(dat(:));
D(eye(size(D))==1) = 1;
DS = 1-D; % dissimilarity
y = squareform(DS);
z = linkage(y,'complete');
    
clustID = cluster(z,'maxclust',nclust);
[~,~,ic] = unique(clustID);
clustN = accumarray(ic,1);

[idSort,is] = sort(clustID);
DSsort = DS(:,is);
DSsort = DSsort(is,:);
Dsort = D(:,is);
Dsort = Dsort(is,:);


% silhouette
% Rousseeuw 1978, J comp and Applied Math, Silhouettes: a graphical aid to the interpretation and validation of cluster analysis
% page 55-57

tmps = [];
for ic1=1:nclust
    for ic2=1:nclust
       sel1 = idSort==ic1; 
       sel2 = idSort==ic2; 
       tmp = DSsort(sel1,sel2);

       tmps(ic1,ic2) = mean(tmp(:));

       foo=1;
    end
end

a = diag(tmps);
b = tmps;
b(eye(size(b))==1)=Inf;
b = min(b,[],2);

sel = a>b;
silh = 1 - a ./ b;
silh(sel) = b(sel)./a(sel) - 1;

%% likley funtional label
func = cell(ncell,1);
clust2func = {};
if 1
    MU = [];
    SE = [];
    for ic=1:nclust
        sel = clustID==ic;
        tmp = Bnorm(sel,:);
        [mu,se] = avganderror(tmp,'mean',1);
      
        % label the cluster type
        s = sign(mu);
        if mu(end)>0 && abs(mu(end-1) ./ mu(end))<0.2
            lab = 'outcome';
        elseif s(end-1) ~= s(end)
            lab = 'rpe';
        elseif s(end) == s(end-1)
            lab = 'integrate';
        else
            lab = ['clust' num2str(ic)];
        end

        %update
        clust2func{ic} = lab;
        func(sel) = {lab};
        MU(ic,:) = mu;
        SE(ic,:) = se;

    end 
end

%% output
out = [];
out.clust = clust;
out.clustalgo = clustalgo;
out.nclust = nclust;
out.K = K;
out.eva = eva;
out.eva_all = eva_all;

out.clustID = clustID;
out.clustN = clustN;
out.clustN_boot = clustN_boot;
out.Dorig = dat;
out.D = D;
out.DS = DS;
out.DSsort = DSsort;
out.Dsort = Dsort;
out.silh = silh;

out.func = func;
out.clust2func = clust2func;
out.mu = MU;
out.se = SE;
out.Bnorm = Bnorm;

