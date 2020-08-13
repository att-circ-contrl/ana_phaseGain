function [glm_all,clustOut,theseAreas,theseAreasStr] = get_glm_cluster(glm_all,clustOut)
% [glm_all,clustOut,theseAreas,theseAreasStr] = get_glm_cluster(glm_all)
% [glm_all,clustOut,theseAreas,theseAreasStr] = get_glm_cluster(glm_all,clustOut)
%
% Performs clustering on beta weights and assigns a functional label
%   1: cluster significant cells for each area independently
%   2: assign label
%   3: assign label to NON-CODING cells, based on distance to coding cells
%
% if clustOut not provided, clusters first

% Copyright 2020, Benjamin Voloh
% Distributed under a GNU GENERAL PUBLIC LICENSE

%% prep
if nargin<2 || ~isstruct(clustOut)
    doCluster = 1;
else
    doCluster = 0;
end

% useful
theseAreas = {{'46','8a','8'},{'ACC'},{'CD','VS'}};
theseAreasStr = {'lpfc','acc','str'};
 
% useful
area = [glm_all.area];
selsig = [glm_all.p_all]<0.05;
    

%% CLUSTER
if doCluster
    [glm_all(:).clustID] = deal(nan);

    clustOut = [];
    for ia=1:numel(theseAreas)
        fprintf('============================================\n %s\n',theseAreasStr{ia})
        sel = ismember(area,theseAreas{ia}) & selsig;
        
        tmp = glm_all(sel);
        [tmp,tmpout] = cluster_glm(tmp);
        
        tmpout.area = theseAreas{ia};
        tmpout.areaStr = theseAreasStr(ia);
        glm_all(sel) = tmp;
        clustOut = cat(1,clustOut,tmpout);

        foo=1;
    end 
end

%% assign best label for non-coding cells
tmp = glm_all(~selsig);
tmp = nearest_funcLabel(tmp,clustOut);
glm_all(~selsig) = tmp;

%% ------------------------------------------------
% FUNC

% wrapper for cellFUnction_cluster. prepares data
function [glm_all,out] = cluster_glm(glm_all)
% [glm_all,out] = cluster_glm(glm_all)

% prep beta weights
B = cat(2,glm_all(:).B)';

Bnorm = B(:,3:end); % assumes regressors for outcome -5...0 and intercept
sel = Bnorm(:,end)<0;
Bnorm(sel,:) = -Bnorm(sel,:);

out = cellFunction_cluster(Bnorm);

for id=1:numel(glm_all)
    glm_all(id).funcLabel = out.func{id};
    glm_all(id).clustID = out.clustID(id);
end

% assignes nnearest label to non-codinng cells
function glm_all = nearest_funcLabel(glm_all,clustOut)
% [glm_all,out] = nearest_funcLabel(glm_all,clustOut)
%
% assigns nonn-codinng cells to nearest group of coding cells

B = nan(numel(glm_all),5);
for id=1:numel(glm_all)
    
    % figure out which area it is
    area = glm_all(id).area;
    
    tmp = {clustOut.area};
    iarea = cellfun(@(x) any(ismember(x,area)),tmp);
    
    % find the nearest neighbor
    b = glm_all(id).B(3:end);
    if b(end)<0 % flip sign
        b = -b;
    end
    b = b ./ max(abs(b));
    B(id,:) = b;
    
    mu = clustOut(iarea).mu;
    
    c = [];
    for ic=1:size(mu,1)
        c(ic) = getCosineSimilarity(b,mu(ic,:));
    end
    [~,imx] = max(c);

    % assign a putative label
    func = clustOut(iarea).clust2func{imx};
    glm_all(id).funcLabel = func;
end
