function workerIndices = splitjobs(type,ind,njob,varargin)
% workerIndices = splitjobs(type,ind,njob)
% workerIndices = splitjobs(type,ind,njob,w)
% 
% workerIndices: njobx1 cell array of indices, splitting indices "ind" by
%                   method "type"
%
% type: 'none' - jobs split as they appear
%         'spt' - shortest processsing time (specify "w")
%         'lpt' - longet processing time (specify "w")
%         - additionally, can concatante 'random' to "type", and jobs will be randomized first
% ind: indices to split
% njob: number of jobs/workers that will be used
% w: weight of each job (expected processing, memory load..., default=uniform weight)

% init
workerIndices = cell(1,njob);

%mayeb we want to randomize
if contains(type,'random')
    ind = ind(randperm(numel(ind)));
    type = erase(type,'random');
end

% nothing special, split among X workers
if strcmpi(type,'') || strcmpi(type,'none')
    for in=1:njob
        tmp = worker_indices(in,njob,numel(ind));
        workerIndices{in} = ind(tmp);
    end
end

% sorted...
if strcmpi(type,'spt') || strcmpi(type,'lpt')
    %check
    if numel(varargin)==0
        error('Specify weight vector')
    end
    w = varargin{1}; 
    if numel(w)~=numel(ind)
        error('size of weight and index vectors must match')
    end
    
    % sort
    if strcmpi(type,'spt')
        [~,isort] = sort(w,'ascend');
    elseif strcmpi(type,'lpt')
        [~,isort] = sort(w,'descend');
    end
    
    % distribute
    idx = 0;
    while numel(isort)>0
        idx = mod(idx,njob) + 1;
        workerIndices{idx} = [workerIndices{idx},ind(isort(1))];
        isort(1) = [];
    end
end


% ======================================================================
function workerIndices = worker_indices(p,P,N)
% ind=get_index(p,P,N)
% Eg: p=current worker ID, P=# of workers, N=# of datasets

st = ceil( (p-1)./P * N ) + 1;
fn = ceil( p/P * N );

workerIndices = st:fn;