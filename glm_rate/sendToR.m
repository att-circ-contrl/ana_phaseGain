function out = sendToR(dat,varargin)
% out = sendToR(dat)
% out = sendToR(dat,rFuncPath, rFunc)
%
% dat: (eg for glmnet wrapper) 
%   X
%   y
%   family
%   nfold
%   useMinLambda
%   verbose: default=0
%
% optional:
%       rFunc: rFunction to call
%       rFuncPath: parent diredctory of rFunc
% 
% inspired by code from Runyan et al 2017

% Copyright 2020, Benjamin Voloh
% Distributed under a GNU GENERAL PUBLIC LICENSE

% set R paths
rootpath = set_decoding_paths(0);
if contains(computer,'WIN') % PC
    rpath = 'C:\Program Files\R\R-3.6.0\bin\Rscript';
else % MAC
    rpath = '/usr/local/bin/Rscript';
end

% the R funtion to call
if numel(varargin) > 0
    rFuncPath = varargin{1};
    rFunc = varargin{2};
else
    rFuncPath = checkpath( [rootpath '/git/ana_acc_codes/_new/ana_glm_r'] );
    rFunc = 'cvglmnetFromMatlab';
end

fprintf('***  calling R  ***\npath: %s\nfunc: %s\n',rFuncPath,rFunc);

% check verbose
dat = checkfield(dat,'verbose',0);
verbose = dat.verbose;
dat = rmfield(dat,'verbose');

% prepare tmp file names
tmpname = tempname;
name_in = [tmpname '_in.mat'];
name_out = [tmpname '_out.mat'];

% save input data and call R
save(name_in,'-struct','dat')

commandStr = sprintf('"%s" %s/%s.R "%s"',rpath,rFuncPath,rFunc,tmpname);
if verbose
    [status,result] = system(commandStr,'-echo');
else
    [status,result] = system(commandStr);
end

% if succesful, load the output data and update with info
if status==0
    %load back the data
    outtmp = load(name_out);

    % add the original data
    out = outtmp;
    f = fieldnames(dat);
    for ii = 1:length(f)
        if ~isfield(out,f{ii})
            out.(f{ii}) = dat.(f{ii});
        else
            warning('field "%s" already exists, wont add original...',f{ii})
        end
    end

    doClean(name_in,name_out)
else
    out = [];
    doClean(name_in,name_out)
end

foo=1;

%% FUNCS

% clean temp files
function doClean(name_in,name_out)
try, delete(name_in), end
try, delete(name_out), end