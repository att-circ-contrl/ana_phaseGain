function glm_phase = cat_glm_phase(indivPath,glm_all,cfg)
% glm_phase = cat_glm_phase(indivPath,glm_all,cfg)
%
% selection of phase-binned GLM data. Concatenate to big sturcture
%
% indivPath: path to phase-binned GLM results for each cell-LFP pair
% glm_all: rate-based GLM structure
% cfg: configuration with fields
%       assertEncoding: only use rate-encoding cells (default=0)
%       assertLock: only use phasae-locking cell-LFP paairs  (default=0)
%       assertData: only use data that meets minimum criteria (default=0)
%       monk: 2 letter monnkey identifier (default='all')
%       foi: if assertLock=1, [min max] frequency to determine locking
%       theseData: skip all selection and retrieve this data. Nx2 cell
%       array of {dataset name, lfp channel}

% Copyright 2020, Benjamin Voloh
% Distributed under a GNU GENERAL PUBLIC LICENSE


% settings
cfg = checkfield(cfg,'assertEncoding',0);
cfg = checkfield(cfg,'assertLock',0);
cfg = checkfield(cfg,'assertData',0);
cfg = checkfield(cfg,'monk','all');
cfg = checkfield(cfg,'foi',nan);
cfg = checkfield(cfg,'theseData',[]); % Nx2 cell array, names and lfp labels

theseData = cfg.theseData;

% main
glm_phase = [];
for id=1:numel(glm_all)
    name = glm_all(id).name;
    name = [name(1:end-7) 'regPhase.mat'];
    fprintf('%g: %s\n',id,name)

    if ~strcmp(cfg.monk,'all') && ~strncmp(name,cfg.monk,2)
        continue
    end
        
    if isempty(theseData)
        thisFile = [indivPath '/' name];
        
        if exist(thisFile)
            load(thisFile)
        else
            continue
        end

        % assert stuff about the cells
        if cfg.assertEncoding && ~out.isencoding
            warning('not encoding, skipping...')
            continue
        end

        if cfg.assertData
            ti = glm_all(id).taskInfo;
            nblock = sum(ti.learnedBlock_moveavg);
            flag = nblock<2;

            if flag
                warning('indufficient data, skipping...')
                continue
            end
        end
    
        % locking?
        if cfg.assertLock
            foi = cfg.foi;
            islock = cellfun(@(x) any(x>=foi(1) & x<=foi(2)),out.pk_freq)'; 
        else
            islock = true(size(out.badChannel));
        end

        % channel selection
        bad = out.sameChannel' | out.badChannel | ~islock;

        if all(bad)
            warning('no channels left after selection...')
            continue
        end
    else
        if ~ismember({name},theseData(:,1))
            warning('not one of the requested data...')
            continue
        end
        
        % select the channels
        seldat = strcmp(theseData(:,1),name);
        ch = theseData(seldat,2);
        
        bad = ~ismember(out.lfplabel,ch);
    end

    % downsample
    out.Bphase(bad,:,:,:,:) = [];
    out.warn(bad,:,:,:,:) = [];
    out.iswarn(bad,:,:,:,:) = [];
    out.meanAngle(bad,:,:,:,:) = [];
    out.deviance(bad,:,:,:,:) = [];
    out.R2(bad,:,:,:,:) = [];
    out.phaseRangeCentre(bad,:,:,:,:) = [];
    out.meanAngleRand(bad,:,:,:,:) = [];
    
    out.countRand(bad,:,:,:,:) = [];
    out.BphaseRand(bad,:,:,:,:) = [];
    out.phaseRangeRand(bad,:,:,:,:) = [];
    out.phaseRangeCentreRand(bad,:,:,:,:) = [];
    out.devianceRand(bad,:,:,:,:) = [];
    out.R2rand(bad,:,:,:,:) = [];

    out.count(bad,:,:,:,:) = [];
    out.phaseRange(bad,:,:,:,:) = [];
    
    out.lfplabel(bad) = [];
    out.area_local(bad) = [];
    out.area_lfp(bad) = [];
    out.badChannel(bad) = [];
    out.sameChannel(bad) = [];
    out.ipk(bad) = [];
    out.pk_freq(bad) = [];
    
    out.ppcdat.ipk_orig(bad) = [];
    out.ppcdat.pk_orig(bad) = [];
    out.ppcdat.ipk(bad) = [];
    out.ppcdat.pk(bad) = [];
    out.ppcdat.pk_freq(bad) = [];
    out.ppcdat.prominence_orig(bad) = [];
    out.ppcdat.prominence(bad) = [];
    out.ppcdat.ppc(bad,:) = [];
    out.ppcdat.ral(bad,:) = [];
   
    % store
    glm_phase = cat(1,glm_phase,out);
end