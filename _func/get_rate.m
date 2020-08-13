function trialRate = get_rate(spkTime,spkTrial,trials,toi,returnCount)
% r = get_rate(spkTime,spkTrial,trials,toi,returnCount)
% 
% spike rate in toi. inspired by ft_spike_rate

% Copyright 2020, Benjamin Voloh
% Distributed under a GNU GENERAL PUBLIC LICENSE

selt = spkTime >= toi(1) & spkTime < toi(2);
seltrl = ismember(spkTrial,trials);

tr = spkTrial(selt & seltrl);
trialBins   = sort([trials-0.5; trials+0.5]);
trialRate   = histc(tr(:),trialBins);
trialRate   = trialRate(1:2:end-1); % the uneven bins correspond to the trial integers

%update
if ~returnCount
    trialRate = trialRate ./ diff(toi);
end

if isrow(trialRate); trialRate = trialRate'; end
