function [phbin,smpInR,A_out,mu] = get_phaseBins_equalSamples(A,nbin)
% [R,smpInR,a_out,mu] = get_phaseBins_equalSamples(a)
%
% Splits angles in A into "nbin" phase bins. Phase bin limits adjusted such
% that approx equal number of smaples falls in each bin

% prep
bad = isnan(A);
A_out = A;
mu = 0;

%useulf
step = 0.001; % decreasing this provides better resolution, but is slower
nsmp = sum(~bad);

% loop until all samples are accounted for
thresh = ceil(nsmp ./ nbin); % target numebr of samples per bin
phbin = [];
smpInR = nan(size(A_out));
atmp = A_out;
for ib=1:nbin
    % first bin is centred on zero, and bin limits are symmetric
    % For each subsequent bin, extend limit to the right
    if ib==1
        lim = [-step step];
    else
        lim = [-pi -pi];
    end
    
    flag = 1;
    iter = 0;
    while flag
        iter=iter+1;
        if ib==1
            lim = lim + [-step step];
        else
            lim(2) = lim(2) + step;
        end

        % find all samples that fall in the new bin
        sel = atmp>=lim(1) & atmp<=lim(2);
        smpInR(sel) = ib;
        s = sum(sel);

        % stop if we've accounted for enough samples, or weve extended the
        % limit too far
        if s >= thresh || sum(~isnan(smpInR))==nsmp || s==nsmp || lim(2)>2*pi
            flag = 0;
        end
    end
    
    % update limits
    if ib==1
    	phbin = lim;
    else
        phbin(end+1) = phbin(end)+diff(lim);
    end
    
    % push all samples to the left of the lim
    atmp = wrapToPi(atmp - lim(2)-pi);

end

% fix the phase bin
phbin(end) = phbin(1);
phbin = wrapToPi(phbin);

