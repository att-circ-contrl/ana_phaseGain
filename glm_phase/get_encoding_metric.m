function [Bstat,Bstat_orig] = get_encoding_metric(X,outputType)
% [Bstat,Bstat_orig] = get_encoding_metric(X,outputType)
%
% Calculates the encoding metric. X can be different input (eg beta
% weights, R2). X should be Nx2x[...] matrix (can add more dimensions too). Based
% on the input, outputType defines the output
%
% outputType:
%       0: default, 3 diff beta calculations
%           [Nx3x...] matrix, in order {'outcome','rpe','rewHist'}
%       1: for data of form chan X 1 X bin X freq (eg R2)
% output is in order {'outcome','rpe','rewHist'}

if nargin<2
    outputType = 0;
end

% prep the stats
if outputType==0
    tmp1 = X(:,2,:,:,:);
    tmp2 = diff(X(:,1:2,:,:,:),[],2);
    tmp3 = sum(X(:,1:2,:,:,:),2);
elseif outputType==1
    tmp1 = X(:,1,:,:,:);
    tmp2 = X(:,1,:,:,:);
    tmp3 = X(:,1,:,:,:);
end

Bstat_orig = cat(2,tmp1,tmp2,tmp3);
Bstat = abs(Bstat_orig);
