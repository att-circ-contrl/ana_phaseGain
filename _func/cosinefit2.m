function [A, b, r] = cosinefit2(x, y, b, flag)

%[A,b,r] = cosinefit(x,y,b,flag);
%
%x is the independent variable i.e. the phase
%y is the dependent variable
%
%A is the amplitude
%b is the phase
%r is the proportion of explained variance
%
%if the third input argument is empty, the optimal phase
%will be estimated, otherwise the phase will be clamped at
%b
%if the fourth input argument is 1, feedback will be 
%given about the fit of the regression

% Copyright 2020, Thilo Womelsdorf and Benjamin Voloh
% Distributed under a GNU GENERAL PUBLIC LICENSE

A = NaN;  r = NaN; 
if isempty(b), b = NaN; end % BV: switched, now this comes first
if sum(isnan(y)>0), return, end


x = reshape(x,1,length(x));
y = reshape(y,1,length(y));

y  = y - mean(y);

n  = length(x);
S  = sum(y.*sin(x));
C  = sum(y.*cos(x));
dS = sum(sin(2.*x));
dC = sum(cos(2.*x));

if isnan(b)%isempty(b),
  nom   = S*n + S*dC - C*dS; 
  denom = C*n - S*dS - C*dC;
  if denom > 0
    b = atan(nom/denom);
  elseif denom < 0
    b = atan(nom/denom) + pi;
  elseif denom == 0
    b = 0.5 * pi;
  end
end

A = 2 * (C*cos(b)+S*sin(b)) / (dS*cos(2*b)+dC*sin(2*b)+n);

r = 1 - sum( (y - A*cos(x - b)).^2) / sum(y.^2); %proportion of explained variance

if flag == 1
  figure;
  plot(x, y); hold on;
  plot(x, A*cos(x - b), 'r');
end
