function F = fastmcd(P)

% F = fastmcd(P)
%
% Sample trimming using fast-MCD

n = size(P, 1);
p = size(P, 2);

if n < 2
  warning(sprintf('Input only contains %d samples\n', n));
  F = P;
  return;
end

h = floor((n+p+1)/2) + 1;

numStarts = 100;

D_list = zeros(numStarts, 1);
hInd_list = zeros(numStarts, h);

for i=1:numStarts

  % Get h unique random sample indices
  hInd = randomInd(n, h);
  
  % Perform 3 c-step iterations
  [hInd_prime, D] = csteps(P, hInd, 3);
  
  D_list(i) = D;
  hInd_list(i, :) = hInd_prime';

end

numCandidates = 10;

[D_list, si] = sort(D_list);
candidatesInd = si(1:numCandidates);

D_list = zeros(numCandidates, 1);
hInd_list = hInd_list(candidatesInd, :);

for i=1:numCandidates
  
  % Obtain the previously estimated subset
  hInd = hInd_list(i, :);
  
  % Perform 500 c-step iterations
  [hInd_prime, D] = csteps(P, hInd, 800);
  
  D_list(i) = D;
  hInd_list(i, :) = hInd_prime';

end

% Create the distribution estimate from the best candidate subsets
[D_list, si] = sort(D_list);
hInd = hInd_list(si(1), :);

X = P(hInd, :);
[mu, sigma] = getparams(X);

% Filter outliers using the estimated distribution
% by thresholding the Mahalanobis distances
dist = (P - repmat(mu, n, 1)) / sqrtm(sigma);
dist = dist .^ 2;
dist = sum(dist, 2);
dist = sqrt(dist);
inlierInd = dist <= 3.0;
F = P(inlierInd, :);

return;

%--------------------------------------------------------------------------


function [hInd_prime, D] = csteps(P, hInd, maxIters)

h = length(hInd);
n = size(P, 1);

% Compute initial distribution parameters
X = P(hInd, :);
[mu, sigma] = getparams(X);
D = det(sigma);
Dold = D;

numIters = 0;

while 1
  
  % Compute squared Mahalanobis distances
  dist = (P - repmat(mu, n, 1)) / sqrtm(sigma);
  dist = dist .^ 2;
  dist = sum(dist, 2);
  
  [dist, si] = sort(dist);
  hInd_prime = si(1:h);
  X = P(hInd_prime, :);
  [mu, sigma] = getparams(X);
  Dold = D;
  D = det(sigma);
  
  % Check convergence
  numIters = numIters + 1;
  deltaD = abs(Dold - D);
  if numIters >= maxIters | deltaD < 1e-6
    break;
  end
  
end

return;

%--------------------------------------------------------------------------

function [mu, sigma] = getparams(X)

mu = sum(X, 1) ./ size(X, 1);
sigma = cov(X);

return;
