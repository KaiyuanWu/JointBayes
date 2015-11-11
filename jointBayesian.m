function [ Smu, Sepsilon ] = jointBayesian(X, labels)
%JOINTBAYESIAN Summary of this function goes here
%   Detailed explanation goes here
%   Assume X = [X1, X2, X3, ..., Xn]
%     labels = [y1, y2, y3, ..., yn]
tic
nsamples = length(labels);
nfeatures = size(X,1);

% Make sure labels are nice
[classes, ~, labels] = unique(labels);
nc = length(classes);
% Rearrange the data according to their labels
X0 = zeros(size(X,1), size(X,2));
labels0 = zeros(size(labels,1), size(labels,2));

% Current Index
nn = 1;
for c = 1:nc
    index = labels == c;
    X0(:,nn:nn + sum(index)-1) =  X(:, labels == c);
    labels0(nn:nn + sum(index)-1) =  c*ones(1, sum(labels == c));
    nn = nn + sum(index);
end
X = X0;
labels = labels0;
clear X0 labels0

% Count number of instances of each class
m =  diff([0 find(diff(sort(labels))) length(labels)]);
% Make sure numbers of instances of each class are nice
munique = unique(m);
% fprintf('Number of different counts (per person):\t%d\n', length(munique));
% for i=1:length(munique)
%     fprintf('%d ', munique(i));
% end
% fprintf('\n');

G = containers.Map({0},{[]});

% Mean of each class
mu = zeros(nfeatures, nc);
% Epsilon of each samples
epsilon = zeros(nfeatures, nsamples);
% Sum of each class
sumX = zeros(nfeatures, nc);

% Init mu and epsilon
for c = 1:nc
   mu(:, c) = mean(X(:,labels == c),2);
   sumX(:, c) = sum(X(:, labels == c), 2);
   epsilon(:, labels == c) = bsxfun(@minus, X(:,labels == c),mu(:, c) ); 
end

Smu = cov(mu');
Sepsilon = cov(epsilon');

oldSepsilon = Sepsilon;

% Maximum number iterations
max_iter = 50;
toc
fprintf('start iteration\n');
for iter = 1:max_iter
   %E Step
   F = inv(Sepsilon);
   temp = Smu*F;
   for im = 1:length(munique)
       G(munique(im)) = -(munique(im)*Smu + Sepsilon)\temp;
   end
   
   temp = Sepsilon*F;
   for c = 1:nc
      mu(:,c) = Smu*(F + m(c)*G(m(c)))*sumX(:, c);
%       fprintf('size(Sepsilon) = (%d, %d), size(F) = (%d, %d)\n', size(Sepsilon,1), size(Sepsilon,2), size(F,1), size(F,2));
%       fprintf('size(X) = (%d, %d), size(sumX) = (%d, %d)\n', size(X(:, labels == c), 1), size(X(:, labels == c),2), size(sumX(:, c),1), size(sumX(:, c),2));
      
      epsilon(:, labels == c) = temp*X(:, labels == c) + Sepsilon*G(m(c))*repmat(sumX(:, c),[1, sum(labels == c)]);
%       toc
%       fprintf('E step for class %d\n', c);
   end
   
   Smu = mu*mu'/nc;
   Sepsilon = epsilon*epsilon'/nsamples;
   
   fprintf('%d %f\n',iter,norm(Sepsilon - oldSepsilon)/norm(Sepsilon));
   toc;
   if norm(Sepsilon - oldSepsilon)/norm(Sepsilon)<1e-6
       break;
   end;
   oldSepsilon = Sepsilon; 
end
end

