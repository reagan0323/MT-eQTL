function lik= jointlik(X,Delta,Sigma,P,index)
% This function returns a N*1 vector of the joint likelihood of each
% observed data and one corresponding configuration specified by index,
% i.e. f(x,gamma)=p(gamma)f(x|gamma)
% Summing up the index from 1 to 2^K will lead to the full likelihood of
% the observed data X.
% Note: for simplicity, we always assume mu=0 in our model.
%
% Input: 
%           X       N*K matrix, Fisher transformed Pearson correlation
%           Delta   K*K positive definite matrix
%           Sigma   K*K positive definite matrix
%           P       2^K*1 vector, containing prior probability for each configuration
%           index   scalar, ranging from 1 to 2^K, where 1 corresponds to
%                   (0,...,0), 2 corresponds to (0,...,1), ... , 2^K
%                   corresponds to (1,...,1)
%
% Output:
%           lik     N*1 vector, containing the joint likelihood for the
%                   corresponding observation and the configuration
%
% copyright $Gen Li, 2013.10.30$

% check dimension
[~,K]=size(X); 

% precalculate some critical values
gamma=double(num2str(dec2bin(0:(2^K)-1))=='1'); % 2^K*K, ordered 0-1 series
gamma=gamma(index,:); % 1*K
sigma_star=Delta+Sigma.*(gamma'*gamma);
sigma_hfiv=sigma_star^(-0.5);
sigma_dethfiv=(det(sigma_star))^(-0.5);
clear sigma_star;

% calculate the weight vector
w=(1/(2*pi)^(K/2))*(P(index)*sigma_dethfiv); % scalar

% calculate likelihood components
X = X*sigma_hfiv;
X = X.^2; % reuse X to save memory
lik=exp(-0.5*sum(X,2))*w;
