function Post= MAposterior2_1(X,Delta,Mu,Sigma,P)
% This fcn returns posterior probs for all 2^K classes
% When K is large, this will lead to out-of-memory problem
%
% 7/22/12, by Gen Li
%
% INPUT
%       
%       X    n*K fisher transformed sample correlation matrix
%            each row is an observation, each column is corresponding with
%            each tissue. n can be smaller than the total number of
%            cis-pairs.
%
%       Delta     K*K covariance matrix (param) for conditional distribution 
%                 of fisher transformed sample correlations, the diagnal values
%                 should be 1/(DF-1), outside check
%
%       Mu   K*1 general mean vector (param)
%            The mixture Gaussian means are derived from this Mu
%       
%       Sigma   K*K general variance vector (param)
%               The mixture Gaussian variances are derived from this Sigma
%
%       P    (2^K)*1 mass vector (param)
%            It sums up to 1 and determines the weight of mixture Gaussian.
%            The corresponding 0-1 series are in order of binary number,
%            i.e. P(1)~(0,0,...,0)=dec2bin(0), P(2)~(0,0,...,1)=dec2bin(1),
%            P(2^K)~(1,1,...,1)=dec2bin(2^K-1)
%       
%
% OUTPUT
%
%       Post    n*(2^K) posterior probability matrix
%                   Each row corresponds to a sample, each col is for a
%                   class. Each row adds up to 1.
%

% NOTE: almost the same asa MAlikelihood2, except the last step

[n,K]=size(X); 

% precalculate some critical values
gamma=double(num2str(dec2bin(0:(2^K)-1))=='1'); % 2^K*K, ordered 0-1 series
mu_star=bsxfun(@times,gamma,Mu');% 2^K*K, mu_star vectors
gamma=mat2cell(gamma,ones(1,2^K),K); % convert to cell type
mu_star=mat2cell(mu_star,ones(1,2^K),K); 
sigma_hfiv=cell(2^K,1); % Sigma_star^(-0.5)
sigma_dethfiv=zeros(2^K,1); % |Sigma_star|^(-0.5)
for l=1:2^K
    sigma_star=Delta+Sigma.*(gamma{l}'*gamma{l});
    sigma_hfiv{l}=sigma_star^(-0.5);
    sigma_dethfiv(l)=(det(sigma_star))^(-0.5);
end;
clear sigma_star;

% calculate the weight vector
w=(1/(2*pi)^(K/2))*(P.*sigma_dethfiv); % 2^K*1

% calculate likelihood components
PostProb=zeros(n,2^K);
for l=1:2^K
    Xm = bsxfun(@minus, X, mu_star{l});
    Xm = Xm*sigma_hfiv{l};
	Xm = Xm.^2; % reuse Xm to save memory
    PostProb(:,l)=exp(-0.5*sum(Xm,2))*w(l);% likelihood
end;

% calculate posterior
Post=bsxfun(@rdivide,PostProb,sum(PostProb,2));

clear PostProb