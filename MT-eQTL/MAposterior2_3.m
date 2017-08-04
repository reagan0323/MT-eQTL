function PostTerm= MAposterior2_3(X,Delta,Mu,Sigma,P,index)
% this fcn returns the a N*1 vector with each entry to be f(X,gamma) for
% the index category
% namely, 2^K of this function outputs lead to MAposterior2_2
% index=1..2^K

[~,K]=size(X); 

% precalculate some critical values
gamma=double(num2str(dec2bin(0:(2^K)-1))=='1'); % 2^K*K, ordered 0-1 series
gamma=gamma(index,:); % 1*K
mu_star=bsxfun(@times,gamma,Mu'); % 1*K, mu_star vectors
sigma_star=Delta+Sigma.*(gamma'*gamma);
sigma_hfiv=sigma_star^(-0.5);
sigma_dethfiv=(det(sigma_star))^(-0.5);
clear sigma_star;

% calculate the weight vector
w=(1/(2*pi)^(K/2))*(P(index)*sigma_dethfiv); % scalar

% calculate likelihood components
% size(X)
% size(mu_star)
Xm = bsxfun(@minus, X, mu_star);
Xm = Xm*sigma_hfiv;
Xm = Xm.^2; % reuse Xm to save memory
PostTerm=exp(-0.5*sum(Xm,2))*w;

