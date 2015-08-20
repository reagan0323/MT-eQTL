function [Delta_opt,Sigma_opt,P_opt,ll_rec]=EM_ref(X,DF,Delta_ini,Sigma_ini, P_ini)
% This function refines the parameter estimation for MT-eQTL model based on
% the modified EM algorithm proposed in the paper. More specifically,
% parameter P will be estimated in length (i.e., all 2^K entries) in each
% iteration, while Delta and Sigma will be estimated mainly based on the
% two extreme configurations.
%
% Input
%
%       X        n*K matrix, Fisher transformed Pearson correlation data
%
%       DF       K*1 vector, degrees of freedom for corresponding tissues
%       
%       Delta_ini, Sigma_ini, P_ini are initial estimates, potentially from
%       EM_ini.m
%
%
%
%
% Output
%
%       Delta_opt K*K positive definite matrix, optimal value estimated by
%                 EM algorithm. The diagonal values are 1/(DF-1)
%
%       Sigma_opt K*K positive definite matrix
%
%       P_opt     2^K*1 vector, optimal mass from EM algorithm
%
%       ll_rec    vector, recording the log likelihood of X in each iteration
%
%
% Need to call:
%       jointlik.m
%
% copyright $Gen Li, 2013.10.30$
%



% check dimension
[N,K]=size(X);
if K~=length(DF)
    error('Inconsistent number of tissues in X and DF');
end;

% initial estimation
Delta=Delta_ini;
Sigma=Sigma_ini; % corr(Sigma)=1, good start for our model
P=P_ini;

% initial likelihood
ll_opt=-inf;
likelihood=zeros(N,1);
for i=1:(2^K)
    likelihood=likelihood+ jointlik(X,Delta,Sigma,P,i); 
end;
ll=sum(log(likelihood)); 
ll_rec=ll;


% initial display
disp(['Memory-Efficient version of EM for ',num2str(K),...
            '-Tissue Model Parameter Estimation']);
disp(['iter        time        log(L)']);
disp(['0                       ',num2str(ll)]);



% stopping rule (stop when ll reach maxniter, or adjacent likelihood difference is less than 0.01)
iter=1;
maxniter=1000;
maxflag=0; % indicate whether a local maximizer is found

% EM iteration
while iter<maxniter
    tic   
    ll_old=ll; 

	% E step (see draft)
    % calculate the conditional expectation (gamma|X) of the joint log 
    % likelihood of X and gamma, where gamma is the latent group label.
    % P is separated from Delta and Sigma
    
    % M step (redundant calculation)
        
        % P: calculate the Post Prob, only store first column and last column
        tempP=P;
        Pi1=jointlik(X,Delta,Sigma,P,1)./likelihood;
        tempP(1)=mean(Pi1);
        Pi2K=jointlik(X,Delta,Sigma,P,2^K)./likelihood;
        tempP(2^K)=mean(Pi2K);
        for i=2:(2^K-1)
            tempP(i)=mean(jointlik(X,Delta,Sigma,P,i)./likelihood); 
        end;
        P=tempP; % only update P afterwards to maintain sum to 1
        
        % Delta
        temp_X=bsxfun(@times,X,sqrt(Pi1)); %(0,0,...)
        Delta=(temp_X'*temp_X)/sum(Pi1);
        clear temp_X;
        Delta=Delta.* ( sqrt((1./(DF-1))./diag(Delta)) * sqrt((1./(DF-1))./diag(Delta))' ); % scale Delta as a whole, maintain corrlation

        % Sigma
        temp_X=bsxfun(@times,X,sqrt(Pi2K));  %(1,1,...)
        Sigma=(temp_X'*temp_X)/sum(Pi2K)-Delta;
        clear temp_X;
        [V,D]=eig(Sigma);
        Sigma=V*max(D,0)*V'; % psd
        
    % calculate new likelihood for the new parameters
    likelihood=zeros(N,1);
    for i=1:(2^K)
        likelihood=likelihood+ jointlik(X,Delta,Sigma,P,i); 
    end; 
    ll=sum(log(likelihood)); % scalar, exact log likelihood of X
    ll_rec=[ll_rec,ll];
    if ll>ll_opt
        Delta_opt=Delta;
        Sigma_opt=Sigma;
        P_opt=P;
        ll_opt=ll;
    else
        maxflag=1;
    end;
    
    if ll<=ll_old + 0.01 && ll>=ll_old-0.01% stopping rule
        break;
    else 
        iter=iter+1;
    end;

    time=toc;
    disp(sprintf('%-4.0f        %-4.2f       %-10.2f', iter-1, time, ll));

end;

% if EM doesn't converge after maxniter iterations
if iter==maxniter && maxflag==0
    disp(['EM algorithm did not meet the stopping criterion after ',...
        num2str(maxniter),' iterations!']);
    Delta_opt=Delta;
    Sigma_opt=Sigma;
    P_opt=P;
end;

