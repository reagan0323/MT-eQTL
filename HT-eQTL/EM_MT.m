function [Delta_opt,Sigma_opt,P_opt,ll_rec]=EM_MT(X,DF,paramstruct)
% This function does parameter estimation for MT-eQTL model based on
% the modified EM algorithm proposed in the paper. More specifically,
% parameter P will be estimated in length (i.e., all 2^K entries) in each
% iteration, while Delta and Sigma will be estimated mainly based on the
% two extreme configurations.
% 
% Input
%
%       X        n*K matrix, Fisher transformed Pearson correlation data;
%                   or Zstat (corresp DF=ones(K,1)*2) 
%
%       DF       K*1 vector, degrees of freedom for corresponding tissues
%       
%       paramstruct
%           struct('XXX',xxx,...)
%       llthres  threshold for loglikelihood difference, default=(1E-3)*ll_final;
%       ngthres  threshold for decreasing log likelihood (since it's not
%                  strict EM, slight decrease of loglik is tolerable), defalt=100;
%       maxniter max number of iteration, default=200;
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
% update: 2014.4.21: new stopping rule, stop if loglik drop too much
% update: 2014.6.2: put initial estimation as a subfunction in this one
% update: 2016.12.4: stopping criterion llthres is based on EM_ini ll_final


% check dimension
[N,K]=size(X);
if K~=length(DF)
    error('Inconsistent number of tissues in X and DF');
end;

% initial estimation
[Delta,Sigma,P,ll_final]=EM_ini(X,DF);


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



% stopping rule (stop when ll reach maxniter, or adjacent likelihood difference is less than llthres)
iter=1;
maxflag=0; % indicate whether a local maximizer is found
% customized parameters
llthres=(1E-6)*abs(ll_final);
disp(['log likelihood threshold is ',num2str(llthres)]);
ngthres=100;
maxniter=200;
if nargin > 2 ;   %  then paramstruct is an argument
  if isfield(paramstruct,'llthres') ;    
    llthres = getfield(paramstruct,'llthres') ; 
  end ;
  if isfield(paramstruct,'ngthres') ;    
    ngthres = getfield(paramstruct,'ngthres') ; 
  end ;
  if isfield(paramstruct,'maxniter') ;    
    maxniter = getfield(paramstruct,'maxniter') ; 
  end ;
end;


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
%     ll_rec=[ll_rec,ll];
    if ll>ll_opt
        Delta_opt=Delta;
        Sigma_opt=Sigma;
        P_opt=P;
        ll_opt=ll;
        save('Param_temp.mat','Delta_opt','Sigma_opt','P_opt');
    else
        maxflag=1;
    end;
    
    if abs(ll-ll_old)<= llthres || ll<ll_opt - ngthres
        disp(['EM algorithm converges after ',...
            num2str(iter),' iterations!']);
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



likelihood=zeros(N,1);
for i=1:(2^K)
    likelihood=likelihood + jointlik(X,Delta_opt,Sigma_opt,P_opt,i);
end;
ll=sum(log(likelihood)); 
disp(['Final loglik for EM is ', num2str(ll,'%-10.2f')]);

end % EM







function [Delta_ini,Sigma_ini,P_ini,ll_final]=EM_ini(X,DF)
% This function returns a reasonable initial estimation for our model. It
% uses an approximate EM algorithm where we simplify the mixture of 2^K
% Gaussians to just 2: the non-eQTL distribution and the common eQTL 
% distribution. This is based on the observation that in multi-tissue
% studies, these two configurations are dominant. See details in our draft.
% Although this would not lead to accurate estimation of parameters, the 
% approximation serves as a fast initial estimation for cases where
% the number of tissues is large.
%
% We also comment that this algorithm is very efficient in memory.
%
% Input
%
%       X        N*K matrix, Fisher transformed Pearson correlation data
%
%       DF       K*1 vector, degrees of freedom for corresponding tissues
%
%
%
% Output
%
%       Delta_ini K*K positive definite matrix, optimal value estimated by
%                 EM algorithm. The diagonal values are 1/(DF-1)
%
%       Sigma_ini K*K positive definite matrix
%
%       P_opt     2^K*1 vector, other than P(1) and P(2^K), all entries are
%                 naively interpolated at the very end
%
%       ll_final   last log likelihood, provide guidance of stopping threshold
%
% Need to call:
%       jointlik.m
%
% copyright $Gen Li, 2013.10.30$
% latest update: 2014.4.21: new stopping rule, stop if loglik drop too much





% check dimension
[N,K]=size(X);
if K~=length(DF)
    error('Inconsistent number of tissues in X and DF');
end;

% initial estimation
Delta=0.05*(1./sqrt(DF-1))*(1./sqrt(DF-1))';
Delta(logical(eye(K)))=1./(DF-1); % replace diag value of Delta, corr(Delta)=0.05)
Sigma=6*(1./sqrt(DF-1))*(1./sqrt(DF-1))'; % corr(Sigma)=1, good start for our model
P=[0.9;zeros(2^K-2,1);0.1];


% initial likelihood
ll_opt=-inf;
Like1= jointlik(X,Delta,Sigma,P,1); % f(x,(0000))
Like2K= jointlik(X,Delta,Sigma,P,2^K); % f(x,(1111))
likelihood=Like1+Like2K; % pseudo likelihood
ll=sum(log(likelihood)); 
ll_rec=ll;


% initial display
% disp(['Memory-Efficient version of EM for ',num2str(K),...
%             '-Tissue Model Parameter Estimation']);


% stopping rule (stop when ll reach maxniter, or adjacent likelihood difference is less than 10)
iter=1;
llthres=10;
maxniter=100;
maxflag=0; % indicate whether a local maximizer is found

% EM iteration
while iter<maxniter
    tic   
    ll_old=ll; 

	% E step (see draft)
    % calculate the conditional expectation (gamma|X) of the joint log 
    % likelihood of X and gamma, where gamma is the latent group label.
    % P is separated from Delta and Sigma
    
    % M step
        % P
        Pi1=Like1./likelihood;
        P(1)=mean(Pi1);
        Pi2K=Like2K./likelihood;
        P(2^K)=mean(Pi2K);
        
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
        Like1= jointlik(X,Delta,Sigma,P,1); % f(x,(0000))
        Like2K= jointlik(X,Delta,Sigma,P,2^K); % f(x,(1111))
        likelihood=Like1+Like2K; % pseudo likelihood

    
    % record log likelihood
    ll=sum(log(likelihood)); % scalar, exact log likelihood of X
    ll_rec=[ll_rec,ll];
    if ll>ll_opt
        Delta_ini=Delta;
        Sigma_ini=Sigma;
        P_ini=P;
        ll_opt=ll;
    else
        maxflag=1;
    end;
    
    % check stopping criterion
    if abs(ll-ll_old)<= llthres || ll<ll_opt - 100
        break;
    else 
        iter=iter+1;
    end;

        
    time=toc;
%     disp(sprintf('%-4.0f        %-4.2f       %-10.2f', iter-1, time, ll));

end;


% if EM doesn't converge after maxniter iterations
if iter==maxniter && maxflag==0
    disp(['EM algorithm did not meet the stopping criterion after ',...
        num2str(maxniter),' iterations!']);
    Delta_ini=Delta;
    Sigma_ini=Sigma;
    P_ini=P;
end;


% normalize P to get reasonable initial estimation without zero entries
P_ini=[P_ini(1);ones(2^K-2,1)*0.001;P_ini(2^K)];
P_ini=P_ini/sum(P_ini);
ll_final=ll;
end % EM initial
