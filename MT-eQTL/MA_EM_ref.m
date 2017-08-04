function [Delta_opt,Sigma_opt,P_opt,ll_rec]=MA_EM_ref(X,DF,Delta_ini,Sigma_ini, P_ini, paramstruct)
% This function is similar with MA_EM2. The only difference is that we
% specify the initial value here. This is very efficient for large number
% of tissues. Usually, we use MA_EM_ini to get reasonable initial values.
%
%
% Input
%
%       X        n*K matrix, Fisher transformed data
%
%       DF       K*1 vector, degree of freedom for each tissue
%       
%       Delta_ini, Sigma_ini, P_ini are initial estimates
%
%   paramstruct - a Matlab structure of input parameters
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%    fields            values
%
%    maxniter           maximum number of iteration, default=300
%
%    stoprule          threshold for stopping rule, default=1
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
%       ll_rec    vector, recording the exact log likelihood of X of each iteration
%
%
%
% copyright $Gen Li, 2013.10.29$
% modified 8/3/2017 by Gen Li
%    add maxiter and stoprule

% need to call:
%   MAlikelihood2.m (memory-efficient)


maxniter=300;
stoprule=0.1;
if nargin > 5 ;   %  then paramstruct is an argument
  if isfield(paramstruct,'maxniter') ;    %  then change to input value
    maxniter = getfield(paramstruct,'maxniter') ; 
  end;
  if isfield(paramstruct,'stoprule') ;    %  then change to input value
    stoprule = getfield(paramstruct,'stoprule') ; 
  end;
end ;

% check
[N,K]=size(X);
if K~=length(DF)
    error('Inconsistent number of tissues in X and DF');
end;

% initial estimation
Mu=zeros(K,1);
Delta=Delta_ini;
Sigma=Sigma_ini; % corr(Sigma)=1, good start for our model
P=P_ini;

% initial likelihood
ll_opt=-inf;
likelihood=MAlikelihood2(X,Delta,Mu,Sigma,P); % N*1, exact likelihood of X
ll=sum(log(likelihood)); 
ll_rec=ll;


% initial display
disp(['Memory-Efficient version of EM for ',num2str(K),...
            '-Tissue Model Parameter Estimation']);
disp(['iter        time        log(L)']);
disp(['0                       ',num2str(ll)]);



% stopping rule (stop when ll reach maxniter, or adjacent likelihood difference is less than 0.01)
iter=1;
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
    %
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
    w=(1/(2*pi)^(K/2))*(P.*sigma_dethfiv); % 2^K*1

    % calculate the Post Prob, only store first column and last column
    Xm = bsxfun(@minus, X, mu_star{1});
    Xm = Xm*sigma_hfiv{1};
    Xm = Xm.^2; % reuse Xm to save memory
    Pi1=(exp(-0.5*sum(Xm,2))*w(1))./likelihood;
    P(1)=mean(Pi1);

    Xm = bsxfun(@minus, X, mu_star{2^K});
    Xm = Xm*sigma_hfiv{2^K};
    Xm = Xm.^2; % reuse Xm to save memory
    Pi2K=(exp(-0.5*sum(Xm,2))*w(2^K))./likelihood;
    P(2^K)=mean(Pi2K);

    for l=2:2^K-1
        Xm = bsxfun(@minus, X, mu_star{l});
        Xm = Xm*sigma_hfiv{l};
        Xm = Xm.^2; % reuse Xm to save memory
        P(l)=mean((exp(-0.5*sum(Xm,2))*w(l))./likelihood);% likelihood
    end;
    clear Xm;

    % Delta
    temp_X=bsxfun(@times,X,sqrt(Pi1)); %(0,0,...)
    Delta=(temp_X'*temp_X)/sum(Pi1);
    clear temp_X;
%         Delta(logical(eye(K)))=1./(DF-1); % replace diag value of Delta
    Delta=Delta.* ( sqrt((1./(DF-1))./diag(Delta)) * sqrt((1./(DF-1))./diag(Delta))' ); % scale Delta as a whole, maintain corrlation
    % Sigma
    temp_X=bsxfun(@times,X,sqrt(Pi2K));  %(1,1,...)
    Sigma=(temp_X'*temp_X)/sum(Pi2K)-Delta;
    clear temp_X;
    [V,D]=eig(Sigma);
    Sigma=V*max(D,0)*V'; % psd

    % calculate new likelihood for the new parameters
    likelihood=MAlikelihood2(X,Delta,Mu,Sigma,P); % N*1 new likelihood
    
    
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
    
    if ll<=ll_old + stoprule && ll>=ll_old - stoprule
        break;
    else 
        iter=iter+1;
    end;
        
       
    
    time=toc;
    disp(sprintf('%-4.0f        %-4.2f       %-10.2f', iter-1, time, ll));

end;



% if EM doesn't converge.....
if iter==maxniter && maxflag==0
    disp(['EM algorithm did not meet the stopping criterion after ',...
        num2str(maxniter),' iterations!']);
    Delta_opt=Delta;
    Sigma_opt=Sigma;
    P_opt=P;
end;

