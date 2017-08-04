function [Delta_ini,Sigma_ini,P_ini,ll_rec]=MA_EM_ini(X,DF,paramstruct)
% This further simplify the EM algorithm by NOT calculating full length P
% in each step. Instead, we only care about P(0,...,0) and P(1,..,1) in the
% beginning. This will serve as a good initial estimation for cases where
% the number of tissues is large
%
% Input
%
%       X        n*K matrix, Fisher transformed data
%
%       DF       K*1 vector, degree of freedom for each tissue
%
%
%   paramstruct - a Matlab structure of input parameters
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%       input as: struct('field1',values1,...
%                            'field2',values2,...
%                            'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%
%    fields            values
%
%    limit_memory     indicator of limit_memory version of the EM algorithm
%                     0   not limit_memory version
%                         faster (2X), but takes a lot of memories, b/c it needs
%                         to store a large N*2^K posterior prob matrix (8 bytes per entry)
%                     1   (default)  limit_memory version, 
%                         slower (1X) but very memory efficient, can deal with any K
%                         only need to store N*1 vectors
%
%    display           indicator of screen display
%                      0   no display
%                      1   (default) show running time and afterward log likelihood
%                      of each iteration
%
%    maxniter           maximum number of iteration, default=300
%
%    stoprule          threshold for stopping rule, default=0.1
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
%                  interpolated naively
%
%       ll_rec    vector, recording the exact log likelihood of X of each iteration
%
%
%
% copyright $Gen Li, 2013.10.29$
% modified 8/3/2017 by Gen Li
%    add maxiter and stoprule




limit_memory=1;
display=1;
maxniter=300;
stoprule=0.1;
if nargin > 2 ;   %  then paramstruct is an argument

  if isfield(paramstruct,'limit_memory') ;    %  then change to input value
    limit_memory = getfield(paramstruct,'limit_memory') ; 
  end;
  if isfield(paramstruct,'display') ;    %  then change to input value
    display = getfield(paramstruct,'display') ; 
  end;
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
Delta=0.05*(1./sqrt(DF-1))*(1./sqrt(DF-1))';
Delta(logical(eye(K)))=1./(DF-1); % replace diag value of Delta, corr(Delta)=0.05)
Sigma=6*(1./sqrt(DF-1))*(1./sqrt(DF-1))'; % corr(Sigma)=1, good start for our model
P=[0.9;zeros(2^K-2,1);0.1];


% initial likelihood
ll_opt=-inf;
Like1= MAposterior2_3(X,Delta,Mu,Sigma,P,1); % f(x,(0000))
Like2K= MAposterior2_3(X,Delta,Mu,Sigma,P,2^K); % f(x,(1111))
likelihood=Like1+Like2K;
ll=sum(log(likelihood)); 
ll_rec=ll;


% initial display
if display==1
    if limit_memory==1
        disp(['Memory-Efficient version of EM for ',num2str(K),...
            '-Tissue Model Parameter Estimation']);
    else
        disp(['NON-Memory-Efficient version of EM for ',num2str(K),...
            '-Tissue Model Parameter Estimation']);
    end;
    disp(['iter        time        log(L)']);
    disp(['0                       ',num2str(ll)]);
end;


% stopping rule (stop when ll reach maxniter, or adjacent likelihood difference is less than 0.01)
iter=1;
maxflag=0; % indicate whether a local maximizer is found

% EM iteration
while iter<maxniter 
    tic   
    % record
%     Delta_old=Delta;
%     Sigma_old=Sigma;
%     P_old=P;
    ll_old=ll; 

	% E step (see draft)
    % calculate the conditional expectation (gamma|X) of the joint log 
    % likelihood of X and gamma, where gamma is the latent group label.
    % P is separated from Delta and Sigma
    
    % M step
    %
    if limit_memory==1   % limit_memory version
      
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
        Like1= MAposterior2_3(X,Delta,Mu,Sigma,P,1); % f(x,(0000))
        Like2K= MAposterior2_3(X,Delta,Mu,Sigma,P,2^K); % f(x,(1111))
        likelihood=Like1+Like2K;
    else
        error(['Wrong limit_memory indicator!']);
    end;
    
    
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
    Delta_ini=Delta;
    Sigma_ini=Sigma;
    P_ini=P;
end;


% normalize P
P_ini=[P_ini(1);ones(2^K-2,1)*0.001;P_ini(2^K)];
P_ini=P_ini/sum(P_ini);