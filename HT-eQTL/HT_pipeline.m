function [Delta,Sigma,P,P_sparse,sparse_ind]=HT_pipeline(Zstat)
% This function estimates parameters based on HT-eQTL
% Input is standardized Z-stat matrix X
% Output is all model parameters
% Need to call EM_MT.m
 

[~,K]=size(Zstat);
% fit all pairwire models

for i=1:(K-1)
    for j=(i+1):K
        X=Zstat(:,[i,j]);
        DF=[2;2];
        name=[num2str(i),'_',num2str(j)];
        [Delta_opt,Sigma_opt,P_opt,~]=EM_MT(X,DF); %,struct('llthres',.1)
        save(['Param_MT_',name],'Delta_opt','Sigma_opt','P_opt');
    end;
end;


% assemble Delta and Sigma
% Delta
Delta_HT=zeros(K,K); 
for i=1:(K-1)
    for j=(i+1):K
        name=[num2str(i),'_',num2str(j)];
        load(['Param_MT_',name]);
        Delta_HT(i,j)=Delta_opt(1,2);  % diagonal should be 1
    end;
end;
Delta_HT=Delta_HT+Delta_HT'+eye(K);
[V,D]=eig(Delta_HT); % make it positive definite
[~,Delta]=cov2corr(V*abs(D)*V'); 

% Sigma (ad hoc diag)
Sigma_corr=zeros(K,K);
Sigma_diag=zeros(K,K);
for i=1:(K-1)
    for j=(i+1):K
        name=[num2str(i),'_',num2str(j)];
        load(['Param_MT_',name]);
        Sigma_corr(i,j)=Sigma_opt(1,2)./sqrt(Sigma_opt(1,1)*Sigma_opt(2,2));        
        Sigma_diag(i,j)=Sigma_opt(1,1); % tissue i'th sigma from (i,j)
        Sigma_diag(j,i)=Sigma_opt(2,2); % tissue j'th sigma from (i,j)
    end;
end;
Sigma_corr=Sigma_corr+Sigma_corr'+eye(K);
[V,D]=eig(Sigma_corr); % make it positive definite
[~,Sigma_corr]=cov2corr(V*abs(D)*V');
Sigma_diag=min(Sigma_diag+diag(inf(K,1)),[],2); % ad hoc, take minimum 2MT Sigma Diag as HT Sigma Diag
Sigma=corr2cov(sqrt(Sigma_diag),Sigma_corr);



% Reproduce P (multi probit model)
% estimate correlation structure and final thresholds 
Rho=zeros(K,K); % K-variate correlation matrix
thres=zeros(K,K); % each row for one tissue
for i=1:(K-1)
    for j=(i+1):K
        name=[num2str(i),'_',num2str(j)];
        load(['Param_MT_',name]);

        Pi=P_opt(1)+P_opt(2); % null prob of tissue i
        Pj=P_opt(1)+P_opt(3);
        thresi=norminv(Pi); % threshold for tissue i: <thresi it's null
        thresj=norminv(Pj);
        thres(i,j)=thresi; % threshold for tissue i
        thres(j,i)=thresj;
        myfun=@(rho) (mvncdf([thresi,thresj],[0,0],[1,rho;rho,1])-P_opt(1))^2;
        
        % in practice, one could first run grid search to get a good
        % initial guess and then use fminsearch to refine the est of Rho     
        xx=0.5:0.01:0.99; % 50-grid search, corr from .5 to 1
        yy=zeros(size(xx));
        for tt=1:length(xx)
            yy(tt)=myfun(xx(tt));
        end;
        [~,I]=min(yy);
        Rho(i,j)=fminsearchbnd(myfun,xx(I),0.5,1); % maybe a bit slow but should be more accurate
    end;
end;
Rho=Rho+Rho'+eye(K);
[V,D]=eig(Rho); % make it positive definite
[~,Rho]=cov2corr(V*abs(D)*V');
thres=thres+diag(inf(K,1));
finalthres=min(thres,[],2); 

% estimate pmf based on mvn integral (matlab limit dim<=25)
P_HT=zeros(2^K,1);
gamma=double(num2str(dec2bin(0:(2^K)-1))=='1'); % 2^K*K, ordered 0-1 series
parfor i=1:(2^K)
    config=gamma(i,:);
    xl=-inf(1,K);xu=inf(1,K);
    xl(logical(config))=finalthres(logical(config));
    xu(logical(1-config))=finalthres(logical(1-config));
    P_HT(i)=mvncdf(xl,xu,zeros(1,K),Rho);
end;
sum(P_HT) %0.9995
P=P_HT;

% truncate the smallest pmfs
[sP,I]=sort(P_HT,'ascend');
tau=1E-5; % threshold for every config prior
num=2^K-find(sP>tau,1)+1; % leftover
disp(['Left ', num2str(100*(num)/2^K),'% configs (', num2str(num),'/',num2str(2^K),')']);
P_sparse=P_HT;
P_sparse(I(1:(2^K-num)))=0;
P_sparse=P_sparse/sum(P_sparse); % sum to 1
sparse_ind=sort(I((2^K-num+1):end)); % index of non-zero entries in P_sparse

