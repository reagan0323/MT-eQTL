%%%%%%%%%%%%%%%%%%%%%%%%%
% Main File for MT-eQTL %
%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the main file for MT-eQTL.
% We assume users have two data sets:
% 1. cormat
% 2. df
% The matrix "cormat" of size N*K, where each entry is a gene-SNP 
% Pearson correlation. Each row corresponds to a cis gene-SNP pair 
% (N pairs in total), and each column corresponds to a tissue 
% (K tissues in total). The donors for different tissues may vary. 
% The date set "df" is a K*1 vector recording degrees of freedom in all 
% K tissues.

%% Fisher Transformation
Fcormat=0.5*log((1+cormat)./(1-cormat));
[N,K]=size(Fcormat);

%% Fit MT-eQTL Model

% specify the model you want to fit
tissue=1:K; % e.g., assuming you want to fit a full K-dimensional model

% extract the data
DF=df(tissue);
X=Fcormat(:,tissue);

% apply modified EM algorithm
[Delta_ini,Sigma_ini,P_ini,~]=EM_ini(X,DF);
[Delta_opt,Sigma_opt,P_opt,~]=...
    EM_ref(X,DF,Delta_ini,Sigma_ini, P_ini); % XXX_opt is the optimal parameter estimation

% save the optimal parameters
save('P.mat','P_opt');
save('Sigma.mat','Sigma_opt');
save('Delta.mat','Delta_opt');


%% Discover eQTL (~00000) Using Local False Discovery Rate
% assuming you have obtained parameters for a K-dimensional model
% set lfdr threshold
FDRthres=0.05;

% calculate likelihood
likelihood=zeros(N,1);
for i=1:(2^K)
    likelihood=likelihood+ jointlik(X,Delta_opt,Sigma_opt,P_opt,i); 
	% i is index for config in binary order, e.g. K=4, i=1: 0000, i=2: 0001, i=3: 0010, i=4: 0011, ...  
end;
% testing H0: 0000
lik0=jointlik(X,Delta_opt,Sigma_opt,P_opt,1); % likelihood of f(x,0000)
lfdr=lik0./likelihood;
% for testing other H0, one may just add up several terms in lik0. For example, in K=2 case, for testing H1: (0,1), we have lik0=jointlik(X,Delta_opt,Sigma_opt,P_opt,1)+jointlik(X,Delta_opt,Sigma_opt,P_opt,3)+jointlik(X,Delta_opt,Sigma_opt,P_opt,4)


% adaptive thresholding
[fdr,~]=sort(lfdr); % fdr ascending
num_discv=find(cumsum(fdr)./[1:length(fdr)]' > FDRthres, 1, 'first')-1;
lfdrthres=fdr(num_discv);
clear fdr;
disp(['The number of eQTL in at least one tissues is ',num2str(num_discv), ' at FDR=',num2str(FDRthres),'.']);
disp(['The local false discovery rate threshold is ',num2str(lfdrthres),'.']);
eQTLind= (lfdr<=lfdrthres); % N*1 indicator vector: 1=eQTL in at least one tissue, 0=non-eQTL


% configuration assignment for a gene-SNP pair
lambda=1; % choose a gene-SNP pair, from 1 to N
if (eQTLind(lambda)==0) 
    warning('This gene-SNP pair has no eQTL across tissues under FDR=0.05! Please choose another one.');
end;
gamma=double(num2str(dec2bin(1:(2^K)-1))=='1'); % (2^K-1)*K eQTL configuration matrix
PostProb=zeros(1,2^K);
for i=1:(2^K)
    PostProb(i)=jointlik(X(lambda,:),Delta_opt,Sigma_opt,P_opt,i); 
end;
PostProb=PostProb/likelihood(lambda);
[MaxPostProb,Index]=max(PostProb(2:end),[],2);
disp('The eQTL configuration for the ',num2str(lambda),'th gene-SNP pair is ', num2str(gamma(Index,:)),'.');



% configuration assignment for all eQTL gene-SNP pairs 
% (we do NOT recommend to do this for K>4)
gamma=double(num2str(dec2bin(1:(2^K)-1))=='1'); % (2^K-1)*K eQTL configuration matrix
PostProb=zeros(sum(eQTLind),2^K);
for i=1:(2^K)
    PostProb(:,i)=jointlik(X(eQTLind,:),Delta_opt,Sigma_opt,P_opt,i)./likelihood(eQTLind); 
end;
[MaxPostProb,Index]=max(PostProb(:,2:end),[],2); 
Config=zeros(N,1); % this vector records eQTL configurations for all gene-SNP pairs.
Config(eQTLind)=Index; % 0=(0000), i=gamma(i,:) where i from 1 to 2^K-1
tabulate(Config) % 0=non-eQTL, 1=(0,..,0,1)eQTLs, 2=(0,..,1,0)eQTLs,..., 2^K-1=(1,..,1)eQTLs 


%% Scatter Plot of 2-Tissue eQTL Discoveries (an example)
% Assuming you have estimated the parameters for a 2-tissue model.
% Namely, K=2, X is N*2, and DF is 2*1.
% The parameter estimations:
% Delta_opt is 2*2, Sigma_opt is 2*2, P_opt is (2^2)*1 

% specify FDR threshold
FDRthres=0.05;


% calculate lfdr
likelihood=zeros(N,1);
for i=1:(2^K)
    likelihood=likelihood+ jointlik(X,Delta_opt,Sigma_opt,P_opt,i); 
end;
lik0=jointlik(X,Delta_opt,Sigma_opt,P_opt,1); % likelihood of f(x,0000)
lfdr=lik0./likelihood;


% adaptive thresholding
[fdr,~]=sort(lfdr); % fdr ascending
num_discv=find(cumsum(fdr)./[1:length(fdr)]' > FDRthres, 1, 'first')-1;
lfdrthres=fdr(num_discv);
clear fdr;
disp(['The number of eQTL in at least one tissues is ',num2str(num_discv), ' at FDR=',num2str(FDRthres),'.']);
disp(['The local false discovery rate threshold is ',num2str(lfdrthres),'.']);
eQTLind= (lfdr<=lfdrthres); % N*1 indicator vector: 1=eQTL in at least one tissue, 0=non-eQTL
  
% configuration assignment
PostProb=zeros(sum(eQTLind),2^K);
for i=1:(2^K)
    PostProb(:,i)=jointlik(X(eQTLind,:),Delta_opt,Sigma_opt,P_opt,i)./likelihood(eQTLind); 
end;
[MaxPostProb,Index]=max(PostProb(:,2:end),[],2); 
Config=zeros(N,1); % this vector records eQTL configurations for all gene-SNP pairs.
Config(eQTLind)=Index; % 0=(0..0), i=gamma(i,:) where i from 1 to 2^K-1
colormap= double([Index==1,Index==2,Index==3]); % sharp RGB


% plot eQTL discoveries in different configurations
lim=20;
figure(1);clf;
plot(0,0,'r.');
hold on;
plot(0,0,'g.');
plot(0,0,'b.');
plot(0,0,'w.');
scatter(sqrt(DF(1)-1)*X(eQTLind,1),sqrt(DF(2)-1)*X(eQTLind,2),1,colormap,'filled');
plot([-lim,lim],[0,0],'k-');
plot([0,0],[-lim,lim],'k-');
axis([-lim lim -lim lim]);
xlabel(['z-stat of Tissue 1'],'fontsize',35);
ylabel(['z-stat of Tissue 2'],'fontsize',35);
title(['Gene-SNP Level Scatter Plot, FDR=',num2str(FDRthres)],'fontsize',35);
set(gca,'fontsize',25);
% Plot information
ncis01=sum(Index==1);
ncis10=sum(Index==2);
ncis11=sum(Index==3);
h_legend=legend(['(0,1): ',num2str(ncis01),' Gene-SNP pairs'],... 
                ['(1,0): ',num2str(ncis10),' Gene-SNP pairs'],...
                ['(1,1): ',num2str(ncis11),' Gene-SNP pairs'],...
                    2);
set(h_legend,'FontSize',15);
