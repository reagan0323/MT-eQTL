% This file conducts complete simulation analysis for MT paper 
% 4-tissue study, with true parameters from MT on GTEx pilot data
%
% By Gen Li, 8/3/2017
clear all

%% Load parameters
df=[137;100;119;86]; % in the paper
df=df-2; % in the code
P=[0.7721,0.0202,0.0190,0.0037,0.0104,0.0033,0.0010,0.0107,0.0196,0.0010,0.0008,0.0009,0.0029,0.0085,0.0019,0.1240]';
Delta=[1,0.1347,0.0805,0.1089
    0.1347,1,0.1204,0.1794
    0.0805,0.1204,1,0.1288
    0.1089,0.1794,0.1288,1]; % in the paper
Delta=bsxfun(@rdivide,bsxfun(@rdivide,Delta,sqrt(df-1)),sqrt(df'-1)); % in the code
Sigma=[6.5699, 5.3098, 4.4683, 4.7126
    5.3098, 5.9752, 4.7906, 5.5778
    4.4683, 4.7906, 5.5263, 4.6493
    4.7126, 5.5778, 4.6493, 6.0178]; % in the paper
Sigma=bsxfun(@rdivide,bsxfun(@rdivide,Sigma,sqrt(df-1)),sqrt(df'-1)); % in the code


K=4;
%% Generate data
n=1E7;
% Number for each configuration
gamma=double(num2str(dec2bin(0:(2^K)-1))=='1');
rand('seed',1); % 2013.11.1
nsub=histc(rand(1,n),[0,cumsum(P)']);
nsub=nsub(1:(end-1)); % this is the truth in sim data

randn('seed',1); 
Fcormat=mvnrnd(zeros(1,K),Delta,nsub(1));
for i=2:2^K
    Fcormat=[Fcormat;...
        mvnrnd(zeros(1,K),Delta+diag(gamma(i,:))*Sigma*diag(gamma(i,:)),nsub(i))];
end;
DF=df;
X=Fcormat;
% True eQTL label (just for convenience)
true_gamma=[];
for i=1:2^K
    true_gamma=[true_gamma;kron(ones(nsub(i),1),gamma(i,:))];
end;


%% paramete estimation and comparison
time=tic;
[Delta_ini,Sigma_ini,P_ini,ll_rec]=MA_EM_ini(X,DF,struct('stoprule',0.01,'maxniter',1000));
[Delta_opt,Sigma_opt,P_opt,ll_rec]=MA_EM_ref(X,DF,Delta_ini,Sigma_ini, P_ini,struct('stoprule',0.01,'maxniter',1000));
out=toc(time)
Mu_opt=zeros(K,1);

% compare with true parameters
(100*abs(Delta_opt-Delta)./Delta)
(100*abs(Sigma_opt-Sigma)./Sigma)
100*abs(P_opt-P)./P



%% test different hypothesis
% calculate N*2^K posterior matrix
X=Fcormat;
Post=MAposterior2_1(X,Delta_opt,Mu_opt,Sigma_opt,P_opt);   
size(Post) % 10E7 * 16
FDRthres=0.05;

% test H0=0000
lfdr=Post(:,1);
[fdr,~]=sort(lfdr); % fdr ascending
num_discv=find(cumsum(fdr)./[1:length(fdr)]' > FDRthres, 1, 'first')-1;
fdrthres=fdr(num_discv);
truealt=(sum(true_gamma,2)~=0);
discv=(lfdr<=fdrthres)
num_alt=sum(truealt)
num_discv=sum(discv)
overlap=sum(truealt & discv)
discvratio=overlap/num_alt % TPR
obsfdr=1-overlap/num_discv % FDR



% test H0 0~~~ 
for tisind=1:4; 
    nullind=find(gamma(:,tisind)==0);
    lfdr=sum(Post(:,nullind),2);
    [fdr,~]=sort(lfdr); % fdr ascending
    num_discv=find(cumsum(fdr)./[1:length(fdr)]' > FDRthres, 1, 'first')-1;
    fdrthres=fdr(num_discv);
        truealt=(true_gamma(:,tisind)==1);
    discv=(lfdr<=fdrthres);
    num_alt=sum(truealt)
    num_discv=sum(discv)
    overlap=sum(truealt & discv)
    discvratio=overlap/num_alt
    obsfdr=1-overlap/num_discv
end;


% test H0 ~=1111
nullind=1:15;
lfdr=sum(Post(:,nullind),2);
[fdr,~]=sort(lfdr); % fdr ascending
num_discv=find(cumsum(fdr)./[1:length(fdr)]' > FDRthres, 1, 'first')-1;
fdrthres=fdr(num_discv);
    truealt=(sum(true_gamma,2)==4);
discv=(lfdr<=fdrthres);
num_alt=sum(truealt)
num_discv=sum(discv)
overlap=sum(truealt & discv)
discvratio=overlap/num_alt
obsfdr=1-overlap/num_discv





%% increasing line plot
FDRthres=0.05;
ttcand=[1,2,3,4];
num_fd=[]; % num of false discvr
num_discovery=[];
num_truealt=sum(true_gamma(:,ttcand(1))==1);

tt=ttcand;
tt1=tt(1); % target tissue index
for k=1:size(tt,2) % k-tissue model
    % set data
    t=sort(tt(1:k)); % current tissue index for k-tissue model
    loc=find(t==tt1); % locate the target tissue
    DF=df(t);
    X=Fcormat(:,t);

    % load param
    % marg P
    tempgamma=gamma(:,t); % extract columns corresp with model2
    tempindex=sum(bsxfun(@times, tempgamma, 2.^[(k-1):-1:0]),2);
    P_sub = grpstats(P_opt,tempindex,'sum');
    % marg Delta and Sigma
    Sigma_sub=Sigma_opt(t,t);
    Delta_sub=Delta_opt(t,t);
    Mu_sub=zeros(k,1);


    % get lfdr for testing H0 (memory efficient)
    gamma_sub=double(num2str(dec2bin(0:(2^k)-1))=='1'); % 2^k * k
    sumind_0=find(gamma_sub(:,loc)==0); % index set of H0 categories
    fdr=zeros(size(X,1),1);
    for tempind=sumind_0'
        fdr=fdr+MAposterior2_3(X,Delta_sub,Mu_sub,Sigma_sub,P_sub,tempind);
    end;
    fdr=fdr./MAlikelihood2(X,Delta_sub,Mu_sub,Sigma_sub,P_sub);
    [sfdr,~]=sort(fdr); % ascending
    num_discv=find(cumsum(sfdr)./[1:length(sfdr)]'>FDRthres,1,'first')-1;
    if (isempty(num_discv))
        num_discv=length(X);
    end;
    fdrthres=sfdr(num_discv);
    clear sfdr;

    % for simulation study
    eQTLind=(fdr<=fdrthres);
    num_fd=[num_fd,sum(eQTLind & true_gamma(:,tt1)==0)];
    num_discovery=[num_discovery,sum(eQTLind)];
    clear fdr X
end;

FDR=num_fd./num_discovery % FDR 0.0500 0.0499 0.0496 0.0493
TPR=(num_discovery-num_fd)/num_truealt % TPR 0.2753 0.3475 0.3806 0.4045






