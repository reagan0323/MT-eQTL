%% Simulation Study 
%
% Simulate z-statistics from MT model (K=9) with parameters
%  from pilot 9-tissue analysis
% Estimate parameters with MT and HT (test for transf, marg est, multi-probit, etc)
% compare parameter estimates
% Compare eQTL discoveries (total number, overlap, ROC, etc)
%
% By Gen Li 12/4/2016

%% generate Z-stat data matrix
% size
N=1E5; % a reduced size
K=9;
gamma=double(num2str(dec2bin(0:(2^K)-1))=='1'); % each row is a config, 2^K rows

load('MT9_Param_1_2_3_4_5_6_7_8_9');
P_True=P_MT;
Delta_True=Delta_MT;
Sigma_True=Sigma_MT;


% Generate Zstat matrix
nsub=ceil(N*P_True); 
nsub(1)=nsub(1)-sum(nsub)+N; % cancel round error
Zstat=[];
Tlabel=[]; % true eQTL label, 1 to 2^K, config=gamma(Tlabel,:)
for i=1:(2^K)
    tempgamma=gamma(i,:);
    Zstat=[Zstat;...
        mvnrnd(zeros(1,K),Delta_True+Sigma_True.*(tempgamma'*tempgamma),nsub(i))];
    Tlabel=[Tlabel;i*ones(nsub(i),1)];
end;




%% Fit a sequence of MT-eQTL models
MT_Time=zeros(K,1);
for i=1:K
    tstart=tic;
    X=Zstat(:,1:i);
    DF=2*ones(i,1);
    [Delta_opt,Sigma_opt,P_opt,~]=EM_MT(X,DF);
    save(['Param_MT_1to',num2str(i)],'Delta_opt','Sigma_opt','P_opt');
    T=toc(tstart);
    MT_Time(i)=T;
end;
P_MT=P_opt;Delta_MT=Delta_opt;Sigma_MT=Sigma_opt;
save(['Param_MT.mat'],'P_MT','Delta_MT','Sigma_MT');


%% Fit HT-eQTL for 9 tissues
[Delta_HT,Sigma_HT,P_HT,P_HTtc,non0ind_P_HTtc]=HT_pipeline(Zstat);
save(['Param_HT.mat'],'Sigma_HT','Delta_HT','P_HT','P_HTtc','non0ind_P_HTtc');







%% Get naive TBT Pval (H0: N(0,1))
Pval=zeros(N,K);
for t=1:K
    X=Zstat(:,t);
    Pval(:,t)=normcdf(-abs(X),0,1)*2; % two-sided p value based on z~N(0,1) under H0
end;
figure(1);clf;
hist(Pval(:,1),100); % should have enrichment at 0 and uniform otherwise

save(['Pval_TBT'],'Pval');



%% Evaluate Parameter Estimates

% compare Delta (entrywise diff)
temp1=100*(Delta_HT-Delta_True)./Delta_True ;
 quantile(abs(temp1(:)),[.25 .50 .75])
temp2=100*(Delta_MT-Delta_True)./Delta_True ;
 quantile(abs(temp2(:)),[.25 .50 .75])



% compare Sigma
% diagonal
temp1=100*(Sigma_HT-Sigma_True)./(Sigma_True) ;
 quantile(abs(temp1(:)),[.25 .50 .75])
temp2=100*(Sigma_MT-(Sigma_True))./(Sigma_True) ;
 quantile(abs(temp2(:)),[.25 .50 .75])

% corr
[~,temp1]=cov2corr(Sigma_True);
[~,temp2]=cov2corr(Sigma_HT);
[~,temp3]=cov2corr(Sigma_MT);
100*(temp2-temp1)./temp1
100*(temp3-temp1)./temp1





% compare P
% KL divergence 
sum(P_True.*log(P_True./P_MT)) 
sum(P_True.*log(P_True./P_HT)) 
sum(P_True.*log(P_True./P_HTtc))





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare eQTL discoveries
% We consider the following testings, for each we compare MT model with 
% param from True, MT, and HT, and TBT method using ROC
% H0: 0000; ~=1111; 0000&1111; 0~~~;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate eQTL posterior matrices for True, MT and HT
Post_True=MAposterior2_1(Zstat,Delta_True,zeros(K,1),Sigma_True,P_True); 
LFDR_True=[Post_True(:,1), sum(Post_True(:,1:(end-1)),2), Post_True(:,1)+Post_True(:,end)];
for t=1:K
    tempgamma=gamma(:,t);
    nullind=find(1-tempgamma);
    LFDR_True=[LFDR_True,sum(Post_True(:,nullind),2)];
end;
clear Post_True;
size(LFDR_True) % from left to right, lfdr for H0: 0000; ~1111; 0000&1111; 0~~~; ~0~~; ...
save(['LFDR_True'],'LFDR_True');


Post_MT=MAposterior2_1(Zstat,Delta_MT,zeros(K,1),Sigma_MT,P_MT);   
LFDR_MT=[Post_MT(:,1), sum(Post_MT(:,1:(end-1)),2), Post_MT(:,1)+Post_MT(:,end)];
for t=1:K
    tempgamma=gamma(:,t);
    nullind=find(1-tempgamma);
    LFDR_MT=[LFDR_MT,sum(Post_MT(:,nullind),2)];
end;
clear Post_MT;
size(LFDR_MT) % from left to right, lfdr for H0: 0000; ~1111; 0000&1111; 0~~~; ~0~~; ...
save(['LFDR_MT'],'LFDR_MT');



Post_HT=MAposterior2_1(Zstat,Delta_HT,zeros(K,1),Sigma_HT,P_HT);   
LFDR_HT=[Post_HT(:,1), sum(Post_HT(:,1:(end-1)),2), Post_HT(:,1)+Post_HT(:,end)];
for t=1:K
    tempgamma=gamma(:,t);
    nullind=find(1-tempgamma);
    LFDR_HT=[LFDR_HT,sum(Post_HT(:,nullind),2)];
end;
clear Post_HT;
size(LFDR_HT) % from left to right, lfdr for H0: 0000; ~1111; 0000&1111; 0~~~; ~0~~; ...
save(['LFDR_HT'],'LFDR_HT');





%% Compare different methods in various tests
% Test1: H0:0000
testind=1; filename='H0_0000';
testname='in Any Tissue';
Nullind=1; % H0 index in Tlabel, a subset of {1,2,...,2^K}
TBTstat=min(Pval,[],2);% test stat

% % Test 2: H0: ~1111 
% testind=2; filename='H1_1111';
% testname='in All Tissues';
% Nullind=1:(2^K-1); % H0 index in Tlabel, a subset of {1,2,...,2^K}
% TBTstat=max(Pval,[],2);% test stat

% % Test 3: H0: 0000&1111 
% testind=3; filename='H0_0000_1111';
% testname='a Subset of Tissues';
% Nullind=[1,2^K]; % H0 index in Tlabel, a subset of {1,2,...,2^K}
% TBTstat=min(Pval,[],2)+1-max(Pval,[],2);% test stat
% 
% % Test 4~K+3: H0: 0~~~~ 
% % for t=3+(1:K)
% t=4;
% testind=t; 
% filename=['H0_Tissue_',num2str(t-3)];
% testname=['in Tissue ',num2str(t-3)];
% Nullind=find(~gamma(:,t-3)); % H0 index in Tlabel, a subset of {1,2,...,2^K}
% TBTstat=Pval(:,t-3);% test stat
% 




% Below are universal for different tests
Tlabel_t=~ismember(Tlabel,Nullind); % 0/1 vector, 1=true positive
n1_t=sum(Tlabel_t); % number of true positives
n0_t=N-n1_t; % number of true negatives

% TBT (p value)
[~,I]=sort(TBTstat,'ascend'); % sort min-pval-stat
clear TBTstat;
label_Pval=Tlabel_t(I); % sort 0-1 label accordingly
clear I;
x_Pval=cumsum(1-label_Pval)./n0_t; % false positive rate: #FP/#TotalNeg (not FDR)
y_Pval=cumsum(label_Pval)./n1_t; % true positive rate: #TP/#TotalPos
clear label_Pval

% True
fdr=LFDR_True(:,testind); 
[~,I]=sort(fdr,'ascend'); % ascending
clear fdr;
label=Tlabel_t(I); % sort 0-1 label accordingly
clear I;
x_True=cumsum(1-label)./n0_t;
y_True=cumsum(label)./n1_t;
clear label;

% MT
fdr=LFDR_MT(:,testind); 
[~,I]=sort(fdr,'ascend'); % ascending
clear fdr;
label=Tlabel_t(I); % sort 0-1 label accordingly
clear I;
x_MT=cumsum(1-label)./n0_t;
y_MT=cumsum(label)./n1_t;
clear label;

% HT
fdr=LFDR_HT(:,testind); 
[~,I]=sort(fdr,'ascend'); % ascending
clear fdr;
label=Tlabel_t(I); % sort 0-1 label accordingly
clear I;
x_HT=cumsum(1-label)./n0_t;
y_HT=cumsum(label)./n1_t;
clear label;

save(['ROC_coordinates_',filename],'x_Pval','y_Pval','x_True','y_True','x_MT','y_MT','x_HT','y_HT');



%% ROC curve figures -- paper
filename='H0_0000';
% filename='H1_1111';
% filename='H0_0000_1111';
% filename='H0_Tissue_1';
load(['ROC_coordinates_',filename]);
% Plot ROC Curves
figure();clf;
plot(x_True,y_True,'k-','linewidth',3);
hold on;
plot(x_HT,y_HT,'r-.','linewidth',3);
plot(x_MT,y_MT,'b--','linewidth',3);
plot(x_Pval,y_Pval,'g:','linewidth',3);
plot([0,1],[0,1],'k--','linewidth',2);
h=legend('Oracle','HT-eQTL','MT-eQTL','TBT','location','SouthEast'); % "oracle" means we use the true parameters and the true generative model (MT) to detect eQTL, using lfdr approach
set(h,'fontsize',35);
set(gca,'fontsize',30);
xlabel('False Positive Rate','fontsize',35);
ylabel('True Positive Rate','fontsize',35);
title(['(a) Any eQTL Detection'],'fontsize',35);
% title(['(b) Common eQTL Detection'],'fontsize',35);
% title(['(c) Tissue-Specific eQTL Detection'],'fontsize',35);
% title(['(d) Single-Tissue eQTL Detection'],'fontsize',35);
orient landscape;
print('-dpdf',['ROC_',filename]);
