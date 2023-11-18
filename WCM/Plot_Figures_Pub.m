function Plot_Figures_Pub

%% UBIQUITOUS CODE

% Path where figures will be saved
figpath='figures/';
mkdir(figpath)
% Path where matfiles are saved
matpath='matfiles/';

% Ubiquitous color data
colororder=[0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

% RunPrep
[dataS,dataG]=RunPrep;

set(0,'DefaultAxesFontSize',8)


%% Figure 2B
% Example stochastic data

%Get data
St=load(strcat(matpath,'STE_1D.mat')); 
xoutG_all=St.xoutG_all;
xoutS_all=St.xoutS_all;
tout=St.tout_all;
[~,dataG]=RunPrep;
I=Plot_GetInfo;
tERK=I.tERK;
    
%Plot genes; MAPK1, MAPK3
g1=figure; 
subplot(4,1,1)
plot(tout/3600,xoutG_all(:,116),'Color',colororder(1,:)) %ERK1
axis tight
xlim([0 24])
set(gca,'XTick',0:4:24)
subplot(4,1,2)
plot(tout/3600,xoutG_all(:,117),'Color',colororder(2,:)) %ERK2
axis tight
xlim([0 24])
set(gca,'XTick',0:4:24)

%Plot mRNA over time; MAPK1, MAPK3
subplot(4,1,3)
plot(tout/3600,xoutG_all(:,[398,399]))
axis tight
xlim([0 24])
set(gca,'XTick',0:4:24)

%Plot protein over time stochastic; ERK
subplot(4,1,4)
purple=mean([1 0 0;0 0 1]);
plot(tout/3600,sum(xoutS_all(:,find(tERK)),2),'Color',purple)
axis tight
xlim([0 24])
set(gca,'XTick',0:4:24)

%% Figure 2C
% CV Histogram

[dataS,~]=RunPrep;
St=load(strcat(matpath,'RandomPopCells_100cells.mat'));
for m=1:length(St.cells0)
    xoutG=St.cells0{m}.xoutG_all;
    xoutS=St.cells0{m}.xoutS_all;
    Obs(:,m) = GetObservables(xoutS(1,:),dataS.VxPARCDL,dataS.kS(3));
    mRNAs(:,m)=xoutG([283:287,292:294,314:337,340:345,347:348,350:415,418:422]);
end
s=std(Obs');
m=mean(Obs');
stds=s./m;
h1=figure; histogram(stds(find(stds>0)),'FaceColor','none','EdgeColor','r')
hold on
s=std(mRNAs');
m=mean(mRNAs');
stds=s./m;
histogram(stds(and(stds>0,~isnan(stds))),'FaceColor','none','EdgeColor','b')

%% Figure S2B
% Test time steps; phenotypic test

timesteps=[0.1,0.5,1,5,10]*60;
St=load(strcat(matpath,'TestTimeStep.mat'));
cells=St.cells;
for k=1:length(timesteps)-1
    for i=1:length(cells)
        [x,m]=min(find(cells{i,k}.xoutS_all(:,69)>25));
        if x
            ms(k,i)=cells{i,k}.tout_all(x);
        else ms(k,i)=NaN;
        end
    end
end
%plot
tt=figure;
boxplot(ms'/3600)
ylim([0, Inf])
set(gca,'XTickLabels',timesteps)

%% Figure S2C
% Effect of EIF4E on Translation Rate

k50E=100;
EIF4E=0:1000;
vTL=[];

for i=1:length(EIF4E)
    vTL(i)=EIF4E(i)/(k50E+EIF4E(i));
end

a=figure; 
plot(EIF4E,vTL,'k-','LineWidth',2)

%% Figure S2D
% Ribosome levels doubling in ~20 hours

%Get data
St=load(strcat(matpath,'Ribosome_dynamics.mat')); 
xoutS_all=St.xoutS_all;
datasheet='master_MCF10A.xlsx';
[cellparams,~,~]=xlsread(datasheet,'cellparams','','basic');
VolumeofCell=cellparams(1);
[~,~,~,mpc2nmcf_Vc,~,~]=CalcVolumeParams(VolumeofCell);

% Plot and convert to molecules per cell
a=figure; plot(tout_all/3600,xoutS_all(:,1)*(1/mpc2nmcf_Vc),'k-','LineWidth',2)
axis tight

%% Figure S2E
% Global correlation between all proteins stimulated with growth factors
% Need to run Fig 6F-D to get this data.

% Get data
files={'CellCycle2S_E+I_400cells.mat'};
St=load(strcat(matpath,files{1})); 
cells=St.cells;

% Compute correlations between all proteins
[dataS,dataG]=RunPrep;
VxPARCDL=dataS.VxPARCDL;
Vc=dataS.kS(3);
Obs_mat=[];
Obs_mat_vcorr=[];
for i=1:length(cells)
    xout=cells{i}.xoutS_all;
    Obs = GetObservables_matrix(xout,VxPARCDL,Vc);
    Obs_mat_vcorr(:,i)=Obs(end,:)*cells{i}.xoutS_all(end,1);
    Obs_mat(:,i)=Obs(end,:);
    disp(i)
end

corr_mat_vcorr=corrcoef(Obs_mat_vcorr');
% plot without ones that are zero
indnonzero=find(nansum(corr_mat_vcorr)~=0);
a=figure; imagesc(corr_mat_vcorr(indnonzero,indnonzero),[-1 1])
colorbar

set(0,'DefaultAxesFontSize',8)
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])

%% Figure S2F

datasheet='data_experimental.xlsx';
[nums,~,~]=xlsread(datasheet,'ShiProtein','','basic');

figure; loglog(nums(:,1),nums(:,2),'.','MarkerSize',10)
ylim([500 1000000])
xlim([500 1000000])
refline(1,0)

%% Figure 3A
% Ligand-Receptor Cooperativity

hillequation=@(x,L) (x(2).*(L.^x(1)))./((x(3).^x(1))+(L.^x(1)));

set(0,'DefaultAxesFontSize',8)
    
% LIGAND-MATCHED RECEPTOR EXPRESSION
%Import data
St=load(strcat(matpath,'Cooperativity1D_sameR.mat'));
conds=St.conds;
%conds is a cell array matrix. Each row is a different growth factor and
%each column is a different dose of that growth factor.

% Cooperativity - D
indsS = 156:162; 
indsSP = 189:194;
indsSPP = 195:214;
indsSPSP = 215:225;

colors=[[0 0 1];[0 0 1];[0 0 0];[1 0 0];[0 0 1];[0 0 1];[0 0 1]];
dlmwrite(strcat(figpath,'ns.txt'),'')

a=figure;
for i=1:size(conds,1) %Different ligands (rows)   
SPP=[];
SP=[];
SPSP=[];
S=[];
for k=1:size(conds,2) %Different doses
    SPP(:,k)=sum(conds{i,k}.xoutS_all(:,indsSPP),2);
    SP(:,k)=sum(conds{i,k}.xoutS_all(:,indsSP),2);
    SPSP(:,k)=sum(conds{i,k}.xoutS_all(:,indsSPSP),2);
    S(:,k)=sum(conds{i,k}.xoutS_all(:,indsS),2);
end
%Plot cooperativity
lig=S(1,:);
S=S(end,:);
SP=SP(end,:);
SPP=SPP(end,:);
SPSP=SPSP(end,:);

% Dose Response
ns_prod_sum = SP + SPP + 2*SPSP;
%To calculate the Hill coefficient (n) of the dose-response curve.
options = optimoptions('lsqcurvefit','MaxFunEvals',5000,'MaxIter',5000,'TolFun',1e-14,'TolX',1e-14);
[ns_out,~] = lsqcurvefit(@(x,L) (x(2).*(L.^x(1)))./((x(3).^x(1))+(L.^x(1))),[1;max(ns_prod_sum);1],lig,ns_prod_sum,[0,0,0],[inf,inf,inf],options);
ns = ns_out(1);
subplot(3,3,i)
semilogx(lig,ns_prod_sum,'.','Color',colors(i,:),'MarkerSize',10)
hold on
semilogx(lig,hillequation(ns_out,lig),'-','Color',colors(i,:),'LineWidth',1)
text(min(lig)*2,max(ns_prod_sum)*0.8,strcat('n=',num2str(roundn(ns,-1))),'FontSize', 6)
dlmwrite(strcat(figpath,'ns.txt'),roundn(ns,-1),'-append')
axis tight
xlim([0.001,1000])
set(gca,'XTick',[0.001,0.1,10,1000])
end

%% Figure 3B-D

% ppERK simulation
I=Plot_GetInfo;
tppERK=I.tppERK;
inds=tppERK;
a1=figure;
for n=1:3 %three conditions
StD=load(strcat(matpath,'GFdoseresponse1D_EI_',num2str(n),'.mat'));
condsD=StD.condsD;
colors=colororder;%cool(4);
% timepoints plots Determ
subplot(3,1,n)
    for i=1:length(condsD)
        x=condsD{i}.tout_all/60;
        y=sum(condsD{i}.xoutS_all(:,inds).*repmat([dataS.VxPARCDL(inds)/dataS.kS(3)]',length(x),1),2);
        plot(x,y,'Color',colors(i,:));
        hold on
        xlim([0 360])
        ylim([0 120])
        set(gca,'XTick',0:120:360)
    end
end
% pERK uWestern data
datasheet='data_experimental.xlsx';
[nums,~,~]=xlsread(datasheet,'uWestern_091516','','basic');
p_avg=nums(:,20);
p_std=nums(:,21);
time=[0,5,30,180,360]; %minutes
colors=colororder;
% rows=doses or conditions, columns=timepoints
p_mat_avg=[p_avg(1:15),p_avg(17:17+14),p_avg(33:33+14),p_avg(49:49+14)];
p_mat_std=[p_std(1:15),p_std(17:17+14),p_std(33:33+14),p_std(49:49+14)];
p_mat_avg_norm=p_mat_avg./repmat(p_mat_avg(13,:),15,1);
p_mat_std_norm=p_mat_std./repmat(p_mat_avg(13,:),15,1);
p_mat_avg_norm=[ones(size(p_mat_avg_norm,1),1),p_mat_avg_norm];
p_mat_std_norm=[repmat(p_mat_std_norm(13,1),size(p_mat_std_norm,1),1),p_mat_std_norm];
egfrows=1:4;
insrows=5:8;
bcrows=9:12;
% Timecourses
a2=figure;
numsets=3;
for k=1:numsets
    if k==1; r=egfrows; end
    if k==2; r=insrows; end
    if k==3; r=bcrows; end
    subplot(numsets,1,k)
    for i=1:length(r)
        x=time;
        yavg=p_mat_avg_norm(r(i),:);
        ystd=p_mat_std_norm(r(i),:);
        errorbar(x,yavg,ystd,ystd,'Color',colors(i,:));
        hold on
        plot(x,yavg,'.-','MarkerSize',10,'Color',colors(i,:));
        xlim([0 360])
        ylim([0 5])
        set(gca,'XTick',0:120:360)
        grid on
    end
end


% ppAKT simulation
I=Plot_GetInfo;
tppAKT=I.tppAKT;
inds=tppAKT; 
b1=figure;
for n=1:3 %three conditions
StD=load(strcat(matpath,'GFdoseresponse1D_EI_',num2str(n),'.mat'));
condsD=StD.condsD;
colors=colororder;
% timepoints plots Determ
subplot(3,1,n)
    for i=1:length(condsD)
        x=condsD{i}.tout_all/60;
        y=sum(condsD{i}.xoutS_all(:,inds).*repmat([dataS.VxPARCDL(inds)/dataS.kS(3)]',length(x),1),2);
        plot(x,y,'Color',colors(i,:));
        hold on
        xlim([0 360])
        ylim([0 35])
        set(gca,'XTick',0:120:360)
    end
end
% pAKT uWestern data
datasheet='data_experimental.xlsx';
[nums,~,~]=xlsread(datasheet,'uWestern_091516','','basic');
p_avg=nums(:,42);
p_std=nums(:,43);
time=[0,5,30,180,360]; %minutes
colors=colororder;
% rows=doses or conditions, columns=timepoints
p_mat_avg=[p_avg(1:15),p_avg(17:17+14),p_avg(33:33+14),p_avg(49:49+14)];
p_mat_std=[p_std(1:15),p_std(17:17+14),p_std(33:33+14),p_std(49:49+14)];
p_mat_avg_norm=p_mat_avg./repmat(p_mat_avg(13,:),15,1);
p_mat_std_norm=p_mat_std./repmat(p_mat_avg(13,:),15,1);
p_mat_avg_norm=[ones(size(p_mat_avg_norm,1),1),p_mat_avg_norm];
p_mat_std_norm=[repmat(p_mat_std_norm(13,1),size(p_mat_std_norm,1),1),p_mat_std_norm];
egfrows=1:4;
insrows=5:8;
bcrows=9:12;
% Timecourses
b2=figure;
numsets=3;
for k=1:numsets
    if k==1; r=egfrows; end
    if k==2; r=insrows; end
    if k==3; r=bcrows; end
    subplot(numsets,1,k)
    for i=1:length(r)
        x=time;
        yavg=p_mat_avg_norm(r(i),:);
        ystd=p_mat_std_norm(r(i),:);
        errorbar(x,yavg,ystd,ystd,'Color',colors(i,:));
        hold on
        plot(x,yavg,'.-','MarkerSize',10,'Color',colors(i,:));
        xlim([0 360])
        ylim([0 50])
        set(gca,'XTick',0:120:360)
        grid on
    end
end


% mTOR1D
% 4EBP1 simulation
inds=727; %pEIF4EBP1
c1=figure;
for n=1:3 %three conditions
StD=load(strcat(matpath,'GFdoseresponse1D_EI_',num2str(n),'.mat'));
condsD=StD.condsD;
colors=colororder;
% timepoints plots Determ
subplot(3,1,n)
    for i=1:length(condsD)
        x=condsD{i}.tout_all/60;
        y=sum(condsD{i}.xoutS_all(:,inds).*repmat([dataS.VxPARCDL(inds)/dataS.kS(3)]',length(x),1),2);
        plot(x,y,'Color',colors(i,:));%,'LineWidth',0.3)
        hold on
        xlim([0 360])
        ylim([0 10])
        set(gca,'XTick',0:120:360)
    end
end
% 4EBP1 uWestern data
datasheet='data_experimental.xlsx';
[nums,~,~]=xlsread(datasheet,'uWestern_110216','','basic');
p_avg=nums(:,37);
p_std=nums(:,38);
time=[0,5,30,180,360]; %minutes
colors=colororder;
% rows=doses or conditions, columns=timepoints
p_mat_avg=[p_avg(1:15),p_avg(17:17+14),p_avg(33:33+14),p_avg(49:49+14)];
p_mat_std=[p_std(1:15),p_std(17:17+14),p_std(33:33+14),p_std(49:49+14)];
p_mat_avg_norm=p_mat_avg./repmat(p_mat_avg(13,:),15,1);
p_mat_std_norm=p_mat_std./repmat(p_mat_avg(13,:),15,1);
p_mat_avg_norm=[ones(size(p_mat_avg_norm,1),1),p_mat_avg_norm];
p_mat_std_norm=[repmat(p_mat_std_norm(13,1),size(p_mat_std_norm,1),1),p_mat_std_norm];
egfrows=1:4;
insrows=5:8;
bcrows=9:12;
% Timecourses
c2=figure;
numsets=3;
for k=1:numsets
    if k==1; r=egfrows; end
    if k==2; r=insrows; end
    if k==3; r=bcrows; end
    subplot(numsets,1,k)
    for i=1:length(r)
        x=time;
        yavg=p_mat_avg_norm(r(i),:);
        ystd=p_mat_std_norm(r(i),:);
        errorbar(x,yavg,ystd,ystd,'Color',colors(i,:));
        hold on
        plot(x,yavg,'.-','MarkerSize',10,'Color',colors(i,:));
        xlim([0 360])
        ylim([0 3.5])
        set(gca,'XTick',0:120:360)
        grid on
    end
end

%% Figure 3E
% DNA Damage Deterministic Simulations

colors=[0 0 1;1 0 0];
a=figure;
for k=1:3
St=load(strcat(matpath,'DNADamage1D_',num2str(k),'.mat'));
condsD=St.condsD;
subplot(3,1,k)
    for i=1:length(condsD)
        plot(condsD{i}.tout_all/3600,condsD{i}.xoutS_all(:,3),'Color',colors(i,:),'LineWidth',1.2)
        axis tight
        xlim([0 24])
        set(gca,'XTick',0:6:24)
        hold on
    end
end

%% Figure 3F
% Stochastic DNA Damage Unit Testing (Bar Plots)

% Get Data
St=load(strcat(matpath,'DNADamage1S.mat'));
cells=St.cells;
for i=1:size(cells,1)
    for k=1:size(cells,2)
        xoutIN=cells{i,k}.xoutS_all(:,3);
        tout=cells{i,k}.tout_all/3600;
        [pks,locs,w]=findpeaks(xoutIN,tout,'MinPeakHeight',400,'MinPeakDistance',5);
        peaks{i,k}=[pks,w];
    end
    damages(i)=cells{i,k}.xoutS_all(1,31);
end

% Plot box plots
for i=1:size(peaks,1)
    avar=[];
    for k=1:size(peaks,2) %Cells
        avar(k,:)=mean(peaks{i,k},1);
        np(k)=length(peaks{i,k}(:,1));
    end
    hei(:,i)=avar(:,1);
    wid(:,i)=avar(:,2);
    nps(:,i)=np;
end
b2=figure; errorbar(damages,mean(nps),std(nps),'-sk','MarkerSize',5,'MarkerFaceColor','k'); set(gca,'XTick',damages); xlim([min(damages),max(damages)]); 
b3=figure; errorbar(damages,mean(hei),std(hei),'-sk','MarkerSize',5,'MarkerFaceColor','k'); set(gca,'XTick',damages); xlim([min(damages),max(damages)]); ylim([0 1500])
b4=figure; errorbar(damages,mean(wid),std(wid),'-sk','MarkerSize',5,'MarkerFaceColor','k'); set(gca,'XTick',damages); xlim([min(damages),max(damages)]); ylim([0 5])

%% Figure 3G
% TRAIL dose response Time Course

% Left
b1=figure;
St=load(strcat(matpath,'Apoptosis1D_dr.mat'));
conds=St.condsD;
colors=jet(length(conds));
tds=[];
doses=[];
for i=1:length(conds)
    x=conds{i}.tout_all/3600;
    y=conds{i}.xoutS_all(:,106);
    if ~rem(i,2); plot(x,y,'Color',colors(i,:),'LineWidth',2); hold on; end
    xlim([0 72])
    td=x(min(find(y>100)));
    if ~isempty(td); tds(i)=td; else tds(i)=72; end
    doses(i)=conds{i}.xoutS_all(1,84);
end

% Right
b2=figure;
dosesngperml=doses*2.597402597402597e+01;
FO = fit(log10(dosesngperml'), tds', 'exp2');
h=plot(FO,log10(dosesngperml),tds,'r.');
legend off; xlabel(''); ylabel('');
xlim([-2 3])
set(gca,'XTick',[-2 0 3])
set(gca,'XTickLabels',{'10^{-2}','10^{0}','10^{3}'});%[1E-2,1E-1,1E0,1E1,1E3])
set(h(1),'MarkerSize',10)
set(h(2),'Color','r')

%% Figure 3H
% TRAIL dose response Percent Cell Death

St=load(strcat(matpath,'Apoptosis1S_dr.mat'));
conds=St.condsS;
dead_cells=[];
doses=[];
for i=1:size(conds,1)
    for k=1:size(conds,2)
        x=conds{i,k}.tout_all/3600;
        y=conds{i,k}.xoutS_all(:,106);
        %plot(x,y,'Color',colors(i,:))
        hold on
        if sum(y>100)
            dead_cells(i,k)=1;
        else dead_cells(i,k)=0;
        end
    end
    doses(i)=conds{i,k}.xoutS_all(1,84);
end
dosesngperml=doses*2.597402597402597e+01;

percent_dead_cells=sum(dead_cells,2)/size(conds,2)*100;
x=log10(dosesngperml');
y=100-percent_dead_cells;
[params,x_vector,y_vector]=sigm_fit(x,y);
%Experimental data (Luis)
[traildeath_exp,~,~]=xlsread('data_experimental.xlsx','TRAILdr','','basic');
plot(log10(traildeath_exp(:,1)),traildeath_exp(:,2),'k*','MarkerSize',5)
hold on
plot(x_vector,y_vector,'r-')
plot(x,y,'r.','MarkerSize',10)
ylim([0 100])
xlim([0 3])
set(gca,'XTick',[0 1 2 3])
set(gca,'XTickLabels',{'10^{0}','10^{1}','10^{2}','10^{3}'});

%% Figure 3I
% Effect of ppERK/ppAKT and PUMA/NOXA on time-to-death

colors=[0,0,0;0 0 1;1 0 0];    
a=figure;
St=load(strcat(matpath,'Apoptosis1D.mat'));
conds=St.conds;
for i=1:length(conds)
    x=conds{i}.tout_all/3600;
    y=conds{i}.xoutS_all(:,106);
    plot(x,y,'Color',colors(i,:),'LineWidth',2)
    hold on
end
xlim([0 60])

%% Figure 3J
% Dynamics of Cell Cycle Proteins

set(0,'DefaultAxesFontSize',8)

% IMPORT DATA
St=load(strcat(matpath,'CellCycle1D.mat')); 
conds=St.conds;

sp=[46,50,59,69]; %Cyclin D, E, A, and B, active.
a=figure;
for k=1:length(sp)
subplot(2,2,k)
for i=1:length(conds)
    x=conds{i}.tout_all/3600;
    y=conds{i}.xoutS_all(:,sp(k));    
    plot(x,y,'Color',colororder(i,:),'LineWidth',1.2)
    axis tight
    xlim([0 max(x)])
    set(gca,'XTick',0:24:96)
    hold on
end
    if k==1; set(gca,'YTick',0:40:120); end
    if k==2; set(gca,'YTick',0:20:60); end 
end

%% Figure S3A
% Ligand-Receptor Cooperativity Receptors All Highly Expressed

hillequation=@(x,L) (x(2).*(L.^x(1)))./((x(3).^x(1))+(L.^x(1)));

set(0,'DefaultAxesFontSize',8)

% LIGAND-MATCHED RECEPTOR EXPRESSION
%Import data
St=load(strcat(matpath,'Cooperativity1D_sameR.mat'));
conds=St.conds;
%conds is a cell array matrix. Each row is a different growth factor and
%each column is a different dose of that growth factor.

% Cooperativity - D
indsS = 156:162; 
indsSP = 189:194;
indsSPP = 195:214;
indsSPSP = 215:225;

colors=[[0 0 1];[0 0 1];[0 0 0];[1 0 0];[0 0 1];[0 0 1];[0 0 1]];
dlmwrite(strcat(figpath,'ns.txt'),'')

a=figure;
for i=1:size(conds,1) %Different ligands (rows)   
SPP=[];
SP=[];
SPSP=[];
S=[];
for k=1:size(conds,2) %Different doses
    SPP(:,k)=sum(conds{i,k}.xoutS_all(:,indsSPP),2);
    SP(:,k)=sum(conds{i,k}.xoutS_all(:,indsSP),2);
    SPSP(:,k)=sum(conds{i,k}.xoutS_all(:,indsSPSP),2);
    S(:,k)=sum(conds{i,k}.xoutS_all(:,indsS),2);
end
%Plot cooperativity
lig=S(1,:);
S=S(end,:);
SP=SP(end,:);
SPP=SPP(end,:);
SPSP=SPSP(end,:);

% Dose Response
ns_prod_sum = SP + SPP + 2*SPSP;
%To calculate the Hill coefficient (n) of the dose-response curve.
options = optimoptions('lsqcurvefit','MaxFunEvals',5000,'MaxIter',5000,'TolFun',1e-14,'TolX',1e-14);
[ns_out,~] = lsqcurvefit(@(x,L) (x(2).*(L.^x(1)))./((x(3).^x(1))+(L.^x(1))),[1;max(ns_prod_sum);1],lig,ns_prod_sum,[0,0,0],[inf,inf,inf],options);
ns = ns_out(1);
subplot(7,1,i)
semilogx(lig,ns_prod_sum,'.','Color',colors(i,:),'MarkerSize',10)
hold on
semilogx(lig,hillequation(ns_out,lig),'-','Color',colors(i,:),'LineWidth',1)
text(min(lig)*2,max(ns_prod_sum)*0.8,strcat('n=',num2str(roundn(ns,-1))),'FontSize', 6)
dlmwrite(strcat(figpath,'ns.txt'),roundn(ns,-1),'-append')
axis tight
xlim([0.001,1000])
set(gca,'XTick',[0.001,0.1,10,1000])
end



% MCF10A CONTEXT
%Import Data
St=load(strcat(matpath,'Cooperativity1D.mat'));
conds=St.conds;
%conds is a cell array matrix. Each row is a different growth factor and
%each column is a different dose of that growth factor.

% Cooperativity - D
indsS = 156:162; 
indsSP = 189:194;
indsSPP = 195:214;
indsSPSP = 215:225;

colors=[[0 0 1];[0 0 1];[0 0 0];[1 0 0];[0 0 1];[0 0 1];[0 0 1]];
dlmwrite(strcat(figpath,'ns.txt'),'')

a=figure;
for i=1:size(conds,1) %Different ligands (rows)   
SPP=[];
SP=[];
SPSP=[];
S=[];
for k=1:size(conds,2) %Different doses
    SPP(:,k)=sum(conds{i,k}.xoutS_all(:,indsSPP),2);
    SP(:,k)=sum(conds{i,k}.xoutS_all(:,indsSP),2);
    SPSP(:,k)=sum(conds{i,k}.xoutS_all(:,indsSPSP),2);
    S(:,k)=sum(conds{i,k}.xoutS_all(:,indsS),2);
end
%Plot cooperativity
lig=S(1,:);
S=S(end,:);
SP=SP(end,:);
SPP=SPP(end,:);
SPSP=SPSP(end,:);

% Dose Response
ns_prod_sum = SP + SPP + 2*SPSP;
%To calculate the Hill coefficient (n) of the dose-response curve.
options = optimoptions('lsqcurvefit','MaxFunEvals',5000,'MaxIter',5000,'TolFun',1e-14,'TolX',1e-14);
[ns_out,~] = lsqcurvefit(@(x,L) (x(2).*(L.^x(1)))./((x(3).^x(1))+(L.^x(1))),[1;max(ns_prod_sum);1],lig,ns_prod_sum,[0,0,0],[inf,inf,inf],options);
ns = ns_out(1);
subplot(7,1,i)
semilogx(lig,ns_prod_sum,'.','Color',colors(i,:),'MarkerSize',10)
hold on
semilogx(lig,hillequation(ns_out,lig),'-','Color',colors(i,:),'LineWidth',1)
text(min(lig)*2,max(ns_prod_sum)*0.8,strcat('n=',num2str(roundn(ns,-1))),'FontSize', 6)
dlmwrite(strcat(figpath,'ns.txt'),roundn(ns,-1),'-append')
axis tight
xlim([0.001,1000])
set(gca,'XTick',[0.001,0.1,10,1000])
end

%% Figure S3B
% Phosphorylated receptor dynamics

%Get data
St=load(strcat(matpath,'Phos_Receptor_Dynamics.mat')); 
xoutS_all=St.xoutS_all;
tout_all=St.tout_all;

%Plot
inds_Rint=[291:352,382:410,498:526];
inds_R=[229:261,353:381,411:497,527:642];
a=figure;
plot(tout_all/60,sum(xoutS_all(:,inds_R),2),'k-','LineWidth',2)
hold on
plot(tout_all/60,sum(xoutS_all(:,inds_Rint),2),'k--','LineWidth',2)
axis tight
set(gca,'XTick',0:20:120)

%% Figure S3C
% ppERK and ppAKT Dynamic Dose Response

I=Plot_GetInfo;
tppERK=I.tppERK;
tppAKT=I.tppAKT;

St=load(strcat(matpath,'GFDoseResponse.mat'));
cells=St.cells;

colors=jet(size(cells,2));

a=figure;
for i=1:size(cells,1)
    subplot(size(cells,1),1,i)
    for k=1:size(cells,2)
        x=cells{i,k}.tout_all/60;
        y=sum(cells{i,k}.xoutS_all(:,tppERK).*repmat([dataS.VxPARCDL(tppERK)/dataS.kS(3)]',length(x),1),2);
        plot(x,y,'Color',colors(k,:),'LineWidth',1.2)
        hold on
        axis tight
        ylim([0 130])
        xlim([0 120])
        set(gca,'XTick',[0:40:120])
    end
end

b=figure;
for i=1:size(cells,1)
    subplot(size(cells,1),1,i)
    for k=1:size(cells,2)
        x=cells{i,k}.tout_all/60;
        y=sum(cells{i,k}.xoutS_all(:,tppAKT).*repmat([dataS.VxPARCDL(tppAKT)/dataS.kS(3)]',length(x),1),2);
        plot(x,y,'Color',colors(k,:),'LineWidth',1.2)
        hold on
        axis tight
        ylim([0 40])
        xlim([0 120])
        set(gca,'XTick',[0:40:120])
    end
end

%% Figure S3E
% Stochastic DNA Damage Unit Testing (Time courses)

St=load(strcat(matpath,'DNADamage1S.mat'));
cells=St.cells;
b1=figure;
for i=1:size(cells,1)
    subplot(size(cells,1),1,i)
    for k=1:size(cells,2)
        xoutIN=cells{i,k}.xoutS_all(:,3);
        tout=cells{i,k}.tout_all/3600;
        [pks,locs,w]=findpeaks(xoutIN,tout,'MinPeakHeight',400,'MinPeakDistance',5);
        plot(tout,xoutIN)
        hold on
        ylim([0 1500])
        peaks{i,k}=[pks,w];
    end
    damages(i)=cells{i,k}.xoutS_all(1,31);
end

%% Figure S5

% A
% Experimental Data
St=load([matpath,'CellCycle_E+I_HeteroExp_AllCells.mat']);
cells=St.AllCells;
cells=cells(1:end-2); %20 and 21 are outliers
x=[];
y=[];
for i=1:length(cells)
    x=(cells{i}(:,1)-1)*3-15;
    y(:,i)=cells{i}(:,2);
end
a=figure; plot(x,y)
axis tight
xlim([0 40])
   
% Simulation Data
St=load([matpath,'CellCycle_E+I_HeteroSims.mat']);
cells=St.cells;
x=[];
y=[];
youts=zeros(length(cells{1}.tout_all),length(cells));
inds_tppERK=[151   155   226   227   228   676   718   729   730   733   738   740   742   750   757   758];
b=figure;
for i=1:length(cells)
    x=cells{i}.tout_all/60;
    y=sum(cells{i}.xoutS_all(:,inds_tppERK).*repmat([dataS.VxPARCDL(inds_tppERK)/dataS.kS(3)]',length(x),1),2);
    youts(1:length(x),i)=y;

    plot(x,y)
    hold on
end
axis tight
xlim([0 40])



% B (controls)
St1=load([matpath,'CellCycle_E+I_HeteroExp_AllCells.mat']);
St2=load([matpath,'CellCycle_E+I_Hetero_CtrlMock.mat']);
St3=load([matpath,'CellCycle_E+I_Hetero_CtrlMEKi100nM.mat']);

c=figure;
% E+I
subplot(1,3,1)
cells=St1.AllCells;
cells=cells(1:end-2); %20 and 21 are outliers
x=[];
y=[];
for i=1:length(cells)
    x=(cells{i}(:,1)-1)*3-15;
    y(:,i)=cells{i}(:,2);
end
plot(x,y)
hold on
plot(x,mean(y,2),'--k','LineWidth',2)
axis tight
xlim([0 40])
ylim([0.5 1.7])

% Mock
subplot(1,3,2)
xIN=(St2.x-1)*3-15;
yIN=St2.y;
plot(xIN,yIN)
hold on
plot(xIN,mean(yIN,2),'--k','LineWidth',2)
axis tight
xlim([0 40])
ylim([0.5 1.7])


% 100nM PD (MEKi)
subplot(1,3,3)
xIN=(St3.x-1)*3-15;
yIN=St3.y;
plot(xIN,yIN)
hold on
plot(xIN,mean(yIN,2),'--k','LineWidth',2)
axis tight
xlim([0 40])
ylim([0.5 1.7])

%% Figure 4
% DNA Damage Simulations

% Time courses
ttds_all=[];
for i=1:3 %Number of files
a(i)=figure;
if i==1; St=load(strcat(matpath,'DNADamageSerum')); end
if i==2; St=load(strcat(matpath,'DNADamageNoSerum')); end
if i==3; St=load(strcat(matpath,'DNADamageSerum_p21mut')); end
cells=St.cells;
inds2plot=[3,59,106];
ttds=[];
for k=1:length(inds2plot) %files/conditions
    indspecies=inds2plot(k);
    subplot(3,1,k)
    numcells=length(cells);
    for n=1:numcells %cells
        xoutS_all=cells{n}.xoutS_all;
        tout_all=cells{n}.tout_all;
        y=xoutS_all(:,indspecies);
        plot(tout_all/3600,y)
        axis tight 
        hold on
        if k==3
        % Calc ttd
        time=find(xoutS_all(:,106)>100);
        if ~isempty(time)
            ttds(n,1)=tout_all(min(time))/3600;
        else ttds(n,1)=Inf;
        end
        end
    end
    if k==1; ylim([0 2000]); end
    if k==2; ylim([0 40]); end
    if k==3; end
    if i==1
    end
    if i==2
    xlim([0 72]);
    set(gca,'XTick',0:24:72)
    set(gca,'XTickLabels',0:24:72)
    end    
end
ttds_all(:,i)=ttds;
end
ttds_all(:,[1,3])=ttds_all(:,[1,3])-72;


% Bar graphs
% Prep
d=[];
d(1,:)=sum(ttds_all<=24)/numcells*100;
d(2,:)=sum(ttds_all<=48)/numcells*100;
d(3,:)=sum(ttds_all<=72)/numcells*100;
dsem=sqrt((d.*(100-d))/numcells);
%%% WITH GFs
[death1_exp,~,~]=xlsread('data_experimental.xlsx','death1','','basic');
death_exp_avg=death1_exp(9,13:15);
death_exp_std=death1_exp(9,19:21);
death_exp_sem=death_exp_std/sqrt(numcells);
% Experimental
b1(1)=figure;
bp=barwitherr(death_exp_sem, death_exp_avg);
set(bp,'FaceColor',[0.5,0.5,0.5]);
set(gca,'YTick',[0:20:100])
set(gca,'XTickLabels','')
axis tight
ylim([0 100])
grid on
% Simulation
b1(2)=figure;
bp=barwitherr(dsem(1:3,1),d(1:3,1));
set(gca,'YTick',[0:20:100])
set(gca,'XTickLabels','')
axis tight
ylim([0 100])
grid on
%%% WITHOUT GFs
[death2_exp,~,~]=xlsread('data_experimental.xlsx','death2','','basic');
death_exp_avg=death2_exp(1,9:10);
death_exp_std=death2_exp(1,11:12);
death_exp_sem=death_exp_std/sqrt(numcells);
% Experimental
b2(1)=figure;
bp=barwitherr(death_exp_sem, death_exp_avg);
set(bp,'FaceColor',[0.5,0.5,0.5]);
set(gca,'YTick',[0:20:100])
set(gca,'XTickLabels','')
axis tight
ylim([0 100])
grid on
% Simulation
b2(2)=figure;
bp=barwitherr(dsem(1:2,2),d(1:2,2));%set(bp,'FaceColor',[0.5,0.5,0.5]);
set(gca,'YTick',[0:20:100])
set(gca,'XTickLabels','')
axis tight
ylim([0 100])
grid on


% Panel D; p21 mut
% With mutations that do not allow for cell cycle arrest (p21 mut, ATR cant
% upregulate chk1)
c=figure;
errorbar(d(:,1),dsem(:,1),'Color',colororder(1,:),'LineWidth',2);
hold on
errorbar(d(:,3),dsem(:,3),'Color',colororder(2,:),'LineWidth',2);
set(gca,'XTickLabels','')
axis tight
ylim([0 100])

%% Figure 5A
% Apoptosis Treatments

% Stochastic TRAIL dose response - TIMECOURSE
a1=figure;
for k=[1 2 3 4 5] %files/conditions
subplot(5,1,k)
St=load(strcat(matpath,'Apoptosis2S_',num2str(k))); 
cells=St.cells;
color=colororder(k,:);
numcells=length(cells);
    for n=1:numcells %cells
        if ~isempty(cells{n})
        xoutS_all=cells{n}.xoutS_all;
        tout_all=cells{n}.tout_all;
        y=xoutS_all(:,106);
        ind=length(y);
        plot(tout_all(1:ind)/3600,y(1:ind),'Color',color)
        axis tight
        ylim([0 400])
        xlim([0 80])
        hold on
        % Calc ttd
        time=find(xoutS_all(:,106)>100);
        if ~isempty(time)
            ttds(n,k)=tout_all(min(time));
        else ttds(n,k)=Inf;
        end
        end
    end
disp(k)
end
ttds=ttds/3600; %translate to hours
% Cell viability dose response @X hours
[death1_exp,~,~]=xlsread('data_experimental.xlsx','death1','','basic');
numreplicates=3;
death1_exp_avg=death1_exp(2:6,13:15);
death1_exp_std=death1_exp(2:6,19:21);
death1_exp_sem=death1_exp_std/sqrt(numreplicates);
% Bar plots simulation
a2=figure;
for i=1:size(ttds,2)
    d24=sum(ttds(:,i)<=24)/numcells;
    d48=sum(ttds(:,i)<=48)/numcells;
    d72=sum(ttds(:,i)<=72)/numcells;

    d24sem=sqrt((d24*(1-d24))/numcells);
    d48sem=sqrt((d48*(1-d48))/numcells);
    d72sem=sqrt((d72*(1-d72))/numcells);

    d=[d24;d48;d72]*100;
    errY=[d24sem;d48sem;d72sem]*100;
    
    subplot(5,1,i)
    bp=barwitherr(errY, d);
    set(bp(1),'FaceColor',colororder(i,:));

    axis tight
    colormap(colororder(1:2,:))
    set(gca,'XTickLabels','')
    ylim([0 100])
    grid on
    set(gca,'YTick',[0 25 50 75 100])
end
% Bar plots experimental
a3=figure;
for i=1:size(ttds,2)
    d=[[death1_exp_avg(i,1)];[death1_exp_avg(i,2)];[death1_exp_avg(i,3)]];
    errY=[[death1_exp_sem(i,1)];[death1_exp_sem(i,2)];[death1_exp_sem(i,3)]];
    subplot(5,1,i)
    bp=barwitherr(errY, d);
    set(bp(1),'FaceColor',[0.5,0.5,0.5]);
    axis tight
    colormap(colororder(1:2,:))
    set(gca,'XTickLabels','')
    ylim([0 100])
    grid on
    set(gca,'YTick',[0 25 50 75 100])
end

%% Figure 5C
% BIM/BAD Mechanism

dead=[];
for k=1:2
    St=load(strcat(matpath,'Apoptosis2S_Mech2_',num2str(k))); 
    cells=St.cells;
    for i=1:length(cells)
        if cells{i}.xoutS_all(end,106)>100
            dead(i,k)=1;
        else dead(i,k)=0;
        end

    end
end

figure
bar([sum(dead(:,2))/length(dead)*100,sum(dead(:,1))/length(dead)*100])
axis tight
ylim([0 100])

%% Figure 5D
% Histogram of Correlation Between Time to Death and Protein Levels

StS=load(strcat(matpath,'Apoptosis2S_400cells_5.mat'));
cells=StS.cells;
% To get observables
[dataS,dataG]=RunPrep;
for i=1:length(cells)
    y=cells{i}.xoutS_all;
    cells{i}.Obs=GetObservables_matrix(y,dataS.VxPARCDL,dataS.kS(3));
end
ym=[];
y1=[];
for i=1:length(cells)
    x=cells{i}.tout_all/3600;
    ym(i,:)=nanmean(cells{i}.Obs);
    y1(i,:)=cells{i}.Obs(1,:);
    ttd(i)=x(end);
end
inds2use=find(mean(ym)>1E-6); %Remove obs that are essentially zero.
rm=corr(ym(:,inds2use),ttd');
r1=corr(y1(:,inds2use),ttd');
d=figure;
histogram(r1,5);
axis tight
xlim([-1 1])
set(gca,'XTick',[-1 -0.5 0 0.5 1])

%% Figure 5E
% Apoptosis ROC Curves

% Get data
StS=load(strcat(matpath,'Apoptosis2S_400cells_5.mat'));
cells=StS.cells;
% To get observables
[dataS,dataG]=RunPrep;
for i=1:length(cells)
    y=cells{i}.xoutS_all;
    cells{i}.Obs=GetObservables_matrix(y,dataS.VxPARCDL,dataS.kS(3));
end

% To get predictors and responses
Resp=NaN(length(cells),1);
for i=1:length(cells)
    x=cells{i}.tout_all/3600;
    ttds=x(end); %Simulations stop when PARP=cPARP (cell dead)
    if ttds<=40
        Resp(i)=1;
    else Resp(i)=0;
    end
    
    % Single timepoint
    Pred1(i,:)=cells{i}.Obs(1,:);
    
    % Average across all time
    if Resp(i)==1
        PredMean40(i,:)=mean(cells{i}.Obs); %Observables
    else
        PredMean40(i,:)=mean(cells{i}.Obs(1:length(0:0.25:40),:));
    end
    
    % Average across time
    if ttds<8
        PredMean8(i,:)=mean(cells{i}.Obs);
    else
        PredMean8(i,:)=mean(cells{i}.Obs(1:length(0:0.25:8),:));
    end
    
end

% Define training and test sets
tra_set=1:200;
tes_set=201:400;

% Defining species to include
Pred=PredMean40;
inds2remove=zeros(size(Pred,2),1);
inds2remove=mean(Pred)<1E-4;
inds2remove([7,9,10,11:26,27,47:53,58,71])=1;
keeping=find(~inds2remove);

% Lasso regression for PredMean40
[B,FitInfo] = lasso(PredMean40(tra_set,keeping),Resp(tra_set),'Alpha',1,'CV',5);
%lassoPlot(B,FitInfo,'PlotType','CV');
ind=max(find(sum(logical(B))>4));
top_predictors40=keeping(find(B(:,ind)))

% Lasso regression for PredMean8
[B,FitInfo] = lasso(PredMean8(tra_set,keeping),Resp(tra_set),'Alpha',1,'CV',5);
%lassoPlot(B,FitInfo,'PlotType','CV');
ind=max(find(sum(logical(B))>4));
top_predictors8=keeping(find(B(:,ind)))

% Lasso regression for Pred1
[B,FitInfo] = lasso(Pred1(tra_set,keeping),Resp(tra_set),'Alpha',1,'CV',5);
%lassoPlot(B,FitInfo,'PlotType','CV');
ind=max(find(sum(logical(B))>4));
top_predictors1=keeping(find(B(:,ind)))

% Check top >4 predictors versus BIM and Bcl2 as sole predictors
figure;
for i=1:3
    if i==1; pred=PredMean40(:,top_predictors40); end
    if i==2; pred=PredMean8(:,top_predictors8); end
    if i==3; pred=Pred1(:,top_predictors1); end
    
    mdlSVM = fitcsvm(pred(tra_set,:),Resp(tra_set),'Standardize',true);
    [label,score] = predict(mdlSVM,pred(tes_set,:));
    [Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(Resp(tes_set),score(:,2),1);
    plot(Xsvm,Ysvm)
    hold on 
end
% add refline
rl=refline(1,0);
rl.LineStyle='--';
rl.Color='k';
% format plot
set(gca,'XTick',0:0.2:1)
set(gca,'YTick',0:0.2:1)
legend('Top Predictors Mean 40h','Top Predictors Mean 8h','Top Predictors Initial')

figure;
for i=1:3
    if i==1; pred=PredMean40(:,[37,46]); end
    if i==2; pred=PredMean8(:,[37,46]); end
    if i==3; pred=Pred1(:,[37,46]); end
    
    mdlSVM = fitcsvm(pred(tra_set,:),Resp(tra_set),'Standardize',true);
    [label,score] = predict(mdlSVM,pred(tes_set,:));
    [Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(Resp(tes_set),score(:,2),1);
    plot(Xsvm,Ysvm)
    hold on 
end
% add refline
rl=refline(1,0);
rl.LineStyle='--';
rl.Color='k';
% format plot
set(gca,'XTick',0:0.2:1)
set(gca,'YTick',0:0.2:1)
legend('BIM and Bcl2 Mean 40h','BIM and Bcl2 Mean 8h','BIM and Bcl2 Initial')

% Use BIM and Bcl2 as sole predictors for final figure
e=figure;
for i=1:3
    if i==1; pred=PredMean40(:,[37,46]); end
    if i==2; pred=PredMean8(:,[37,46]); end
    if i==3; pred=Pred1(:,[37,46]); end
    
    mdlSVM = fitcsvm(pred(tra_set,:),Resp(tra_set),'Standardize',true);
    [label,score] = predict(mdlSVM,pred(tes_set,:));
    [Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(Resp(tes_set),score(:,2),1);
    AUCsvm
    plot(Xsvm,Ysvm)
    hold on 
end
% add refline
rl=refline(1,0);
rl.LineStyle='--';
rl.Color='k';
% format plot
set(gca,'XTick',0:0.2:1)
set(gca,'YTick',0:0.2:1)

%% Figure 5F

colors=jet(81);
f(1)=figure;
files={'Apoptosis2S_2.mat','Apoptosis2S_1.mat','Apoptosis2S_3.mat','Apoptosis2S_4.mat','Apoptosis2S_5.mat','Apoptosis2S_400cells_5.mat'};
for k=1:length(files)
StS=load(strcat(matpath,files{k}));
cells=StS.cells;
% To get observables
[dataS,dataG]=RunPrep;
for i=1:length(cells)
    y=cells{i}.xoutS_all;
    cells{i}.Obs=GetObservables_matrix(y,dataS.VxPARCDL,dataS.kS(3));
end

for i=1:length(cells)
    t=cells{i}.tout_all;
    figure(f(1))
    pl=plot(cumsum(cells{i}.Obs(:,37)),cumsum(cells{i}.Obs(:,46)),'-','LineWidth',0.2);        
    pl.Color=[colors(round(t(end)/3600)+1,:),0.7];
    hold on
end
end

%% Figure S6A

colors=jet(81);
f(1)=figure;
files={'Apoptosis2S_TRAIL_1ng.mat','Apoptosis2S_TRAIL_01ng.mat'};
for k=1:length(files)
StS=load(strcat(matpath,files{k}));
cells=StS.cells;
% To get observables
[dataS,dataG]=RunPrep;
for i=1:length(cells)
    y=cells{i}.xoutS_all;
    cells{i}.Obs=GetObservables_matrix(y,dataS.VxPARCDL,dataS.kS(3));
end

for i=1:length(cells)
    t=cells{i}.tout_all;
    figure(f(1))
    pl=plot(cumsum(cells{i}.Obs(:,37)),cumsum(cells{i}.Obs(:,46)),'-','LineWidth',0.2);        
    pl.Color=[colors(round(t(end)/3600)+1,:),0.7];
    hold on
end
end

%% Figure S6B
% Details available upon request.

%% Figure 6A

a=figure;
[sphase_exp,~,~]=xlsread('data_experimental.xlsx','sphase','','basic');
sphase_exp_percent=sphase_exp(1:5,3);
sphase_exp_stds=sphase_exp(1:5,4);
sphase_exp_sems=sphase_exp_stds./sqrt([6;2;3;2;2]);% Number of replicates per condition
bp=barwitherr(sphase_exp_sems,sphase_exp_percent,'stacked');
set(bp,'FaceColor',[0.5,0.5,0.5])
axis tight
ylim([0 55])
grid on
set(gca,'XTickLabel',[])  
    
%% Figure 6B   

% Stochastic TRAIL dose response - TIMECOURSE
files={'CellCycle2S_E.mat','CellCycle2S_I.mat','CellCycle2S_E+I.mat','CellCycle2S_E+I+MEKi','CellCycle2S_E+I+AKTi'};
St=load(strcat(matpath,files{1})); 
cells=St.cells;
b1=figure;
numcells=length(cells);
cellsinsphase=zeros(numcells,length(files)); %# of cells
deadcells=zeros(numcells,length(files)); %# of cells
for k=1:length(files)
St=load(strcat(matpath,files{k})); 
cells=St.cells;
ind2plot=[59]; %cycA
ind2define=[50,59,69]; %cyclin E,A,B
subplot(length(files),1,k)
for i=1:length(cells)
    if ~isempty(cells{i})
    y=cells{i}.xoutS_all(:,ind2plot);
    plot(cells{i}.tout_all/3600,y);
    hold on
    axis tight
    ylim([0 50])
    xlim([0 30]);
    set(gca,'XTick',[0:6:30]);
    %See if in S-phase at 24 hrs
    timeind=find(cells{i}.tout_all/3600==24);
    if sum(cells{i}.xoutS_all(timeind,ind2define))>20 %definition of S-phase, when CYCA/cdk2>XnM
        cellsinsphase(i,k)=1;
    end
    if isempty(timeind)
        deadcells(i,k)=1;
    end 
    
    end
end
end

% bar plots for S-phase
% Exp data
[sphase_exp,~,~]=xlsread('data_experimental.xlsx','sphase','','basic');
sphase_exp_percent=sphase_exp(1:5,3);
sphase_exp_stds=sphase_exp(1:5,4);
sphase_exp_sems=sphase_exp_stds./sqrt([6;2;3;2;2]);% Number of replicates per condition

% Plot
numcellsNOTdead=(sum(deadcells)*-1)+numcells;
b2=figure;
sphase_sim=(sum(cellsinsphase)./numcellsNOTdead)*100;
sphase_sim_sem=sqrt(sphase_sim.*(100-sphase_sim)./numcellsNOTdead);
for i=1:length(sphase_sim)
    subplot(length(sphase_sim),1,i)
    d=[sphase_sim(i);sphase_exp_percent(i)];
    errY=[sphase_sim_sem(i);sphase_exp_sems(i)];
    bp=barwitherr(diag(errY),diag(d),'stacked');
    set(bp(2),'FaceColor',[0.5,0.5,0.5])
    axis tight
    ylim([0 65])
    grid on
    set(gca,'XTickLabel',[])  
end

%% Figure 6C

% Cyclin D uWestern data
datasheet='data_experimental.xlsx';
[nums,~,~]=xlsread(datasheet,'uWestern_110216','','basic');
cd_avg=nums(:,17);
cd_std=nums(:,18);
time=[0,5,30,180,360]; %minutes
colors=colororder;%cool(4);
% rows=doses or conditions, columns=timepoints
cd_mat_avg=[cd_avg(1:15),cd_avg(17:17+14),cd_avg(33:33+14),cd_avg(49:49+14)]; 
cd_mat_std=[cd_std(1:15),cd_std(17:17+14),cd_std(33:33+14),cd_std(49:49+14)];
cd_mat_avg_norm=cd_mat_avg./repmat(cd_mat_avg(13,:),15,1);
cd_mat_std_norm=cd_mat_std./repmat(cd_mat_avg(13,:),15,1);
cd_mat_avg_norm=[ones(size(cd_mat_avg_norm,1),1),cd_mat_avg_norm];
cd_mat_std_norm=[repmat(cd_mat_std_norm(13,1),size(cd_mat_std_norm,1),1),cd_mat_std_norm];
egfrows=1:4;
insrows=5:8;
bcrows=9:12;
% Bar plot of Cyclin D at 360 minutes EXPERIMENT
c1=figure;
r=bcrows;
numsets=3;
for k=1:numsets
    if k==1; r=egfrows(end); end
    if k==2; r=insrows(end); end
    if k==3; r=bcrows(end); end
    yavg_ends(k)=cd_mat_avg_norm(r,end);
    ystd_ends(k)=cd_mat_std_norm(r,end);
end
bwe=barwitherr(ystd_ends,yavg_ends);
axis tight
ylim([1 10])
set(gca,'XTick',[])
grid on
set(bwe,'FaceColor',[0.5 0.5 0.5])
% Bar plot for cyclin D at 360 minutes SIMULATION
c2=figure;
for n=1:3 %three conditions
StD=load(strcat(matpath,'GFdoseresponse1D_EI_',num2str(n),'.mat'));
condsD=StD.condsD;
inds=[44    45    46    47    80];
cycdquant(n)=sum(condsD{4}.xoutS_all(end,inds),2); %Assuming the timecourse is 6 hours
end
bar(cycdquant)
axis tight
ylim([1 75])
set(gca,'XTick',[])
grid on

%% Figure 6D

% ERK AND AKT Deterministic
d1=figure;
I=Plot_GetInfo;
tppERK=I.tppERK;
tppAKT=I.tppAKT;
indsE=tppERK;
indsA=tppAKT;
for n=1:3 %three conditions
    StD=load(strcat(matpath,'GFdoseresponse1D_EI_',num2str(n),'.mat'));
    condsD=StD.condsD;
    colors=colororder;%cool(4);
    x=condsD{4}.tout_all/60;
    y=sum(condsD{4}.xoutS_all(:,indsE).*repmat([dataS.VxPARCDL(indsE)/dataS.kS(3)]',length(x),1),2);
    subplot(2,1,1); 
    
    if n==3
        plot(x,y,'--','Color',colors(n,:));
    else
        plot(x,y,'Color',colors(n,:));
    end
    
    xlim([0 180])
    ylim([0 120])
    hold on
    x=condsD{4}.tout_all/60;
    y=sum(condsD{4}.xoutS_all(:,indsA).*repmat([dataS.VxPARCDL(indsA)/dataS.kS(3)]',length(x),1),2);
    subplot(2,1,2); 
    
    if n==3
        plot(x,y,'--','Color',colors(n,:));
    else
        plot(x,y,'Color',colors(n,:));
    end

    hold on
    xlim([0 180])
    ylim([0 35])
end

%% Figure 6E

set(0,'DefaultAxesFontSize',8)

% Heatmap
St=load(strcat(matpath,'CellCycle_erkaktcycd.mat')); 
xoutS_alls_cases=St.xoutS_alls_cases;
xoutS_alls=St.xoutS_alls;

% Calculate metrics
erks=[];
akts=[];
cyds=[];
for i=1:size(xoutS_alls,1)
    for k=1:size(xoutS_alls,2)
        erks(i,k)=mean(xoutS_alls{i,k}(:,718));
        akts(i,k)=mean(xoutS_alls{i,k}(:,697));
        cyds(i,k)=mean(sum(xoutS_alls{i,k}(:,[44    45    46    47    80]),2));
    end
end

% Make heatmap
erks_a=mean(erks,2);
akts_a=mean(akts)';

figure;
h=surf(erks_a(1:7),akts_a(1:7),cyds(1:7,1:7));
h.EdgeColor='none';
h.FaceColor='interp';


% Add cases
hold on
for i=1:length(xoutS_alls_cases)-2
    xc=mean(xoutS_alls_cases{i}(:,718));
    yc=mean(xoutS_alls_cases{i}(:,697));
    d=[xc,yc,mean(sum(xoutS_alls_cases{i}(:,[44,45,46,47,80]),2))];
    plot(d(1),d(2),'.','MarkerSize',40,'Color',colororder(i+3,:))
end
axis tight
colorbar

%% Figure 6F

% To get responses
StS=load(strcat(matpath,'CellCycle2S_E+I_400cells.mat'));
ind2define=[50,59,69]; %cyclin E,A,B
cells=StS.cells;
Resp=NaN(length(cells),1);
for i=1:length(cells)
    x=cells{i}.tout_all/3600;
    y=sum(cells{i}.xoutS_all(:,ind2define),2);
    if sum(y>20); Resp(i)=1; else Resp(i)=0; end;
end
% plot
I=Plot_GetInfo;
tppERK=I.tppERK;
tppAKT=I.tppAKT;
f1=figure;
f2=figure;
for k=1:2;
if k==1; inds=tppERK; end
if k==2; inds=tppAKT; end
figure(f1)
subplot(2,1,k);
Resp_pos=NaN(length(cells{2}.tout_all),length(cells));
Resp_neg=NaN(length(cells{2}.tout_all),length(cells));
[dataS,~]=RunPrep;
ys=[];
for i=1:length(cells)
    x=cells{i}.tout_all/3600;
    y=sum(cells{i}.xoutS_all(:,inds).*repmat([dataS.VxPARCDL(inds)/dataS.kS(3)]',length(x),1),2);
    if Resp(i); Resp_pos(1:length(x),i)=y; pl=plot(x,y,'Color','r','LineWidth',0.1); pl.Color=[1,0,0,1]; else Resp_neg(1:length(x),i)=y; pl=plot(x,y,'Color','b','LineWidth',0.1); pl.Color=[0,0,1,1]; end
    hold on
    ys(i)=y(1);
    yhist(i)=trapz(y);
end
% plot averages
plot(x, nanmean(Resp_neg,2),'--k','LineWidth',1)
hold on
plot(x, nanmean(Resp_pos,2),'-k','LineWidth',1)
if k==1; ylim([0 100]); end
if k==2; ylim([0 30]); end
% bar plots
figure(f2)
subplot(2,1,k)
x=ys(find(Resp));
y=ys(find(~Resp));
barwitherr([std(x),std(y)],[mean(x),mean(y)])
set(gca,'XTick',[])
axis tight
[h,p,ci,stats]=ttest2(ys(find(Resp)),ys(find(~Resp)));
disp(p)
if k==1; ylim([0 35]); end
if k==2; ylim([0 0.75]); end
end

%% Figure 6G
% ROC curves

% To get observables and predictors
StS=load(strcat(matpath,'CellCycle2S_E+I_400cells.mat'));
cells=StS.cells;
% get responses
Resp=NaN(length(cells),1);
for i=1:length(cells)
    x=cells{i}.tout_all/3600;
    y=cells{i}.xoutS_all(:,59); %cell cycle
    if sum(y>10); Resp(i)=1; else Resp(i)=0; end
end
% get predictions
Pred=[];
[dataS,dataG]=RunPrep;
for i=1:length(cells)
    y=cells{i}.xoutS_all;
    cells{i}.Obs=GetObservables_matrix(y,dataS.VxPARCDL,dataS.kS(3));
    Pred(i,:)=cells{i}.Obs(1,:); %predictors
end

% Lasso regression and ROC curves
% Defining important species
inds2remove=zeros(size(Pred,2),1);
inds2remove=mean(Pred)<1E-4;
inds2remove([7,9,10,11:26,27,47:53,58,71])=1; %Observables: cell cycle proteins
keeping=find(~inds2remove);

% Defining training and test sets
tra_set=1:200;
tes_set=201:400;

% Lasso regression
[B,FitInfo] = lasso(Pred(tra_set,keeping),Resp(tra_set),'Alpha',1,'CV',5);
%lassoPlot(B,FitInfo,'PlotType','CV');
ind=max(find(sum(logical(B))>1));
top_predictors=keeping(find(B(:,ind)))

% ROC curves
g=figure;
% BRaf, CRaf
keeping_new=[75,76];
pred=Pred(:,keeping_new);
mdlSVM = fitcsvm(pred(tra_set,:),Resp(tra_set),'Standardize',true);
[label,score] = predict(mdlSVM,pred(tes_set,:));
[Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(Resp(tes_set),score(:,2),1);
AUCsvm
plot(Xsvm,Ysvm,'LineWidth',1.2)
% BRaf
keeping_new=[76];
pred=Pred(:,keeping_new);
mdlSVM = fitcsvm(pred(tra_set,:),Resp(tra_set),'Standardize',true);
[label,score] = predict(mdlSVM,pred(tes_set,:));
[Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(Resp(tes_set),score(:,2),1);
AUCsvm
hold on; plot(Xsvm,Ysvm,'LineWidth',1.2)
% CRaf
keeping_new=[75];
pred=Pred(:,keeping_new);
mdlSVM = fitcsvm(pred(tra_set,:),Resp(tra_set),'Standardize',true);
[label,score] = predict(mdlSVM,pred(tes_set,:));
[Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(Resp(tes_set),score(:,2),1);
AUCsvm
hold on; plot(Xsvm,Ysvm,'LineWidth',1.2)
% ERK and AKT
keeping_new=[95,86]; %ERK and AKT
pred=Pred(:,keeping_new);
mdlSVM = fitcsvm(pred(tra_set,:),Resp(tra_set),'Standardize',true);
[label,score] = predict(mdlSVM,pred(tes_set,:));
[Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(Resp(tes_set),score(:,2),1);
AUCsvm
hold on; plot(Xsvm,Ysvm,'LineWidth',1.2)
% add refline
rl=refline(1,0);
rl.LineStyle='--';
rl.Color='k';
% format plot
set(gca,'XTick',0:0.2:1)
set(gca,'YTick',0:0.2:1)

%% Figure S7
% AKT levels investigation

filenames={'CellCycle_EI_25.mat','CellCycle_E_25.mat','CellCycle_EAKTi1_25.mat','CellCycle_EAKTi2_25.mat'};
I=Plot_GetInfo;
tppERK=I.tppERK;
tppAKT=I.tppAKT;
yERK=[];
yAKT=[];
Resp=[];
[dataS,~]=RunPrep;

for k=1:length(filenames)
    StS=load(strcat(matpath,filenames{k}));
    ind2define=[50,59,69]; %cyclin E,A,B
    cells=StS.cells;
    for i=1:length(cells)
        x=cells{i}.tout_all/3600;
        y=sum(cells{i}.xoutS_all(:,ind2define),2);
        if sum(y>20); Resp(i,k)=1; else Resp(i,k)=0; end;
    end
    for i=1:length(cells)
        x=cells{i}.tout_all/3600;
        yERK(i,k)=mean(sum(cells{i}.xoutS_all(:,tppERK).*repmat([dataS.VxPARCDL(tppERK)/dataS.kS(3)]',length(x),1),2))';
        yAKT(i,k)=mean(sum(cells{i}.xoutS_all(:,tppAKT).*repmat([dataS.VxPARCDL(tppAKT)/dataS.kS(3)]',length(x),1),2))';
    end
end

% number of cells
%#3
c=figure;
bar((sum(Resp)./length(Resp))*100,'k')
axis tight
%ylim([0 70])
set(gca,'XTick',[])

% total ppERK and ppAKT levels
% #2
b=figure;
bar(mean(yERK),'k')
axis tight
%ylim([0 30])
set(gca,'XTick',[])

% #1
a=figure;
bar(mean(yAKT),'k')
axis tight
%ylim([0 10])
set(gca,'XTick',[])


% for cycling cells
for i=1:length(filenames)
    yERK_C(i)=mean(yERK(logical(Resp(:,i)),i));
    yAKT_C(i)=mean(yAKT(logical(Resp(:,i)),i));
end
%#4
d=figure;
bar(yERK_C,'k')
axis tight
ylim([0 50])
set(gca,'XTick',[])

%% Figure 7A
% Amplification/Deletion Test

% Define gene names for labels
gene_names={'MYC'
'CTNNB1'
'PTEN'
'AKT1'
'AKT2'
'PDPK1'
'RICTOR'
'MTOR'
'GSK3B'
'TSC1'
'TSC2'
'MAPK1'
'MAPK3'
'RHEB'
'RPTOR'
'RPS6KB1'
'RPS6KB2'
'EIF4EBP1'
'SOS1'
'EIF4E'
'IRS1'
'IRS2'
'INSR'
'EGFR'
'ERBB2'
'ERBB3'
'ERBB4'
'SPRY2'
'CBL'
'PLCG1'
'PLCG2'
'PIK3CA'
'PIK3CB'
'PIK3CG'
'PIK3CD'
'PIK3R1'
'PIK3R2'
'PIK3R3'
'PIK3R4'
'RASGRP1'
'RASGRP3'
'NRAS'
'KRAS'
'HRAS'
'NF1'
'RAF1'
'BRAF'
'MAP2K1'
'MAP2K2'
'DUSP6'
'RPS6KA1'
'RPS6KA2'
'RPS6KA3'
'RPS6KA4'
'DUSP1'
'FOS'
'JUN'
'Wild-Type'};
    

% Get data
filenames=dir(strcat(matpath,'mut_*.mat'));

% Plot
numcps=[];

for i=1:length(filenames)
    St=load(strcat(matpath,filenames(i).name)); 
    cells=St.cells;
    numcp=0;
    for k=1:length(cells)  
    indsprolif=[50,59,69]; %cyclin E,A,B
    timeind=find(cells{k}.tout_all/3600==24);
    if sum(cells{k}.xoutS_all(timeind,indsprolif))>20
        numcp=numcp+1;
    end
    end
    numcps(i)=(numcp/length(cells))*100;
end

[numcps_sort,inds_sort]=sort(numcps,'ascend');

set(0,'DefaultAxesFontSize',8)

a=figure; h=barh(diag(numcps_sort),'k','stacked');
hold on
set(h(23),'FaceColor',[0.7 0.7 0.7])
set(gca,'YTick',1:length(numcps))
set(gca,'YTickLabels',gene_names(inds_sort))
axis tight

%% Figure 7B-D
% MEK Mechanism Investigation

%Get data and plot
I=Plot_GetInfo;
tppERK=I.tppERK;
num1=[];
num2=[];
for l=1:3
if l==1; St=load(strcat(matpath,'mut_91.mat')); end
if l==2; St=load(strcat(matpath,'mut_none.mat')); end
if l==3; St=load(strcat(matpath,'MEKmechanism.mat')); end
cells=St.cells;
ind=tppERK;
totaltime=0:900:3600*24;
numinterest=nan(length(totaltime),length(cells));
for i=1:length(cells)
    x=sum(cells{i}.xoutS_all(:,ind),2);
    numinterest(1:length(x),i)=x;
    num1(l,i)=sum(cells{i}.xoutS_all(1,[666,667:670,759,760,765]));
    num2(l,i)=sum(cells{i}.xoutS_all(1,[734,735,761,762,766,767])); 
end

% calculate percent proliferating
numcp=0;
for k=1:length(cells)  
indsprolif=[50,59,69]; %cyclin E,A,B
timeind=find(cells{k}.tout_all/3600==24);
if sum(cells{k}.xoutS_all(timeind,indsprolif))>20
    numcp=numcp+1;
end
end
numcps(l)=(numcp/length(cells))*100;
meansERK(l)=mean(nanmean(numinterest,2));
stdsERK(l)=mean(nanstd(numinterest,0,2));
hold on
end

axis tight
xlim([0 10])

b1=figure; bar(numcps,'k')
set(gca,'XTick',[]);
ylim([0 70])

b2=figure; barwitherr(std(num2,0,2),mean(num2,2),'k')
set(gca,'XTick',[]);

b3=figure; bar(meansERK,'k')
set(gca,'XTick',[]);
ylim([30 37])

set(0,'DefaultAxesFontSize',8)

%% Figure 8
% u87 application

% Apoptosis Treatments
txt{1}={'Apoptosis2S_10A_1.mat',...
      'Apoptosis2S_10A_2.mat',...
      'Apoptosis2S_10A_3.mat',...
      'Apoptosis2S_10A_4.mat'};
txt{2}={'Apoptosis2S_u87_1.mat',...
      'Apoptosis2S_u87_2.mat',...
      'Apoptosis2S_u87_3.mat',...
      'Apoptosis2S_u87_4.mat'};
  
figure;
d48=[];
d48sem=[];
ttds=[];
for s=1:length(txt) 
for k=1:4 %files/conditions
St=load([matpath,txt{s}{k}]); 
cells=St.cells;
numcells(s,k)=length(cells);
    for n=1:numcells(s,k) %cells
        if ~isempty(cells{n})
        xoutS_all=cells{n}.xoutS_all;
        tout_all=cells{n}.tout_all;
        y=xoutS_all(:,106);
        % Calc ttd
        time=find(xoutS_all(:,106)>100);
        if ~isempty(time)
            ttds(n,k)=tout_all(min(time));
        else ttds(n,k)=Inf;
        end
        end
    end
disp(k)
end
ttds=ttds/3600; %translate to hours
d48(s,:)=sum(ttds<=48)./length(ttds);
d48sem(s,:)=sqrt((d48(s,:).*(1-d48(s,:)))./numcells(s,:));
end
d=d48*100;
errY=d48sem*100;
a=figure;
bpa=barwitherr(errY', d');
set(bpa(1),'FaceColor',[0,0,1]);
set(bpa(2),'FaceColor',[1,0,0]);
set(gca,'XTickLabels','')
axis tight
ylim([0 100])
grid on


% Bar plots experimental u87
figure;
[deathu87_exp,~,~]=xlsread('data_experimental.xlsx','u87data','','basic');
numreplicates=3;
death1_exp_avg=deathu87_exp(1,:);
death1_exp_std=deathu87_exp(2,:);
death1_exp_sem=death1_exp_std/sqrt(numreplicates);
b=figure;
bpb=barwitherr(death1_exp_sem, death1_exp_avg);
set(bpb(1),'FaceColor',[0.5,0.5,0.5]);
axis tight
colormap(colororder(1:2,:))
set(gca,'XTickLabels','')
ylim([0 100])
grid on


% mechanism comparison
S10A=load([matpath,'Apoptosis2S_10A_3.mat']);
SU87=load([matpath,'Apoptosis2S_u87_3.mat']);
ind=[719,720,754,768];
c=figure;
d=figure;
for i=1:length(S10A.cells)
    x1=S10A.cells{i}.tout_all/3600;
    y1=S10A.cells{i}.xoutG_all(:,337);
    y11=sum(S10A.cells{i}.xoutS_all(:,ind),2);
    figure(c)
    plot(x1,y1,'b')
    hold on
    figure(d)
    plot(x1,y11,'b')
    hold on
    
    x2=SU87.cells{i}.tout_all/3600;
    y2=SU87.cells{i}.xoutG_all(:,337);
    y22=sum(SU87.cells{i}.xoutS_all(:,ind),2);
    figure(c)
    plot(x2,y2,'r') 
    figure(d)
    plot(x2,y22,'r') 
end
figure(c); axis tight; ylim([0 300])
figure(d); axis tight

end %FUNCTION
