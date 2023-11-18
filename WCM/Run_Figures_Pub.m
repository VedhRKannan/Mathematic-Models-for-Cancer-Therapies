function Run_Figures_Pub

%% UBIQUITOUS CODE
%Path where matfiles will be saved
matpath='matfiles/';
mkdir(matpath)
%Number of species
NumSpecies=775;


%% Figure 2B
flagD=0;
th=24;
STIM=zeros(NumSpecies,1);
[tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],[],[]);
tout_all=tout_all;
xoutG_all=xoutG_all;
xoutS_all=xoutS_all;

%Save out data
txt=strcat(matpath,'STE_1D','.mat');
save(txt,'-v7.3','xoutG_all','xoutS_all','tout_all');

%% Figure 2C

numcells=100;
filename=[matpath,'RandomPopCells_100cells.mat'];
RunRandomPopCells(numcells,filename)

%% Figure S2B
% Test time step; phenotypic test

timesteps=[0.1,0.5,1,5,10]*60;
numcells=30;
cells=[];
flagD=0;
STIM=zeros(775,1);
STIM(156:162)=[10,10,10,10,10,10,1721]; %ligs; full growth stimulus
th=50;
for k=1:length(timesteps)
    for i=1:numcells
        [dataS,~]=RunPrep;
        dataS.ts=timesteps(k);
        [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],dataS,[]);
        
        inds=1:30:length(tout_all);
        cells{i,k}.tout_all=tout_all(inds);
        cells{i,k}.xoutG_all=xoutG_all(inds,:);
        cells{i,k}.xoutS_all=xoutS_all(inds,:);
    end
end

%Save out data
txt=strcat(matpath,'TestTimeStep','.mat');
save(txt,'-v7.3','cells');

%% Figure S2C
% Effect of EIF4E on Translation Rate

%No simulations

%% Figure S2D
% Ribosome levels doubling in ~20 hours

STIM=zeros(775,1);
STIM(156:162)=[100,100,100,100,100,100,100]; %Full growth simulus
[tout_all,~,xoutS_all]=RunModel(1,24,STIM,[],[],[],[]);

txt=strcat(matpath,'Ribosome_dynamics','.mat');
save(txt,'-v7.3','tout_all','xoutS_all');

%% Figure S2E
% Global correlation between all proteins stimulated with growth factors

% Need to run Fig 6F-D to get this data.

%% Figure 3A
% Ligand-Receptor Cooperativity Receptors All Highly Expressed

% In addition to simulations from Figure 3A:

[dataS,dataG]=RunPrep;
dataS_up=dataS;

kA77=0;%apoptosis 1
kA87=0;%apoptosis 2
Rt=dlmread('initialized/i_Rt.txt');
EIF4Efree=dlmread('initialized/i_EIF4Efree.txt');
kDDbasal=dlmread('initialized/i_kDDbasal.txt');
% modifying data.S structure
dataS_up.kS(1)=Rt;
dataS_up.kS(2)=EIF4Efree;
dataS_up.kS(703)=kA77;
dataS_up.kS(713)=kA87;
dataS_up.kS(450)=kDDbasal;
% Turn off synthesis and degradation reactions    
dataS_up.kS(12:14)=0; %ribosome synthesis
dataS_up.kS(17:157)=0; %kTL
dataS_up.kS(299:400)=0; %kTLCd
dataS_up.kS(2050:2655)=0; %kXd
dataS_up.kS([1044,1045])=0; %kP(173:174)
dataS_up.kS(2656:2658)=0; %kE(1:3)    
% Turn off anything beyond receptors
dataS_up.kS(1062:1177)=0; %kRP1,kRP2,kRP3,kRP4

options=odeset('RelTol',1e-6);

xoutS=dlmread('initialized/i_xoutF.txt'); xoutS=xoutS(end,:);

Obs = GetObservables(xoutS,dataS.VxPARCDL,dataS.kS(3));
Ramount=Obs(54);
xoutS(163:188)=0;


doses=[0.001,0.01,0.1,1,10,100,1000];
indsligs=156:162;
th=1;

for i=1:length(indsligs) %ligands

if indsligs(i)==156; xoutS([163])=Ramount; end
if indsligs(i)==157; xoutS([167,168])=Ramount; end
if indsligs(i)==158; xoutS([171])=Ramount; end
if indsligs(i)==159; xoutS([172])=Ramount; end
if indsligs(i)==160; xoutS([173])=Ramount; end
if indsligs(i)==161; xoutS([174])=Ramount; end
if indsligs(i)==162; xoutS([175])=Ramount; end

tstart=0;
tend=3600*1;
tspan=tstart:60:tend;
[tout_all,xoutS_all]=ode15s(@createODEs,tspan,xoutS,options,dataS_up);

xoutS_new=xoutS_all(end,:);

    STIMS=zeros(775,length(doses));
    STIMS(indsligs(i),:)=doses; %insert doses of that particular ligand
    for k=1:size(STIMS,2) %doses
        STIM=STIMS(:,k);
        
        STIMtoch=STIM(1:end-1);
        xout_update=xoutS_new;
        xout_update(logical(STIMtoch))=STIMtoch(logical(STIMtoch));
        
        tstart=0;
        tend=3600*1;
        tspan=tstart:60:tend;
        [tout_all,xoutS_all]=ode15s(@createODEs,tspan,xout_update,options,dataS_up);

        C.tout_all=tout_all;
        C.xoutS_all=xoutS_all;
        conds{i,k}=C; %rows-ligands; columns-doses
    end
    disp(i)
end

txt=strcat(matpath,'Cooperativity1D_sameR','.mat');
save(txt,'-v7.3','conds');

%% Figure 3B-D
% pERK, pAKT, p4EBP1 Dynamics

condsD=[];
for k=1:3
STIMS=zeros(775,4);
indegf=156;
indins=162;
if k==1
STIMS(indegf,1)=0.01;
STIMS(indegf,2)=0.1;
STIMS(indegf,3)=1;
STIMS(indegf,4)=10;
end
if k==2
STIMS(indins,1)=0.17;
STIMS(indins,2)=1.7;
STIMS(indins,3)=17;
STIMS(indins,4)=1721;
end
if k==3
STIMS([indegf,indins],1)=[0.01,0.17];
STIMS([indegf,indins],2)=[0.01,1721];
STIMS([indegf,indins],3)=[10,0.17];
STIMS([indegf,indins],4)=[10,1721];
end

th=6;
flagD=1;
condsD=[];
for i=1:size(STIMS,2)
    STIM=STIMS(:,i);
    [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],[],[]);
    D.tout_all=tout_all;
    D.xoutG_all=xoutG_all;
    D.xoutS_all=xoutS_all;
    condsD{i}=D;    
    disp(i)

end
txt=strcat(matpath,'GFdoseresponse1D_EI_',num2str(k),'.mat');
save(txt,'-v7.3','condsD');
end

%% Figure 3E
% DNA Damage Deterministic Simulations

% Simulate different kinds of DNA damage
% SSB
th=24;
flagD=1;
STIM=zeros(775,1);
STIM(32)=25;
[dataS,dataG]=RunPrep; dataS_up=dataS;
m=[0,0.15];
for i=1:length(m)
    dataS_up.kS(446:448)=dataS.kS(446:448)*m(i);
    [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],dataS_up,[]);
    C.tout_all=tout_all;
    C.xoutS_all=xoutS_all;
    condsD{i}=C;
end
txt=strcat(matpath,'DNADamage1D_1','.mat');
save(txt,'-v7.3','condsD');

% DSB
th=24;
flagD=1;
STIM=zeros(775,1);
STIM(31)=25;
[dataS,dataG]=RunPrep; dataS_up=dataS;
m=[0,5];
for i=1:length(m)
    dataS_up.kS(446:448)=dataS.kS(446:448)*m(i);
    [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],dataS_up,[]);
    C.tout_all=tout_all;
    C.xoutS_all=xoutS_all;
    condsD{i}=C;
end
txt=strcat(matpath,'DNADamage1D_2','.mat');
save(txt,'-v7.3','condsD');

% Both
th=24;
flagD=1;
STIM=zeros(775,1);
STIM(31:32)=25;
[dataS,dataG]=RunPrep; dataS_up=dataS;
m=[0,5];
for i=1:length(m)
    dataS_up.kS(446:448)=dataS.kS(446:448)*m(i);
    [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],dataS_up,[]);
    C.tout_all=tout_all;
    C.xoutS_all=xoutS_all;
    condsD{i}=C;
end
txt=strcat(matpath,'DNADamage1D_3','.mat');
save(txt,'-v7.3','condsD');

%% Figure 3F
% Stochastic DNA Damage Unit Testing

damages=[1 4 10 25];
[dataS,dataG]=RunPrep; dataS_up=dataS;
dataS_up.kS([699,701])=dataS.kS([699,701])*0; %PUMA and NOXA binding BCL2- to prevent death from happening.
th=30;
flagD=0;
for i=1:length(damages)
    STIM=zeros(775,1);
    STIM(31)=damages(i);
    for k=1:20% cells
    [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],dataS_up,[]);
    inds=1:10:length(tout_all);
    C.tout_all=tout_all(inds);
    C.xoutS_all=xoutS_all(inds,:);
    cells{i,k}=C;
    end
end

txt=strcat(matpath,'DNADamage1S','.mat');
save(txt,'-v7.3','cells');

%% Figure 3G
% TRAIL dose response Time Course

traildoses=[0.000385 0.001925 0.00385 0.01925 0.0385 0.1925 0.385 1.9250 3.85 19.25 38.5];% 1,10,100,1000ng/mL   11.5];(300ng/mL)
STIM=zeros(775,1);
th=200;
flagD=1;
for i=1:length(traildoses)
    STIM(84)=traildoses(i);
    [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],[],[]);
    C.tout_all=tout_all;
    C.xoutG_all=xoutG_all;
    C.xoutS_all=xoutS_all;
    condsD{i}=C;
    disp(i)
end

txt=strcat(matpath,'Apoptosis1D_dr','.mat');
save(txt,'-v7.3','condsD');

%% Figure 3H
% TRAIL dose response Percent Cell Death

% Make serum-fed population of random cells
STIM=zeros(775,1);
STIM(156:162)=[3.3,100,100,100,100,100,1721]; %GM sim
filename=[matpath,'RandomPopCells_GM.mat'];
RunRandomPopCells(20,filename,STIM)

% Do TRAIL dose response
traildoses=[0.0385 0.1925 0.385 1.9250 3.85 19.25 38.5];% 1,10,100,1000ng/mL   11.5];(300ng/mL)
STIM=zeros(775,1);
STIM(156:162)=[3.3,100,100,100,100,100,1721]; %GM sim
th=5;
flagD=0;
St=load(filename);
cells0=St.cells0;
cells0=cells0(1:10);
for i=1:length(traildoses)
    for k=1:length(cells0)
        xoutG=cells0{k}.xoutG_all;
        xoutS=cells0{k}.xoutS_all;
        STIM(84)=traildoses(i);
        [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,xoutS,xoutG',[],[]);
        C.tout_all=tout_all;
        C.xoutG_all=xoutG_all;
        C.xoutS_all=xoutS_all;    
        condsS{i,k}=C;
    end
    disp(i)
end
txt=strcat(matpath,'Apoptosis1S_dr','.mat');
save(txt,'-v7.3','condsS');

%% Figure 3I
% Effect of ppERK/ppAKT and PUMA/NOXA on time-to-death

conds=[];
% Add trail to system
STIM=zeros(775,1);
STIM(84)=0.005; %TRAIL, low dose
% No other stim (baseline)
[tout_all1,xoutG_all1,xoutS_all1]=RunModel(1,100,STIM,[],[],[],[]);
C.tout_all=tout_all1;
C.xoutG_all=xoutG_all1;
C.xoutS_all=xoutS_all1;
conds{1}=C;
% Increase ppERK and ppAKT
[dataS,dataG]=RunPrep;
dataS_up=dataS;
%dataS_up.kS([1046,1047,1049,1050,988,992])=0; %pptases
dataS_up.kS([910,913])=dataS_up.kS([910,913])*10;%ppERK
dataS_up.kS([991,995])=dataS_up.kS([991,995])*10;%ppAKT
[tout_all2,xoutG_all2,xoutS_all2]=RunModel(1,100,STIM,[],[],dataS_up,[]);
C.tout_all=tout_all2;
C.xoutG_all=xoutG_all2;
C.xoutS_all=xoutS_all2;
conds{2}=C;
% Add puma and noxa
[dataS,dataG]=RunPrep;
dataS_up=dataS;
xoutG=dataG.x0gm_mpc_D;
xoutG(335:336)=xoutG(335:336)*20; %Add mRNA for PUMA and NOXA
%dataS_up.mExp_nM(53:54)=dataS_up.mExp_nM(53:54)*1000; %Add mRNA for PUMA and NOXA
[tout_all3,xoutG_all3,xoutS_all3]=RunModel(1,100,STIM,[],xoutG,dataS_up,[]);
C.tout_all=tout_all3;
C.xoutG_all=xoutG_all3;
C.xoutS_all=xoutS_all3;
conds{3}=C;


txt=strcat(matpath,'Apoptosis1D','.mat');
save(txt,'-v7.3','conds');

%% Figure 3J
% Dynamics of Cell Cycle Proteins

% Prep
STIM=zeros(NumSpecies,1);
th=96;
flagD=1;
[dataS,dataG]=RunPrep;
dataS_up=dataS;
dataG_up=dataG;

% turn off synthesis and degradation of cyclinD mRNA
pathi='initialized/'; 
kTCleak=dlmread(strcat(pathi,'i_kTCleakF.txt'));
kTCmaxs=dlmread(strcat(pathi,'i_kTCmaxsF.txt'));
dataG_up.kTCleak=kTCleak;
dataG_up.kTCmaxs=kTCmaxs;

dataG_up.kTCleak(10:12)=0; %turn off transcription
dataG_up.kTCmaxs(10:12)=0; %turn off transcription
dataG_up.kTCd(10:12)=0; %turn off mRNA degradation
xoutG=dataG.x0gm_mpc_D;
xoutG_up=xoutG;

m=[1;10;60];
conds=[];

for i=1:length(m)
    xoutG_up(292:294)=xoutG(292:294)*m(i);
    [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],xoutG_up,dataS_up,dataG_up,dataG_up.kTCleak,dataG_up.kTCmaxs);
    C.tout_all=tout_all;
    C.xoutG_all=xoutG_all;
    C.xoutS_all=xoutS_all;
    conds{i}=C;
    disp(i)
end

txt=strcat(matpath,'CellCycle1D','.mat');
save(txt,'-v7.3','conds');

%% Figure S3A
% Ligand-Receptor Cooperativity MCF10A Context

% MCF10A levels    
[dataS,dataG]=RunPrep;
dataS_up=dataS;

kA77=0;
kA87=0;
Rt=dlmread('initialized/i_Rt.txt');
EIF4Efree=dlmread('initialized/i_EIF4Efree.txt');
kDDbasal=dlmread('initialized/i_kDDbasal.txt');
% modifying data.S structure
dataS_up.kS(1)=Rt;
dataS_up.kS(2)=EIF4Efree;
dataS_up.kS(703)=kA77;
dataS_up.kS(713)=kA87;
dataS_up.kS(450)=kDDbasal;
% Turn off synthesis and degradation reactions    
dataS_up.kS(12:14)=0; %ribosome synthesis
dataS_up.kS(17:157)=0; %kTL
dataS_up.kS(299:400)=0; %kTLCd
dataS_up.kS(2050:2655)=0; %kXd
dataS_up.kS([1044,1045])=0; %kP(173:174)
dataS_up.kS(2656:2658)=0; %kE(1:3)    
% Turn off anything beyond receptors
dataS_up.kS(1062:1177)=0; %kRP1,kRP2,kRP3,kRP4

xoutS=dlmread('initialized/i_xoutF.txt'); xoutS=xoutS(end,:);

doses=[0.001,0.01,0.1,1,10,100,1000];
indsligs=156:162;
flagD=1; %Deterministic
th=1;

options=odeset('RelTol',1e-6);

for i=1:length(indsligs) %ligands
    STIMS=zeros(775,length(doses));
    STIMS(indsligs(i),:)=doses; %insert doses of that particular ligand
    for k=1:size(STIMS,2) %doses
        STIM=STIMS(:,k);
        
        STIMtoch=STIM(1:end-1);
        xout_update=xoutS;
        xout_update(logical(STIMtoch))=STIMtoch(logical(STIMtoch));
        
        tstart=0;
        tend=3600*th;
        tspan=[tstart:60:tend];
        [tout_all,xoutS_all]=ode15s(@createODEs,tspan,xout_update,options,dataS_up);

        C.tout_all=tout_all;
        C.xoutS_all=xoutS_all;
        conds{i,k}=C; %rows-ligands; columns-doses
    end
    disp(i)
end

txt=strcat(matpath,'Cooperativity1D','.mat');
save(txt,'-v7.3','conds');

%% Figure S3B
% Phosphorylated receptor dynamics

STIM=zeros(775,1);
STIM(156:162)=[10,0,0,0,0,0,0];
[tout_all,~,xoutS_all]=RunModel(1,2,STIM,[],[],[],[]);

txt=strcat(matpath,'Phos_Receptor_Dynamics','.mat');
save(txt,'-v7.3','tout_all','xoutS_all');

%% Figure S3C
% ppERK and ppAKT Dynamic Dose Response

ligands=156:162;
doses=[0.001,0.01,0.1,1,10,100,1000];

th=3;
flagD=1;
cells=[];

for k=1:length(ligands)
for i=1:length(doses)
    STIM=zeros(775,1);
    STIM(ligands(k))=doses(i);
    [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],[],[]);
    C.tout_all=tout_all;
    C.xoutG_all=xoutG_all;
    C.xoutS_all=xoutS_all;
    cells{k,i}=C;
    disp(i)
end
end

txt=strcat(matpath,'GFDoseResponse','.mat');
save(txt,'-v7.3','cells');

%% Figure S3E
% Stochastic DNA Damage Unit Testing (Time courses)

% Use simulations from Figure 3F

%% Figure S5

% Get random population of cells if necessary
%RunRandomPopCells(30,'matfiles/RandomPopCellsMCF10A_Hetero.mat');

th=2;

% EGF only
STIM=zeros(775,1);
STIM(156:162)=[10,0,0,0,0,0,1721]; %nM

matpath='matfiles/';
[~,dataG]=RunPrep;

flagD=0;
St=load([matpath,'RandomPopCellsMCF10A_Hetero.mat']);
cells0=St.cells0;

for i=1:length(cells0)
    xoutG=cells0{i}.xoutG_all;
    xoutS=cells0{i}.xoutS_all;
    [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,xoutS,xoutG',[],[]);
    C.tout_all=tout_all;
    C.xoutG_all=xoutG_all;
    C.xoutS_all=xoutS_all;
    cells{i,1}=C;
    disp(strcat('cell number =',num2str(i)))
end

filename='CellCycle_E+I_HeteroSims';
txt=strcat(matpath,filename,'.mat');
save(txt,'-v7.3','cells');

%% Figure 4
% DNA Damage Simulations

% Create Random population of cells (if needed)
RunRandomPopCells(25,[matpath,'RandomPopCells.mat']) %for 25 cells

NumSpecies=775;
% WITH GFs
STIM1=zeros(NumSpecies,1);
STIM1(156:162)=[3.3,0,0,0,0,0,1721]; %ligs
th1=72;
th2=72;
STIM2=zeros(NumSpecies,1);
STIM2(775)=100000;
STIM2(156:162)=[3.3,0,0,0,0,0,1721]; %ligs
Run_SimCells2(STIM1,STIM2,th1,th2,'DNADamageSerum')

% WITHOUT GFs
th=72;
STIM=zeros(NumSpecies,1);
STIM(775)=100000;
STIM(156:162)=[0,0,0,0,0,0,0]; %ligs
Run_SimCells(STIM,th,'DNADamageNoSerum')

% WITH GFs, p21/Chk1 mutation
[dataS,~]=RunPrep;
dataS.kS([518,505,507,509])=0; %p21 binding cyclins to zero.
dataS.kS(580)=0; %ATRac activating chk1
STIM1=zeros(NumSpecies,1);
STIM1(156:162)=[3.3,0,0,0,0,0,1721]; %ligs
th1=72;
th2=72;
STIM2=zeros(NumSpecies,1);
STIM2(775)=100000;
STIM2(156:162)=[3.3,0,0,0,0,0,1721]; %ligs
Run_SimCells2(STIM1,STIM2,th1,th2,'DNADamageSerum_p21mut',dataS)

%% Figure 5A
% Apoptosis Treatments

% Create Random population of cells (if needed)
RunRandomPopCells(25,[matpath,'RandomPopCells.mat']) %for 25 cells
%RunRandomPopCells(50,[matpath,'RandomPopCells.mat']) %for 50 cells

% Effect of inhibitors and GFs, TRAIL dose response in SS state.
m=1; %NO STIM
STIM=zeros(NumSpecies,1);
STIM(156:162)=[0,0,0,0,0,0,0]; %ligs
STIM(769)=0; %MEKi
STIM(771)=0; %AKTi
th=80;
filename=strcat('Apoptosis2S_',num2str(m));
Run_SimCells(STIM,th,filename)

m=2;
STIM=zeros(NumSpecies,1);
STIM(156:162)=[3.3,0,0,0,0,0,1721]; %ligs
STIM(769)=0; %MEKi
STIM(771)=0; %AKTi
th=80;
filename=strcat('Apoptosis2S_',num2str(m));
Run_SimCells(STIM,th,filename)

m=3;
STIM=zeros(NumSpecies,1);
STIM(156:162)=[3.3,0,0,0,0,0,1721]; %ligs
STIM(769)=10000; %MEKi
STIM(771)=0; %AKTi
th=80;
filename=strcat('Apoptosis2S_',num2str(m));
Run_SimCells(STIM,th,filename)

m=4;
STIM=zeros(NumSpecies,1);
STIM(156:162)=[3.3,0,0,0,0,0,1721]; %ligs
STIM(769)=0; %MEKi
STIM(771)=10000; %AKTi
th=80;
filename=strcat('Apoptosis2S_',num2str(m));
Run_SimCells(STIM,th,filename)

m=5;
STIM=zeros(NumSpecies,1);
STIM(156:162)=[3.3,0,0,0,0,0,1721]; %ligs
STIM(769)=10000; %MEKi
STIM(771)=10000; %AKTi
th=80;
filename=strcat('Apoptosis2S_',num2str(m));
Run_SimCells(STIM,th,filename)

%% Figure 5C
% BIM/BAD Mechanism

for m=1:2
    
STIM=zeros(NumSpecies,1);

%BAD
if m==1
[dataS,dataG]=RunPrep;
dataS.kS(2052)=0; %ppERK can't bind BAD (so it will dephosphrylate)
dataS.kS(2049)=0; %ppAKT binding BAD
end

%BIM
if m==2
[dataS,dataG]=RunPrep;
dataS.kS(2045)=0; %ppERK can't bind BIM (so it will dephosphrylate)
dataS.kS(1035)=0; %ppAKT binding FOXO
end


STIM(156:162)=[3.3,0,0,0,0,0,1721]; %ligs
th=48;
filename=strcat('Apoptosis2S_Mech2_',num2str(m));
Run_SimCells(STIM,th,filename,dataS)

end

%% Figure 5D-F
% Dual Inhibitor Investigations

% Get random population of cells (if necessary)
%RunRandomPopCells(400,[matpath,'RandomPopCells400.mat']);

% Run sims
m=5;
STIM=zeros(NumSpecies,1);
STIM(156:162)=[3.3,0,0,0,0,0,1721]; %ligs
STIM(769)=10000; %MEKi
STIM(771)=10000; %AKTi
th=80;


flagD=0;
St=load([matpath,'RandomPopCells400.mat']);
cells0=St.cells0;

for i=1:length(cells0)
    
    xoutG=cells0{i}.xoutG_all;
    xoutS=cells0{i}.xoutS_all;
    
    [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,xoutS,xoutG',[],[]);
    
    inds=1:30:length(tout_all);

    C.tout_all=tout_all(inds);
    C.xoutG_all=xoutG_all(inds,:);
    C.xoutS_all=xoutS_all(inds,:);

    cells{i,1}=C;
    
    disp(strcat('cell number =',num2str(i)))

end

filename=strcat('Apoptosis2S_400cells_',num2str(m));
txt=strcat(matpath,filename,'.mat');
save(txt,'-v7.3','cells');

%% Figure S6A

% Get random population of cells
%RunRandomPopCells(25,[matpath,'RandomPopCells.mat']);

% Run sims
STIM=zeros(NumSpecies,1);
STIM(156:162)=[3.3,0,0,0,0,0,1721]; %ligs
STIM(84)=0.038; %TRAIL
th=80;
filename='Apoptosis2S_TRAIL_1ng';
Run_SimCells(STIM,th,filename)

STIM=zeros(NumSpecies,1);
STIM(156:162)=[3.3,0,0,0,0,0,1721]; %ligs
STIM(84)=0.0038; %TRAIL
th=80;
filename='Apoptosis2S_TRAIL_01ng';
Run_SimCells(STIM,th,filename)

%% Figure S6B

% Details available upon request.

%% Figure 6A

% No sims

%% Figure 6B

% Get random population of cells if necessary
%RunRandomPopCells(50,[matpath,'RandomPopCells.mat']);

th=30;

% EGF only
STIM=zeros(NumSpecies,1);
STIM(156:162)=[10,0,0,0,0,0,0];
Run_SimCells(STIM,th,'CellCycle2S_E')

% INS only
STIM=zeros(NumSpecies,1);
STIM(156:162)=[0,0,0,0,0,0,1721];
Run_SimCells(STIM,th,'CellCycle2S_I')

% EGF+INS
STIM=zeros(NumSpecies,1);
STIM(156:162)=[10,0,0,0,0,0,1721];
Run_SimCells(STIM,th,'CellCycle2S_E+I')

% EGF+INS+MEKi
STIM=zeros(NumSpecies,1);
STIM(156:162)=[10,0,0,0,0,0,1721];
STIM(769)=10000;
Run_SimCells(STIM,th,'CellCycle2S_E+I+MEKi')

% EGF+INS+AKTi
STIM=zeros(NumSpecies,1);
STIM(156:162)=[10,0,0,0,0,0,1721];
STIM(771)=10000;
Run_SimCells(STIM,th,'CellCycle2S_E+I+AKTi')

%% Figure 6E

% Get simulations for different ppERK and ppAKT levels
flagD=1;
th=6;
STIM=zeros(775,1);
as=1:10;
es=1:0.5:5.5;
xoutS_alls=[];
[dataS,dataG]=RunPrep;
dataS_up=dataS;
for i=1:length(es)
    dataS_up.kS([910,913])=dataS.kS([910,913])*es(i);%ppERK

    for k=1:length(m)
        dataS_up.kS([985,991,995])=dataS.kS([985,991,995])*as(k);%ppAKT       
        
        [tout_all,~,xoutS_all]=RunModel(flagD,th,STIM,[],[],dataS_up,[]);
        
        disp([mean(xoutS_all(:,718)),mean(xoutS_all(:,697)),mean(sum(xoutS_all(:,[44,45,46,47,80]),2))])
        
        inds=1:30:length(tout_all);
        xoutS_alls{i,k}=xoutS_all(inds,:);
        disp(i)
    end
end

% Get individual cases
NumSpecies=775;
xoutS_alls_cases=[];
flagD=1;
th=6;
% EGF only
STIM=zeros(NumSpecies,1);
STIM(156:162)=[10,0,0,0,0,0,0];
[tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],[],[]);
xoutS_alls_cases{1}=xoutS_all;
% INS only
STIM=zeros(NumSpecies,1);
STIM(156:162)=[0,0,0,0,0,0,1721];
[tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],[],[]);
xoutS_alls_cases{2}=xoutS_all;
% EGF+INS
STIM=zeros(NumSpecies,1);
STIM(156:162)=[10,0,0,0,0,0,1721];
[tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],[],[]);
xoutS_alls_cases{3}=xoutS_all;


txt=strcat(matpath,'CellCycle_erkaktcycd','.mat');
save(txt,'-v7.3','xoutS_alls_cases','xoutS_alls');

%% Figure 6F-G

% EGF+INS 400 cells

% Run if necessary:
%RunRandomPopCells(400,[matpath,'RandomPopCells400.mat']);


STIM=zeros(NumSpecies,1);
STIM(156:162)=[10,0,0,0,0,0,1721];
th=30;

flagD=0;
St=load([matpath,'RandomPopCells400.mat']);
cells0=St.cells0;

for i=1:length(cells0)
    
    xoutG=cells0{i}.xoutG_all;
    xoutS=cells0{i}.xoutS_all;
    
    [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,xoutS,xoutG',[],[]);
    
    inds=1:30:length(tout_all);

    C.tout_all=tout_all(inds);
    C.xoutG_all=xoutG_all(inds,:);
    C.xoutS_all=xoutS_all(inds,:);

    cells{i,1}=C;
    
    disp(strcat('cell number =',num2str(i)))

end

txt=strcat(matpath,'CellCycle2S_E+I_400cells','.mat');
save(txt,'-v7.3','cells');

%% Figure S7
% AKT levels investigation

[dataS,~]=RunPrep;

% E
STIM=zeros(775,1);
STIM(156)=10;
STIM(162)=0;
th=24;
filename='CellCycle_E_25';
Run_SimCells(STIM,th,filename,dataS)

% E+ AKTi 1
STIM=zeros(775,1);
STIM(156)=10;
STIM(162)=1721;
STIM(771)=1;
th=24;
filename='CellCycle_EAKTi1_25';
Run_SimCells(STIM,th,filename,dataS)


% E+ AKTi 2
STIM=zeros(775,1);
STIM(156)=10;
STIM(162)=1721;
STIM(771)=2;
th=24;
filename='CellCycle_EAKTi2_25';
Run_SimCells(STIM,th,filename,dataS)

% E+I
STIM=zeros(775,1);
STIM(156)=10;
STIM(162)=1721;
th=24;
filename='CellCycle_EI_25';
Run_SimCells(STIM,th,filename,dataS)

%% Figure 7A
% Amplification/Deletion Test

%Define indices
inds_genes_oncogenes={58 
59	
60	
61	
71  
72	
73  
74  
75  
76	
77
78
79
80	
82
83	
84
85
86	
88	
89	
90
91	
93
94
95
96	
98	
99	
100	
101	
103
104
105	
106	
107	
116
117	
119	
120	
121
122	
124	
131	
132
133	
140};	

inds_genes_tumorsupressors={68
69	
87	
92	
97	
102	
108	
109	
110   
123};

% Run Sims
numcells=50;

for l=1:3
    if l==1; inds2use=[]; end
    if l==2; inds2use=inds_genes_oncogenes; end
    if l==3; inds2use=inds_genes_tumorsupressors; end
    
    if isempty(inds2use)
        iterations=1;
    else iterations=length(inds2use);
    end
    
for k=1:iterations

[~,dataG]=RunPrep;

if isempty(inds2use)
    ind_current=[];
else ind_current=inds2use{k};
end

xoutG=dataG.x0gm_mpc_D;
pathi='initialized/';
kTCleak=dlmread(strcat(pathi,'i_kTCleakF.txt'));
kTCmaxs=dlmread(strcat(pathi,'i_kTCmaxsF.txt'));

if l==2 %oncogenes
    xoutG(ind_current+141*2)=xoutG(ind_current+141*2)*10;
    kTCleak(ind_current)=kTCleak(ind_current)*10;
    kTCmaxs(ind_current)=kTCmaxs(ind_current)*10;
end

if l==3 %tumor supressors
    xoutG(ind_current+141*2)=xoutG(ind_current+141*2)/10;
    kTCleak(ind_current)=kTCleak(ind_current)/10;
    kTCmaxs(ind_current)=kTCmaxs(ind_current)/10;
end

% Run to steady-state, deterministic
flagD=1;
th=100;
STIM=zeros(775,1);
[tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],xoutG,[],dataG,kTCleak,kTCmaxs);

% Run stochastic simulations
flagD=0;
indsD=[dataG.indsD];
xoutG=dataG.x0gm_mpc;
xoutG([indsD,indsD+141])=dataG.x0gm_mpc_D([indsD,indsD+141]); %Set genes that don't vary to their deterministic amounts
xoutG(141*2+1:end)=xoutG_all(end,141*2+1:end); %Set mRNA to proper amount

cells=[];
for i=1:numcells
    th=24;
    STIM=zeros(775,1);
    [tout_all1,xoutG_all1,xoutS_all1]=RunModel(flagD,th,STIM,xoutS_all(end,:),xoutG,[],dataG,kTCleak,kTCmaxs);
    
    th=24;
    STIM(156:162)=[10,0,0,0,0,0,1721];
    [tout_all2,xoutG_all2,xoutS_all2]=RunModel(flagD,th,STIM,xoutS_all1(end,:),xoutG_all1(end,:)',[],dataG,kTCleak,kTCmaxs);
    
    inds=1:30:length(tout_all2);%th*3600/30; %inds to thin data;
    C.tout_all=tout_all2(inds);
    C.xoutG_all=xoutG_all2(inds,:);
    C.xoutS_all=xoutS_all2(inds,:);
    cells{i}=C;
end

if l==1
    txt=strcat(matpath,'mut_none.mat');
else txt=strcat(matpath,'mut_',num2str(inds2use{k}),'.mat');
end
save(txt,'-v7.3','cells');

end


end

%% Figure 7B-D
% MEK Mechanism Investigation

[dataS,dataG]=RunPrep;

% Increase affinity of MEK for Raf
dataS.kS([890,896,902])=dataS.kS([890,896,902])*10;

% Run to steady-state, deterministic
flagD=1;
th=100;
STIM=zeros(775,1);
[tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],dataS,[]);

numcells=50;
cells=[];
flagD=0;
for i=1:numcells
    th=24;
    STIM=zeros(775,1);
    [tout_all1,xoutG_all1,xoutS_all1]=RunModel(flagD,th,STIM,xoutS_all(end,:),xoutG_all(end,:)',dataS,[]);
    
    th=24;
    STIM(156:162)=[10,0,0,0,0,0,1721];
    [tout_all2,xoutG_all2,xoutS_all2]=RunModel(flagD,th,STIM,xoutS_all1(end,:),xoutG_all1(end,:)',dataS,[]);
    
    inds=1:30:length(tout_all2);%th*3600/30; %inds to thin data;
    C.tout_all=tout_all2(inds);
    C.xoutG_all=xoutG_all2(inds,:);
    C.xoutS_all=xoutS_all2(inds,:);
    cells{i}=C;
end

txt=strcat(matpath,'MEKmechanism','.mat');
save(txt,'-v7.3','cells');

%% Figure 8 Part 1
% MCF10A

% Create Random population of cells (if needed)
RunRandomPopCells(50,[matpath,'RandomPopCells.mat'])
NumSpecies=775;

% Effect of inhibitors and GFs, TRAIL dose response in SS state.
m=1; %NO STIM
STIM=zeros(NumSpecies,1);
STIM(156:162)=[0,0,0,0,0,0,0]; %ligs
STIM(769)=0; %MEKi
STIM(771)=0; %AKTi
th=80;
filename=strcat('Apoptosis2S_10A_',num2str(m));
Run_SimCells(STIM,th,filename)

m=2; %MEKi
STIM=zeros(NumSpecies,1);
STIM(156:162)=[0,0,0,0,0,0,0]; %ligs
STIM(769)=10000; %MEKi
STIM(771)=0; %AKTi
th=80;
filename=strcat('Apoptosis2S_10A_',num2str(m));
Run_SimCells(STIM,th,filename)

m=3; %AKTi
STIM=zeros(NumSpecies,1);
STIM(156:162)=[0,0,0,0,0,0,0]; %ligs
STIM(769)=0; %MEKi
STIM(771)=10000; %AKTi
th=80;
filename=strcat('Apoptosis2S_10A_',num2str(m));
Run_SimCells(STIM,th,filename)

m=4; %MEKi+AKTi
STIM=zeros(NumSpecies,1);
STIM(156:162)=[0,0,0,0,0,0,0]; %ligs
STIM(769)=10000; %MEKi
STIM(771)=10000; %AKTi
th=80;
filename=strcat('Apoptosis2S_10A_',num2str(m));
Run_SimCells(STIM,th,filename)

%% Figure 8 Part 2
% u87

% To run the u87 model, you must uncomment "datasheet_in='master_u87.xlsx';
% ind_u87=1;" in the first section of the RunPrep.m file. Also, beware, running
% RunInitialize (below) will overwrite all files in initialized/ (all files associated with MCF10A cells). 

% % u87 sims
% [dataS,dataG]=RunPrep;
% dataG.tck50as(98,6)=dataG.tck50as(98,6)*2; %CHANGED FROM MCF10A (change doesn't change results, simply allows initialization to proceed)
% dataG.tck50as(98,5)=dataG.tck50as(98,5)*2; %CHANGED FROM MCF10A (change doesn't change results, simply allows initialization to proceed)
% 
% mods=[1,1,1,1,1,1,1,1];
% RunInitialize(mods,dataS,dataG); % Be aware: This will overwrite files in /initialized


% % Create Random population of cells (if needed)
% RunRandomPopCells(25,[matpath,'RandomPopCells.mat'])
% NumSpecies=775;

% Effect of inhibitors and GFs, TRAIL dose response in SS state.
m=1; %NO STIM
STIM=zeros(NumSpecies,1);
STIM(156:162)=[0,0,0,0,0,0,0]; %ligs
STIM(769)=0; %MEKi
STIM(771)=0; %AKTi
th=48;
filename=strcat('Apoptosis2S_u87_',num2str(m));
Run_SimCells(STIM,th,filename)

% m=2; %MEKi
% STIM=zeros(NumSpecies,1);
% STIM(156:162)=[0,0,0,0,0,0,0]; %ligs
% STIM(769)=10000; %MEKi
% STIM(771)=0; %AKTi
% th=48;
% filename=strcat('Apoptosis2S_u87_',num2str(m));
% Run_SimCells(STIM,th,filename)
% 
% m=3; %AKTi
% STIM=zeros(NumSpecies,1);
% STIM(156:162)=[0,0,0,0,0,0,0]; %ligs
% STIM(769)=0; %MEKi
% STIM(771)=10000; %AKTi
% th=48;
% filename=strcat('Apoptosis2S_u87_',num2str(m));
% Run_SimCells(STIM,th,filename)
% 
% m=4; %MEKi+AKTi
% STIM=zeros(NumSpecies,1);
% STIM(156:162)=[0,0,0,0,0,0,0]; %ligs
% STIM(769)=10000; %MEKi
% STIM(771)=10000; %AKTi
% th=48;
% filename=strcat('Apoptosis2S_u87_',num2str(m));
% Run_SimCells(STIM,th,filename)

