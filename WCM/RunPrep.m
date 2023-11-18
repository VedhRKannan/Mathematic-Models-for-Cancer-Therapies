function [dataS,dataG]=RunPrep

% This script loads and calculates many of the intial parameters needed to
% run the model.

%% GET DATA 
datasheet_gp='params.xlsx'; %general parameters
datasheet_in='gene_MCF10A.xlsx'; %genetic inputs (initial gene, mRNA, and protein levels)
%datasheet_in='gene_u87.xlsx'; ind_u87=1; %uncomment to use u87 data

%% Model Essentials
ts=30; %Timestep for simulations;
% General params
datasheet=datasheet_gp;
[cellparams,~,~]=xlsread(datasheet,'cellparams');
VolumeofCell=cellparams(1);
[Vn,Vc,Vm,mpc2nmcf_Vc,mpc2nmcf_Vm,mpc2nmcf_Vn]=CalcVolumeParams(VolumeofCell);
Ve=cellparams(3);
mT0=cellparams(9)*mpc2nmcf_Vc;
PIP2_0=cellparams(10);
M0=cellparams(11);
[tlparams,~,~]=xlsread(datasheet,'tlparams');
Rt=tlparams(1);
kbRi=tlparams(2);
kdR0=tlparams(3);
kbR0=0;
nR=tlparams(4);
k50R=tlparams(5);
kT1=tlparams(6);
kT2=tlparams(7);
kT3=tlparams(8);
kT4=tlparams(9);
k50E=tlparams(10);
[gene,~,~]=xlsread(datasheet,'gene');
kTCd=gene(:,12);
kTLd=gene(:,13);
kTLnat=gene(:,16);
kGin=gene(:,17);
kGac=gene(:,18);
kTLnat(10:12)=kTLnat(10:12)*5;

%Inputs
datasheet=datasheet_in;
[gene,~,~]=xlsread(datasheet);
gExp_mpc=gene(:,1); %gene copy numbers experimental (molecules/cell)
mExp_mpc=gene(:,2); %mRNA expression experimental (molecules/cell)
mExp_nM=mExp_mpc*mpc2nmcf_Vc;


% U87 modifications
if exist('ind_u87','var')
kTLnat(102)=0; %PTEN translation equals zero (therefore PTEN levels equal zero)
end

%Important indices
numberofgenes=141;
numberofTARs=7;
indsCC=40:78;
indsDD=2:27;
indsM=122;
indsPIP2=691;
indsmT=773;


%% PARCDL prep
% PARCDL STOICHIOMETRIC MATRIX
S_PARCDL=csvread('PARCDL_sm18.csv',1,1);
S_PARCDL=sparse(S_PARCDL);
S_TL=S_PARCDL(:,3:3+numberofgenes-1);
S_d=S_PARCDL(:,1844:end);

% RATE CONSTANTS SYNTHESIS AND DEGRADATION
EIF4Efree=(kTLnat(131).*mExp_nM(131))./kTLd(131);
[kTLCd,kTL,kXd,xp_mpc] = PARCDL_rateconstants_sd(k50E,mExp_mpc,kTLnat,kTLd,EIF4Efree,S_TL,S_d);

% RATE CONSTANTS MAIN
[kD,kC,kA,kR,kP,...
    kRP1,kRP2,kRP3,kRP4,kRP5,kRP6,kRP7,kRP8,kRP9,kRP10,kRP11,kRP12,kRP13,kRP14,kRP15,kRP16,kRP17,kRP18,kRP19,kRP20,kRP21,kRP22,kRP23,kRP24,kRP25,kRP26,kRP27,kRP28,kRP29,kRP30,kRP31,kRP32,kRP33,kRP34,...
    kDP,kPA] = PARCDL_rateconstants_main(Vc,Ve,Vm,Vn);
kE(1,1)=kTLCd(102)*mT0; 
kE(2,1)=kTLCd(102);
kE(3,1)=kTLCd(102);

% PARCDL COMPARTMENTS
flagOUT=[1,1,1];
[VxPARCDL,VvPARCDL,VxTL]=PARCDL_compartments(flagOUT,Vc,Ve,Vm,Vn,S_d);

% SET MANUAL INITIAL CONDITIONS
CellCycleSpecies=[80000,0.0023875,3.2308e-05,11012,0.0013746,0.0036083,0.018044,0.0037528,2.5164,8.7989,27.119,114.09,11.28,1412.9,489.7,160.2,552.84,39.644,138.62,52.721,13.158,207.98,6.0486,1087.9,116.34,42.027,420.6,34.408,38.992,6.8625,711.67,7.3241,94.317,265.18,167.41,0.85635,2.0389e-117,88094,0.0013145];
CellCycleSpecies=CellCycleSpecies/10;
DNADamageSpecies=[296.62,6.2458,205.62,2.2305,0,0,6.2458,6.2458,6.2458,6.2457,6.2457,6.2457,6.2457,6.2457,6.2457,6.2457,6.2458,6.2458,6.2457,6.2457,6.2457,6.2457,6.2457,6.2457,6.2457,6.2457];

% SPECIES
pExp_nM=xp_mpc.*(1E9./(VxTL*6.023E23));
x0PARCDL=S_TL*pExp_nM;
x0PARCDL(1)=Rt;
x0PARCDL(indsCC)=CellCycleSpecies;
x0PARCDL(indsDD)=DNADamageSpecies;
x0PARCDL(indsM)=M0;
x0PARCDL(indsPIP2)=PIP2_0;
x0PARCDL(indsmT)=mT0;


%% PARCDL Rate constant MODIFICATIONS

% A
kA(17)=kA(17)/100*10; %pC8 binding C6
kA(14)=kA(14)/10000*10; %C3 binding pC6
kA(13)=kA(13)/10000*10; %C8 cleaving pC3
kA(26)=kA(26)/1000*10; %C8 binding to Bid
kA(23)=kA(23)/1000*10; %C3 binding to PARP
kA([38,42])=kA([38,42])*100000; %Baxm/Bax2 dimerization
kA([39,43])=kA([39,43])*10000; %Bax2/Bax4 dissociation
kA(64)=kA(64)/10000; %Apoptosome cleaving pC3
kA(28)=kA(28)/10000; %C8 cleaving Bid
kA([50,53])=kA([50,53])*1000; %koff CytochromeC Open Mitochondiral Pores. 
kA([51,54])=kA([51,54])/100; %Catalytic constant for above reaction.
kA(3)=kA(3)/10; %Ligand-bound receptor becomming activated.
kA(85)=0.001; %Translocation to mito, BCL2c
kA(86)=0.1; %Translocation to cyto
kA(87)=0; %Basal pC8 cleavage. Set during initialization.
kA(77)=0; %BIM binding Bax

% P
kP(103)=0;
kP(60)=0;
kP([117,121,178,179])=kP([117,121,178,179])*0.1;
kP(154)=kP(154)/10; 
kP(158)=0.001*5; 

% PA
kPA(4)=kPA(4)*5;

% C
kC([59,46,48,50])=kC([59,46,48,50])/5;

% kXd
kXd(2)=kXd(2)*1000; 
kXd([5,8,10,16,19])=kXd([5,8,10,16,19])*100;
kXd(44)=kXd(44)*10;
kXd(46)=kXd(46)*10;
kXd(562)=kXd(562)*10;


%% gm
indsD=[6:9,13:30]; %indices of genes that will not be stochastically regulated.
[xgac_mpc,xgin_mpc,xgac_mpc_D,xgin_mpc_D,kTCleak,AllGenesVec,GenePositionMatrix] = gm_Prep(mExp_mpc,gExp_mpc,kTCd,kGac,kGin,numberofgenes);

% Set kTCmaxs
kTCmax=0.1;
kTCmaxs=ones(numberofgenes,1)*kTCmax;

% Transcriptional Activators
tcnas=ones(numberofgenes,numberofTARs); %Number of genes by number of TARs.
tcnas(10:12,1)=3; %pcFos_cJun act
tcnas(99,1)=3; %pcFos_cJun act on cJun
tcnas(10:12,2)=3; %cMyc act
tcnas([26,53,54],3)=4; %p53ac
tcnas([55,58,59,60,61,63,65,66,127,128,136,140],4)=4; %FOXOnuc
tcnas([68,92,97,98],5)=4; %ppERKnuc
tcnas([68,92,97,98],6)=4; %pRSKnuc
tcnas(100,7)=4; %bCATENIN
 
% Transcriptional Repressors
tcnrs=ones(numberofgenes,numberofTARs);
tcnrs(98,1)=4;

% k50 for TAR Activation
tck50as=zeros(numberofgenes,numberofTARs);
%pcFos_cJun
tck50as(10:12,1)=1.25; %CyclinD
tck50as(99,1)=0.8; %cJun
%cMyc
tck50as(10:12,2)=450; %CyclinD
%p53ac
tck50as(26,3)=50;%p21
tck50as(53:54,3)=1350;%PUMA,NOXA
%FOXO
tck50as(55,4)=45;%19
tck50as([58,59,60,61,63,65,66,127,128,136,140],4)=60; %RTKs
%ppERKnuc
tck50as(68,5)=65; %SPRY2
tck50as([92,97],5)=40; %DUSPs
tck50as(98,5)=20; %cFos
%pRSKnuc
tck50as(68,6)=20; %SPRY2
tck50as([92,97],6)=10; %DUSPs
tck50as(98,6)=5; %cFos
%bCATENINnuc
tck50as(100,7)=250;%cMyc

% k50 for TAR Repression
tck50rs=zeros(numberofgenes,numberofTARs);
tck50rs(98,1)=tck50as(99,1); %cFos

% Convert to molecules per cell
tck50as=tck50as*(1/mpc2nmcf_Vn);
tck50rs=tck50rs*(1/mpc2nmcf_Vn);

% SAVE PARAMETERS IN STRUCTURE
kS=[Rt;...
EIF4Efree;...
Vc;...
Ve;...
Vm;...
Vn;...
kT1;...
kT2;...
kT3;...
kT4;...
k50E;...
kbR0;...
kbRi;...
kdR0;...
nR;...
k50R;...
kTL;...
kTLd;...
kTLCd;...
kD;...
kC;...
kA;...
kR;...
kP;...
kRP1;...
kRP2;...
kRP3;...
kRP4;...
kRP5;...
kRP6;...
kRP7;...
kRP8;...
kRP9;...
kRP10;...
kRP11;...
kRP12;...
kRP13;...
kRP14;...
kRP15;...
kRP16;...
kRP17;...
kRP18;...
kRP19;...
kRP20;...
kRP21;...
kRP22;...
kRP23;...
kRP24;...
kRP25;...
kRP26;...
kRP27;...
kRP28;...
kRP29;...
kRP30;...
kRP31;...
kRP32;...
kRP33;...
kRP34;...
kDP;...
kPA;...
kXd;...
kE];

flagE=1;
dataS=struct;
dataS.ts=ts; %time step
dataS.x0PARCDL=x0PARCDL;
dataS.kS=kS;
dataS.VvPARCDL=VvPARCDL;
dataS.VxPARCDL=VxPARCDL;
dataS.S_PARCDL=S_PARCDL;
dataS.mExp_nM=mExp_nM;
dataS.mMod=mExp_nM;
dataS.pExp_nM=pExp_nM;
dataS.flagE=flagE;

% CALCULATE CELL CYCLE mRNAs
indsccx=[6:9,13:25];
mExp_mpc(indsccx)=17; 
mExp_nM(indsccx)=mExp_mpc(indsccx)*mpc2nmcf_Vc;

%UPDATE CERTAIN DATA SPECIES
dataS.mExp_nM=mExp_nM;
dataS.mMod=mExp_nM;


%% dataG
x0gm_mpc=[xgac_mpc;xgin_mpc;mExp_mpc]; 
x0gm_mpc_D=[xgac_mpc_D;xgin_mpc_D;mExp_mpc];

dataG=struct;
dataG.x0gm_mpc=x0gm_mpc;
dataG.x0gm_mpc_D=x0gm_mpc_D;
dataG.kGin=kGin;
dataG.kGac=kGac;
dataG.kTCleak=kTCleak;
dataG.kTCmaxs=kTCmaxs;
dataG.kTCd=kTCd;
dataG.tcnas=tcnas;
dataG.tcnrs=tcnrs;
dataG.tck50as=tck50as;
dataG.tck50rs=tck50rs;
dataG.GenePositionMatrix=GenePositionMatrix;
dataG.AllGenesVec=AllGenesVec;
dataG.Vn=Vn;
dataG.indsD=indsD;


