function RunInitialize(mods,dataS,dataG)


% This function runs the initialization procedure. It creates text files
% that are imported into RunModel.m. It incorporates many basal activities
% to define a steady-state for a cell in the serum-starved state. See
% Supplemental Information for more details.

% USAGE:
% To run full initialization procedure, select all modules:
% RunInitialize([1,1,1,1,1,1,1,1],[],[])



% PREP
if isempty(dataS)
    [dataS,~]=RunPrep;
end
if isempty(dataG)
    [~,dataG]=RunPrep;
end

dataS_up=dataS;

pathi='initialized/';
mkdir(pathi)

kTL=dataS.kS(17:157);
EIF4Efree=dataS.kS(2);
Vn=dataS.kS(6);
k50E=dataS.kS(11);
numberofgenes=size(dataG.tcnas,1);
mpc2nmcf_Vn=1E9/(Vn*6.023E+23);
indsEIF4E=763;
S_TL=dataS.S_PARCDL(:,3:3+numberofgenes-1);

tck50as=dataG.tck50as;
tck50rs=dataG.tck50rs;
tcnas=dataG.tcnas;
tcnrs=dataG.tcnrs;
kTCd=dataG.kTCd;
x0gm_mpc=dataG.x0gm_mpc;
x0gm_mpc_D=dataG.x0gm_mpc_D;
mExp_mpc=x0gm_mpc(283:423);
kTCmaxs=dataG.kTCmaxs;

options=odeset('RelTol',1e-6,'NonNegative',1:length(dataS.x0PARCDL));



%% Step #1: Basal Ras and PIP Activities; No EIF4E

if mods(1)
disp('Step #1: Basal Ras and PIP Activities; No EIF4E')

disp('kTLAdjust...')
flagR=0; 
flagE=0; 
xoutIN=dataS.x0PARCDL';
kTLIN=kTL;
[toutNEW1,xoutNEW1,kTLNEW,dataS_up] = kTLAdjustWhile(pathi,flagR,flagE,xoutIN,kTLIN,S_TL,dataS,dataS_up);
dataS_up.kS(17:157)=kTLNEW;

EIF4Efreenew=xoutNEW1(end,indsEIF4E);
kTLNEW=kTLNEW*((EIF4Efree/(k50E+EIF4Efree))/(EIF4Efreenew/(k50E+EIF4Efreenew)));
dataS_up.kS(17:157)=kTLNEW;
dataS_up.kS(2)=EIF4Efreenew;

disp('Running to steady-state with KTLNEWs...')
tstart=0;
tend=3600*1000;
tspan=[tstart,tend];
xout_update=xoutNEW1(end,:);
[toutNEW1a,xoutNEW1a]=ode15s(@createODEs,tspan,xout_update,options,dataS_up);

dlmwrite([pathi,'i_kTLNEW_1.txt'],kTLNEW,'Precision',10)
dlmwrite([pathi,'i_xoutNEW1.txt'],xoutNEW1,'Precision',10)
dlmwrite([pathi,'i_toutNEW1.txt'],toutNEW1)

end

%% Step #1: Basal Ras and PIP Activities; Add back EIF4E

if mods(2)
disp('Step #1: Basal Ras and PIP Activities; Add back EIF4E')

xoutNEW1=dlmread([pathi,'i_xoutNEW1.txt']);
kTLIN=dlmread([pathi,'i_kTLNEW_1.txt']);
dataS_up.kS(17:157)=kTLIN;
dataS_up.kS(2)=xoutNEW1(end,indsEIF4E);

xoutIN=xoutNEW1(end,:);
flagR=0;
flagE=1;
[toutNEW2,xoutNEW2,kTLNEW,dataS_up] = kTLAdjustWhile(pathi,flagR,flagE,xoutIN,kTLIN,S_TL,dataS,dataS_up);


dlmwrite([pathi,'i_kTLNEW_2.txt'],kTLNEW,'Precision',10)
dlmwrite([pathi,'i_xoutNEW2.txt'],xoutNEW2,'Precision',10)
dlmwrite([pathi,'i_toutNEW2.txt'],toutNEW2);
end

%% Step #1: Basal Ras and PIP Activities; Add back Ribosome Parameters

if mods(3)
disp('Step #1: Basal Ras and PIP Activities; Add back Ribosome Parameters')

xoutNEW2=dlmread([pathi,'i_xoutNEW2.txt']);
kTLIN=dlmread([pathi,'i_kTLNEW_2.txt']);
dataS_up.kS(17:157)=kTLIN;
dataS_up.kS(2)=xoutNEW2(end,indsEIF4E);

xoutIN=xoutNEW2(end,:);
flagR=1;
flagE=1;
[toutNEW3,xoutNEW3,kTLNEW,dataS_up] = kTLAdjustWhile(pathi,flagR,flagE,xoutIN,kTLIN,S_TL,dataS,dataS_up);

dlmwrite([pathi,'i_kTLNEW_3.txt'],kTLNEW,'Precision',10)
dlmwrite([pathi,'i_xoutNEW3.txt'],xoutNEW3,'Precision',10)
dlmwrite([pathi,'i_toutNEW3.txt'],toutNEW3);
dlmwrite([pathi,'i_Rt.txt'],dataS_up.kS(1),'Precision',10);
dlmwrite([pathi,'i_EIF4Efree.txt'],dataS_up.kS(2),'Precision',10);

end

%% Step #2: Basal Cyclin D Synthesis and p21 Degradation

if mods(4)
disp('Step #2: Basal Cyclin D Synthesis and p21 Degradation')

flagE=1;
xoutNEW3=dlmread([pathi,'i_xoutNEW3.txt']);
kTLIN=dlmread([pathi,'i_kTLNEW_3.txt']);
dataS_up.kS(17:157)=kTLIN;
kbR0=dlmread([pathi,'i_kbR0.txt']);
dataS_up.kS(12)=kbR0;
dataS_up.flagE=flagE;
Rt=dlmread([pathi,'i_Rt.txt']);
EIF4Efree=dlmread([pathi,'i_EIF4Efree.txt']);
dataS_up.kS(1)=Rt;
dataS_up.kS(2)=EIF4Efree;
totalcyclinDfromdata=sum(dataS.pExp_nM(10:12));
totalp21fromdata=dataS.pExp_nM(26);

kC173=dataS.kS(632);
kC82=dataS.kS(541);
xoutIN=xoutNEW3(end,:);
xout_update=xoutIN;
tstart=0;
tend=3600*1000;
tspan=[tstart,tend];

th=0.001;
n=1;
ratio_cd=0.5;
ratio_p21=0.5;

while or(and(ratio_cd<=1+th,ratio_cd>=1-th)<1,and(ratio_p21<=1+th,ratio_p21>=1-th)<1) %Until every ratio is between two limits.

[toutNEW4,xoutNEW4]=ode15s(@createODEs,tspan,xout_update,options,dataS_up);
xout_update=xoutNEW4(end,:);

% Check if cell cycled
if logical(sum(xoutNEW4(:,69)>35));
    disp('Cell cycled with basal levels of cyclinD')
    return
end

ratio_cd=totalcyclinDfromdata/sum(xoutNEW4(end,[44,45,46,47,80]));
if ~and(ratio_cd<=1+th,ratio_cd>=1-th)
    F=1+((ratio_cd-1)/2);
    kC173=kC173*F;
    dataS_up.kS(632)=kC173;
else disp('Threshold Reached!');
end

ratio_p21=totalp21fromdata/sum(xoutNEW4(end,79:83));
if ~and(ratio_p21<=1+th,ratio_p21>=1-th)
    F=1+((ratio_p21-1)/2);
    kC82=kC82/F;
    dataS_up.kS(541)=kC82;
else disp('Threshold Reached!');
end


disp([n, ratio_cd, ratio_p21])
n=n+1;
end


dlmwrite([pathi,'i_kC173.txt'],kC173);
dlmwrite([pathi,'i_kC82.txt'],kC82);
dlmwrite([pathi,'i_xoutNEW4.txt'],xoutNEW4,'Precision',10)
dlmwrite([pathi,'i_toutNEW4.txt'],toutNEW4);

end

%% Step #3: Basal Caspase 8 Cleavage

if mods(5)
disp('Step #3: Basal Caspase 8 Cleavage')

% Import necessary parameters and species
xoutNEW4=dlmread([pathi,'i_xoutNEW4.txt']);
kTLIN=dlmread([pathi,'i_kTLNEW_3.txt']);
kbR0=dlmread([pathi,'i_kbR0.txt']);
kC173=dlmread([pathi,'i_kC173.txt']);
kC82=dlmread([pathi,'i_kC82.txt']);
Rt=dlmread([pathi,'i_Rt.txt']);
EIF4Efree=dlmread([pathi,'i_EIF4Efree.txt']);
dlmwrite([pathi,'i_kA77.txt'],3.162075e-09,'Precision',7);
kA77=dlmread([pathi,'i_kA77.txt']);
% Assign to data structure
dataS_up.kS(17:157)=kTLIN;
dataS_up.kS(12)=kbR0;
dataS_up.kS(632)=kC173;
dataS_up.kS(541)=kC82;
dataS_up.kS(1)=Rt;
dataS_up.kS(2)=EIF4Efree;
dataS_up.kS(709)=kA77;

% Set up
xout_update=xoutNEW4(end,:);

disp('START LOOP')
kA87s=[1E-10 1E-9 1E-8 1E-7 1E-6 1E-5 1E-4 1E-3 1E-2 1E-1];

for i=1:length(kA87s)
    kA87=kA87s(i);
    dataS_up.kS(719)=kA87;
    
    % fix kTL constants
    xoutIN=xout_update;
    flagR=1; 
    flagE=1; 
    [toutNEW5,xoutNEW5,kTLNEW,dataS_up,flagA] = kTLAdjustWhile(pathi,flagR,flagE,xoutIN,kTLIN,S_TL,dataS,dataS_up);

    if ~flagA %If apoptosis didn't happen, save all important parameters to save if cell dies during next iteration.
        xoutNEW5LAST=xoutNEW5;
        toutNEW5LAST=toutNEW5;
        kTLNEWLAST=kTLNEW;
        disp ('Cell did not die. Trying higher number.')        
    end

if flagA;
    dlmwrite([pathi,'i_kA87.txt'],kA87s(i-1),'Precision',10); disp('DONE');
    dlmwrite([pathi,'i_kTLNEW_5.txt'],kTLNEWLAST,'Precision',10)
    dlmwrite([pathi,'i_xoutNEW5.txt'],xoutNEW5LAST,'Precision',10)
    dlmwrite([pathi,'i_toutNEW5.txt'],toutNEW5LAST);
    break
end
    disp(i)
end


end

%% Step #4: Basal DNA Damage

if mods(6)
disp('Step #4: Basal DNA Damage')
    
xoutNEW5=dlmread([pathi,'i_xoutNEW5.txt']); xoutIN=xoutNEW5(end,:);
kTLNEW_5=dlmread([pathi,'i_kTLNEW_5.txt']); kTLIN=kTLNEW_5;
kbR0=dlmread([pathi,'i_kbR0.txt']);
kC173=dlmread([pathi,'i_kC173.txt']);
kC82=dlmread([pathi,'i_kC82.txt']);
kA77=dlmread([pathi,'i_kA77.txt']);
kA87=dlmread([pathi,'i_kA87.txt']);
Rt=dlmread([pathi,'i_Rt.txt']);
EIF4Efree=dlmread([pathi,'i_EIF4Efree.txt']);
% modifying data.S structure
dataS_up.kS(1)=Rt;
dataS_up.kS(2)=EIF4Efree;
dataS_up.kS(12)=kbR0;
dataS_up.kS(17:157)=kTLNEW_5;
dataS_up.kS(632)=kC173;
dataS_up.kS(541)=kC82;
dataS_up.kS(709)=kA77;
dataS_up.kS(719)=kA87;

% Import DNA Damage parameters
BRCA2=xoutIN(28);
MSH6=xoutIN(29);
MGMT=xoutIN(30);
Me=xoutIN(50);
Ma=xoutIN(59);
fixdsb1=dataS.kS(446);
fixmsh=dataS.kS(447);
fixmgmt=dataS.kS(448);
kDDE=dataS.kS(451);
kDEtop=dataS.kS(452);
Etop=dataS.kS(453); %Going to be zero.
kDnSP=dataS.kS(454);
kDkmSP=dataS.kS(455);
kDkmSS=dataS.kS(458);
kDkmDS=dataS.kS(459);

% Set this parameters manually -- zero-order rate constant
% representing basal DNA damage (both single and double stranded breaks)
kDDbasal=1E-6;


% When cells are cycling and there is no Etoposide, hill function for
% s-phase cyclins will be ~1 and Etoposide equation will be 0.
damageDSB_cycling=kDDbasal/(fixdsb1*BRCA2);
damageSSB_cycling=kDDbasal/(fixmsh*MSH6+fixmgmt*MGMT);

if damageDSB_cycling>kDkmDS
    disp('ERROR --- DSB damage is too high, must reduce kDDbasal or increase strength of repair')
end
if damageSSB_cycling>kDkmSS
    disp('ERROR --- SSB damage is too high, must reduce kDDbasal or increase strength of repair')
end

vdamage_on=(kDDbasal+kDDE*(Etop/(Etop+kDEtop)))*(((Me+Ma)^kDnSP)/(((Me+Ma)^kDnSP)+(kDkmSP^kDnSP))); 
damageDSB=vdamage_on/(fixdsb1*BRCA2);
damageSSB=vdamage_on/(fixmsh*MSH6+fixmgmt*MGMT);

% Run simulation to check and get new steady state
dataS_up.kS(450)=kDDbasal;
xoutIN(31)=damageDSB;
xoutIN(32)=damageSSB;
options=odeset('RelTol',1e-6,'NonNegative',1:length(dataS.x0PARCDL));
[toutNEW6,xoutNEW6]=ode15s(@createODEs,[0 1000*3600],xoutIN,options,dataS_up);


dlmwrite([pathi,'i_xoutNEW6.txt'],xoutNEW6,'Precision',10)
dlmwrite([pathi,'i_toutNEW6.txt'],toutNEW6);
dlmwrite([pathi,'i_kDDbasal.txt'],kDDbasal,'Precision',10);

end

%% Save out Final kTLs and Species Values

if mods(7)
    
disp('...Saving out Final Species Data')
xoutF=dlmread([pathi,'i_xoutNEW6.txt']);
toutF=dlmread([pathi,'i_toutNEW6.txt']);
kTLF=dlmread([pathi,'i_kTLNEW_5.txt']);
EIF4Efinal=xoutF(end,indsEIF4E);

dlmwrite([pathi,'i_EIF4Efree.txt'],EIF4Efinal,'Precision',10);
dlmwrite([pathi,'i_xoutF.txt'],xoutF,'Precision',10)
dlmwrite([pathi,'i_toutF.txt'],toutF);
dlmwrite([pathi,'i_kTLF.txt'],kTLF,'Precision',10);
end

%% Step #5: Transcriptional Activators and Repressors (TARs)

if mods(8)
disp('Step #5: Transcriptional Activators and Repressors (TARs)')

% Defining inputs
xoutF=dlmread([pathi,'i_xoutF.txt']);
numberofTARs=7;

pcFos_cJun=xoutF(end,685); %1
cMyc=xoutF(end,686); %2
p53ac=xoutF(end,3); %3
FOXOnuc=xoutF(end,768); %4
ppERKnuc=xoutF(end,676); %5
pRSKnuc=xoutF(end,679); %6
bCATENINnuc=xoutF(end,687); %7

TAs=zeros(numberofgenes,numberofTARs);
TAs([10:12,99],1)=pcFos_cJun;
TAs(10:12,2)=cMyc;
TAs([26,53,54],3)=p53ac;
TAs([55,58,59,60,61,63,65,66,127,128,136,140],4)=FOXOnuc;
TAs([68,92,97,98],5)=ppERKnuc;
TAs([68,92,97,98],6)=pRSKnuc;
TAs(100,7)=bCATENINnuc;
TAs=TAs*(1/mpc2nmcf_Vn);

TRs=zeros(numberofgenes,numberofTARs);
TRs(98,1)=pcFos_cJun;
TRs=TRs*(1/mpc2nmcf_Vn);

TFa=(TAs./tck50as).^tcnas;
TFa(isnan(TFa))=0;
TFr=(TRs./tck50rs).^tcnrs;
TFr(isnan(TFr))=0;
hills=sum(TFa,2)./(1+sum(TFa,2)+sum(TFr,2));

% AP1*cMYC exception:
hills(10:12)=(TFa(10:12,1)./(1+TFa(10:12,1))).*(TFa(10:12,2)./(1+TFa(10:12,2)));

vTCd=kTCd.*mExp_mpc;
xgac_mpc_D=x0gm_mpc_D(1:numberofgenes);
kTCmaxs=kTCmaxs.*logical(mExp_mpc); 
induction=xgac_mpc_D.*kTCmaxs.*hills;

%Check if induction terms do not exceed degradation terms (which should be
%impossible theoretically)
negativecheck=logical(mExp_mpc).*(vTCd-induction);
i2c=find(negativecheck<0); %inds to change.
if ~isempty(i2c)
    disp('WARNING -- Some induction terms exceed degradation terms')
    disp([i2c,vTCd(i2c)-induction(i2c)])
end

% Calculating leak terms
leak=vTCd-induction;
kTCleak=leak./xgac_mpc_D;
kTCleak(isnan(kTCleak))=0; 
kTCleak(isinf(kTCleak))=0;

% Writing out data
dlmwrite([pathi,'i_kTCleakF.txt'],kTCleak,'Precision',10)
dlmwrite([pathi,'i_kTCmaxsF.txt'],kTCmaxs,'Precision',10)

disp('TARs DONE')

end

end %Initialize




%% Adjunct functions

function [toutNEW,xoutNEW,kTLNEW,dataS_up,flagA] = kTLAdjustWhile(pathi,flagR,flagE,xoutIN,kTLIN,S_TL,dataS,dataS_up)

% This function estimates kTL parameters in an interative approach to 
% equate total protein amounts with RNAseq and proteomics data in a
% serum-starved state.

% Setting up initial parameters
NumObs=102;
indspS6K=725;
indsPARP=104;
indskTL=17:157;

x0PARCDL=dataS.x0PARCDL;
VxPARCDL=dataS.VxPARCDL;
Vc=dataS.kS(3);
apop_def=xoutIN(indsPARP)*.5;
kTLIN_up=kTLIN;



if flagR
    Rt=dataS.kS(1);
    kbRi=dataS.kS(13);
    kdR0=dataS.kS(14);
    nR=dataS.kS(15);
    k50R=dataS.kS(16);
    pS6Ki=xoutIN(end,indspS6K);
    f1=((pS6Ki^nR)/(k50R^nR+pS6Ki^nR));
    if f1>0.1
        disp('WARNING -- Tune pS6K levels down!');
    end
    kbR0=(Rt*kdR0)-(kbRi*f1); 
    dlmwrite([pathi,'i_kbR0.txt'],kbR0)
else
    kbRi=0;
    kdR0=0;
    kbR0=0;
end

% modifying data structure
dataS_up.kS(12)=kbR0;
dataS_up.kS(13)=kbRi;
dataS_up.kS(14)=kdR0;
dataS_up.kS(indskTL)=kTLIN_up;
dataS_up.flagE=flagE;


% Initializing While Loop
tstart=0;
tend=3600*1000;
tspan=[tstart, tend];
th=0.01;
dca=27:NumObs;
options=odeset('RelTol',1e-6,'NonNegative',1:length(x0PARCDL));
xout_update=xoutIN;

ratios_new=ones(NumObs,1)*0.5;
n=1;

% START LOOP
while sum(and(ratios_new(dca)<=1+th,ratios_new(dca)>=1-th))<length(dca) %Until every ratio is between two limits.
    
    
[toutNEW,xoutNEW]=ode15s(@createODEs,tspan,xout_update,options,dataS_up);

[ratios_new,~] = GetTLRatios(x0PARCDL',xoutNEW(end,:),VxPARCDL,Vc); 
ratios_TL=Obs2Genes(S_TL,ratios_new); 
F=1+((ratios_TL-1)/2);

if sum(and(ratios_new(dca)<=1+th,ratios_new(dca)>=1-th))<length(dca)
   kTLIN_up=kTLIN_up.*F;
   dataS_up.kS(indskTL)=kTLIN_up;
   xout_update=xoutNEW(end,:);
else disp('Threshold Reached!')
end

% Check if apoptosis happened
if logical(sum(xoutNEW(:,indsPARP)<apop_def))
    flagA=1; kTLNEW=kTLIN_up; disp('Apoptosis Happened During kTLAdjust'); return %apoptosis happened
else flagA=0; %apoptosis didn't happen
end
    
disp([n,sum(and(ratios_new(dca)<=1+th,ratios_new(dca)>=1-th))])
disp([find(~and(ratios_new(dca)<=1+th,ratios_new(dca)>=1-th))]+26)
n=n+1;
end %While

kTLNEW=kTLIN_up;

end %kTLAdjustWhile


function [ratios,dif] = GetTLRatios(xoutIdeal,xoutIs,VxPARCDL,Vc)

% Get ratios to correct translation rate constants in order to account for
% synthesis and degradation of proteins. Need to run the model to steady
% state and pass into this function.

ObsIdeal = GetObservables(xoutIdeal,VxPARCDL,Vc); ObsIdeal(ObsIdeal<=1E-6)=0;
ObsIs = GetObservables(xoutIs,VxPARCDL,Vc); ObsIs(ObsIs<=1E-6)=0; ObsIs(ObsIs==0)=1;

ratios=ObsIdeal./ObsIs;
ratios(ratios==0)=1;

dif=(ObsIdeal-ObsIs)./ObsIdeal;

end


function Obs = GetObservables(xout,VxPARCDL,Vc)

% Input all species and this will return a sum across of all the species
% that comprise by each observable.

ObsMat=csvread('observables_mat_18.csv',1,1);

for i=1:size(ObsMat,2)
    Obs(i)=sum(ObsMat(:,i).*xout'.*(VxPARCDL/Vc));
end

end


function ratios_TL=Obs2Genes(S_TL,ratios_TLC)

ratios_TL=ones(size(S_TL,2),1);
S_TL_obs = S_TL(find(sum(S_TL,2)),:);
for i=1:size(S_TL_obs,2)
    ratio2add=ratios_TLC(find(S_TL_obs(:,i)));
    if ratio2add
        ratios_TL(i)=ratios_TLC(find(S_TL_obs(:,i)));
    else ratios_TL(i)=1;
    end
end

end

