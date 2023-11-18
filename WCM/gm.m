function [xginN,xgacN,AllGenesVecN,xmN,vTC] = gm(flagD,dataG,ts,xoutG,xoutS)

Vn=dataG.Vn;
mpc2nmcf_Vn=1E9/(Vn*6.023E+23);

% Defining Terms
kGin=dataG.kGin;
kGac=dataG.kGac;
kTCleak=dataG.kTCleak;
kTCmaxs=dataG.kTCmaxs;
kTCd=dataG.kTCd;
tcnas=dataG.tcnas;
tcnrs=dataG.tcnrs;
tck50as=dataG.tck50as;
tck50rs=dataG.tck50rs;
GenePositionMatrix=dataG.GenePositionMatrix;
AllGenesVec=dataG.AllGenesVec;
indsD=dataG.indsD;

numberofgenes=size(tck50as,1);
numberofTARs=size(tck50as,2);

% gm species
a=0;
xgac=xoutG(a+1:a+numberofgenes);a=a+numberofgenes;
xgin=xoutG(a+1:a+numberofgenes);a=a+numberofgenes;
xm=xoutG(a+1:a+numberofgenes);

% Gene switching constants
kGin_1=kGin(1);
kGac_1=kGac(1);



%% vTC and vTCd

TAs=zeros(numberofgenes,numberofTARs);
TRs=zeros(numberofgenes,numberofTARs);


% TARs
pcFos_cJun=xoutS(685); %1
cMyc=xoutS(686); %2
p53ac=xoutS(3); %3
FOXOnuc=xoutS(768); %4
ppERKnuc=xoutS(676); %5
pRSKnuc=xoutS(679); %6
bCATENINnuc=xoutS(687); %7
% activators
TAs([10:12,99],1)=pcFos_cJun;
TAs(10:12,2)=cMyc;
TAs([26,53,54],3)=p53ac;
TAs([55,58,59,60,61,63,65,66,127,128,136,140],4)=FOXOnuc;
TAs([68,92,97,98],5)=ppERKnuc;
TAs([68,92,97,98],6)=pRSKnuc;
TAs(100,7)=bCATENINnuc;
TAs=TAs*(1/mpc2nmcf_Vn);
% repressors
TRs(98,1)=pcFos_cJun;
TRs=TRs*(1/mpc2nmcf_Vn);
% make hills
TFa=(TAs./tck50as).^tcnas;
TFa(isnan(TFa))=0;
TFr=(TRs./tck50rs).^tcnrs;
TFr(isnan(TFr))=0;
hills=sum(TFa,2)./(1+sum(TFa,2)+sum(TFr,2));

% With AP1*cMYC exception:
hills(10:12)=(TFa(10:12,1)./(1+TFa(10:12,1))).*(TFa(10:12,2)./(1+TFa(10:12,2)));


% vTC
induced=xgac.*kTCmaxs.*hills;
leak=xgac.*kTCleak;
vTC=leak+induced;

% vTCd
vTCd=kTCd.*xm;


if ~isempty(AllGenesVec)
else
    xginN=[];
    xgacN=[];
    AllGenesVecN=[];
    xmN=[];  
    return
end


%% Poisson Stuff

poff=poisspdf(0,kGin_1*ts); 
pon=poisspdf(0,kGac_1*ts); 

% Generating random numbers and deciding which genes should turn off and on
RandomNumbers=rand(length(AllGenesVec),1);
geneson=logical(AllGenesVec);
genesoff=logical(~AllGenesVec);
ac2in=and(geneson,RandomNumbers>=poff); 
in2ac=and(genesoff,RandomNumbers>=pon); 

% Generating new AllGenesVec and allocating active and inactive genes
AllGenesVecN=AllGenesVec;
AllGenesVecN(ac2in)=0;
AllGenesVecN(in2ac)=1;
xgacN=GenePositionMatrix*AllGenesVecN;
xginN=(xgac+xgin)-xgacN;


% mRNA
Nb=poissrnd(vTC*ts);
Nd=poissrnd(vTCd*ts);

% These genes and mRNAs we don't allow to fluctuate
Nb(indsD)=vTC(indsD)*ts;
Nd(indsD)=vTCd(indsD)*ts;
xgacN(indsD)=xoutG(indsD);
xginN(indsD)=xoutG(indsD+numberofgenes);

% OUTPUT deterministic results instead:
if flagD
    Nb=vTC*ts;
    Nd=vTCd*ts;
    xgacN=xoutG(1:numberofgenes);
    xginN=xoutG(numberofgenes+1:numberofgenes*2);
end

% Finish mRNA
xmN=xm+Nb-Nd;
xmN(xmN<0)=0;





