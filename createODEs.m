function [ydot,flag,new_data,v] = createODEs(t,x,data)


% unpacking
kS=data.kS;
VvPARCDL=data.VvPARCDL;
VxPARCDL=data.VxPARCDL;
S_PARCDL=data.S_PARCDL;
mExp_nM=data.mExp_nM;
mMod=data.mMod;
flagE=data.flagE;

% k
a=0;
u=1;Rt=kS(a+1:a+u);a=a+u;
u=1;EIF4Efree=kS(a+1:a+u);a=a+u;
u=1;Vc=kS(a+1:a+u);a=a+u;
u=1;Ve=kS(a+1:a+u);a=a+u;
u=1;Vm=kS(a+1:a+u);a=a+u;
u=1;Vn=kS(a+1:a+u);a=a+u;
u=1;kT1=kS(a+1:a+u);a=a+u;
u=1;kT2=kS(a+1:a+u);a=a+u;
u=1;kT3=kS(a+1:a+u);a=a+u;
u=1;kT4=kS(a+1:a+u);a=a+u;
u=1;k50E=kS(a+1:a+u);a=a+u;
u=1;kbR0=kS(a+1:a+u);a=a+u;
u=1;kbRi=kS(a+1:a+u);a=a+u;
u=1;kdR0=kS(a+1:a+u);a=a+u;
u=1;nR=kS(a+1:a+u);a=a+u;
u=1;k50R=kS(a+1:a+u);a=a+u;
u=141;kTL=kS(a+1:a+u);a=a+u;
u=141;kTLd=kS(a+1:a+u);a=a+u;
u=102;kTLCd=kS(a+1:a+u);a=a+u;
u=59;kD=kS(a+1:a+u);a=a+u;
u=173;kC=kS(a+1:a+u);a=a+u;
u=87;kA=kS(a+1:a+u);a=a+u;
u=158;kR=kS(a+1:a+u);a=a+u;
u=190;kP=kS(a+1:a+u);a=a+u;
u=29;kRP1=kS(a+1:a+u);a=a+u;
u=29;kRP2=kS(a+1:a+u);a=a+u;
u=29;kRP3=kS(a+1:a+u);a=a+u;
u=29;kRP4=kS(a+1:a+u);a=a+u;
u=29;kRP5=kS(a+1:a+u);a=a+u;
u=29;kRP6=kS(a+1:a+u);a=a+u;
u=29;kRP7=kS(a+1:a+u);a=a+u;
u=29;kRP8=kS(a+1:a+u);a=a+u;
u=29;kRP9=kS(a+1:a+u);a=a+u;
u=29;kRP10=kS(a+1:a+u);a=a+u;
u=29;kRP11=kS(a+1:a+u);a=a+u;
u=29;kRP12=kS(a+1:a+u);a=a+u;
u=29;kRP13=kS(a+1:a+u);a=a+u;
u=29;kRP14=kS(a+1:a+u);a=a+u;
u=29;kRP15=kS(a+1:a+u);a=a+u;
u=29;kRP16=kS(a+1:a+u);a=a+u;
u=29;kRP17=kS(a+1:a+u);a=a+u;
u=29;kRP18=kS(a+1:a+u);a=a+u;
u=29;kRP19=kS(a+1:a+u);a=a+u;
u=29;kRP20=kS(a+1:a+u);a=a+u;
u=29;kRP21=kS(a+1:a+u);a=a+u;
u=29;kRP22=kS(a+1:a+u);a=a+u;
u=29;kRP23=kS(a+1:a+u);a=a+u;
u=29;kRP24=kS(a+1:a+u);a=a+u;
u=29;kRP25=kS(a+1:a+u);a=a+u;
u=29;kRP26=kS(a+1:a+u);a=a+u;
u=29;kRP27=kS(a+1:a+u);a=a+u;
u=29;kRP28=kS(a+1:a+u);a=a+u;
u=29;kRP29=kS(a+1:a+u);a=a+u;
u=29;kRP30=kS(a+1:a+u);a=a+u;
u=29;kRP31=kS(a+1:a+u);a=a+u;
u=29;kRP32=kS(a+1:a+u);a=a+u;
u=29;kRP33=kS(a+1:a+u);a=a+u;
u=16;kRP34=kS(a+1:a+u);a=a+u;
u=4;kDP=kS(a+1:a+u);a=a+u;
u=11;kPA=kS(a+1:a+u);a=a+u;
u=606;kXd=kS(a+1:a+u);a=a+u;
u=3;kE=kS(a+1:a+u);


    
%% Defining Entities

a=0;
xRibosome=x(a+1);a=a+1;
xD=x(a+1:a+38);a=a+38;
xC=x(a+1:a+44);a=a+44;
xA=x(a+1:a+72);a=a+72;
xR=x(a+1:a+73);a=a+73; 
xRP=x(a+1:a+414);a=a+414;
xP=x(a+1:a+130);a=a+130; 
xE=x(a+1:a+2);



%% xD **
p53inac =xD(1);
p53ac =xD(2);
Mdm2 =xD(3);
Wip1 =xD(4);
ATMP =xD(5);
ATRac =xD(6);
Mdm2product1=xD(7);
Mdm2product2=xD(8);
Mdm2product3=xD(9);
Mdm2product4=xD(10);
Mdm2product5=xD(11);
Mdm2product6=xD(12);
Mdm2product7=xD(13);
Mdm2product8=xD(14);
Mdm2product9=xD(15);
Mdm2pro =xD(16);
Wip1product1=xD(17);
Wip1product2=xD(18);
Wip1product3=xD(19);
Wip1product4=xD(20);
Wip1product5=xD(21);
Wip1product6=xD(22);
Wip1product7=xD(23);
Wip1product8=xD(24);
Wip1product9=xD(25);
Wip1pro =xD(26);
BRCA2=xD(27);
MSH6=xD(28);
MGMT=xD(29);
damageDSB=xD(30);
damageSSB=xD(31);
ppAKT_Mdm2=xD(32);
pMdm2=xD(33);
ARF=xD(34);
MDM4=xD(35);
p53ac_MDM4=xD(36);
ATMinac=xD(37);
ATRinac=xD(38);


%% xC **
pRB=xC(1);
pRBp=xC(2);
pRBpp=xC(3);
E2F=xC(4);
Cd=xC(5);
Mdi=xC(6);
Md=xC(7);
Mdp27=xC(8);
Ce=xC(9);
Mei=xC(10);
Me=xC(11);
Skp2=xC(12);
Mep27=xC(13);
Pe=xC(14);
Pai=xC(15);
Pei=xC(16);
Pbi=xC(17);
Ca=xC(18);
Mai=xC(19);
Ma=xC(20);
Map27=xC(21);
p27=xC(22);
Cdh1i=xC(23);
Cdh1a=xC(24);
E2Fp=xC(25);
p27p=xC(26);
Pa=xC(27);
Cb=xC(28);
Mbi=xC(29);
Mb=xC(30);
Cdc20i=xC(31);
Cdc20a=xC(32);
Pb=xC(33);
Wee1=xC(34);
Wee1p=xC(35);
Mbp27=xC(36);
Chk1=xC(37);
pRBc1=xC(38);
pRBc2=xC(39);
p21=xC(40);
Mdp21=xC(41);
Mep21=xC(42);
Map21=xC(43);
Mbp21=xC(44);

%% xA **
L =xA(1);
R =xA(2);
L_R =xA(3);
Ractive =xA(4);
flip =xA(5);
Ractive_flip =xA(6);
pC8 =xA(7);
Ractive_pC8 =xA(8);
C8 =xA(9);
Bar =xA(10);
C8_Bar =xA(11);
pC3 =xA(12);
C8_pC3 =xA(13);
C3 =xA(14);
pC6 =xA(15);
C3_pC6 =xA(16);
C6 =xA(17);
C6_C8 =xA(18);
XIAP =xA(19);
C3_XIAP =xA(20);
PARP =xA(21);
C3_PARP =xA(22);
cPARP =xA(23);
Bid =xA(24);
C8_Bid =xA(25);
tBid =xA(26);
Bcl2c =xA(27);
tBid_Bcl2c =xA(28);
Bax =xA(29);
tBid_Bax =xA(30);
Baxactive =xA(31);
Baxm =xA(32);
Bcl2 =xA(33);
Baxm_Bcl2 =xA(34);
Bax2 =xA(35);
Bax2_Bcl2 =xA(36);
Bax4 =xA(37);
Bax4_Bcl2 =xA(38);
M =xA(39);
Bax4_M =xA(40);
Mactive =xA(41);
CytoCm =xA(42);
Mactive_CytoCm =xA(43);
CytoCr =xA(44);
Smacm =xA(45);
Mactive_Smacm =xA(46);
Smacr =xA(47);
CytoC =xA(48);
Apaf =xA(49);
CytoC_Apaf =xA(50);
Apafactive =xA(51);
pC9 =xA(52);
Apop =xA(53);
Apop_C3 =xA(54);
Smac =xA(55);
Apop_XIAP =xA(56);
Smac_XIAP =xA(57);
C3_Ub =xA(58);
BAD=xA(59);
PUMA=xA(60);
NOXA=xA(61);
Bcl2c_BAD=xA(62);
Bcl2c_PUMA=xA(63);
Bcl2c_NOXA=xA(64);
BIM=xA(65);
BIM_Bax=xA(66);
Bcl2c_BIM=xA(67);
ppERK_BIM=xA(68);
pBIM=xA(69);
ppAKT_BAD=xA(70);
pBAD=xA(71);
ppERK_BAD=xA(72);

%% xR **
E=xR(1);
H=xR(2);
HGF=xR(3);
P=xR(4);
F=xR(5);
I=xR(6);
IN=xR(7);
E1=xR(8);
pE1=xR(9);
E2=xR(10);
pE2=xR(11);
E3=xR(12);
E4=xR(13);
pE4=xR(14);
Ev3=xR(15);
Met=xR(16);
Pr=xR(17);
Fr=xR(18);
Ir=xR(19);
Isr =xR(20);
E1E1=xR(21);
E1E2=xR(22);
E1E3=xR(23);
E1E4=xR(24);
E2E2=xR(25);
E2E3=xR(26);
E2E4=xR(27);
E3E4=xR(28);
E4E4=xR(29);
Met_Met=xR(30);
FrFr=xR(31);
IrIr=xR(32);
Isr_Isr=xR(33);
EE1=xR(34);
HE3=xR(35);
HE4=xR(36);
HGF_Met=xR(37);
PPr=xR(38);
FFr=xR(39);
EE1E2=xR(40);
EE1Ev3=xR(41);
EE1E1=xR(42);
EE1E3=xR(43);
EE1E4=xR(44);
E2HE3=xR(45);
E1HE3=xR(46);
HE3E3=xR(47);
HE3Ev3=xR(48);
HE3E4=xR(49);
E2HE4=xR(50);
HE4Ev3=xR(51);
E1HE4=xR(52);
E3HE4=xR(53);
HE4E4=xR(54);
HGF_Met_Met=xR(55);
PPrPr=xR(56);
FFrFr=xR(57);
IIrIr=xR(58);
IN_Isr_Isr=xR(59);
EE1EE1=xR(60);
EE1HE3=xR(61);
EE1HE4=xR(62);
HE3HE3=xR(63);
HE3HE4=xR(64);
HE4HE4=xR(65);
HGF_Met_HGF_Met=xR(66);
PPrPPr=xR(67);
FFrFFr=xR(68);
IIrIrI=xR(69);
IN_Isr_Isr_IN=xR(70);
E1_ppERK=xR(71);
E2_ppERK=xR(72);
E4_ppERK=xR(73);

%% xRP **
a=0;

SCD = [EE1E2;EE1Ev3;EE1E1;EE1EE1;EE1E3;EE1HE3;EE1E4;EE1HE4;E2HE3;HE3Ev3;E1HE3;HE3E4;HE3HE4;E2HE4;HE4Ev3;E1HE4;E3HE4;HE4E4;HE4HE4;HGF_Met_Met;HGF_Met_HGF_Met;PPrPPr;PPrPr;FFrFFr;FFrFr;IIrIr;IN_Isr_Isr;IIrIrI;IN_Isr_Isr_IN];
pSCD = xRP(a+1:a+29); a=a+33;
Sp_SCD = xRP(a+1:a+29); a=a+29;%[Sp_EE1E2,Sp_EE1Ev3,Sp_EE1E1,Sp_EE1EE1,Sp_EE1E3,Sp_EE1HE3,Sp_EE1E4,Sp_EE1HE4,Sp_E2HE3,Sp_HE3Ev3,Sp_E1HE3,Sp_HE3E4,Sp_HE3HE4,Sp_E2HE4,Sp_HE4Ev3,Sp_E1HE4,Sp_E3HE4,Sp_HE4E4,Sp_HE4HE4,Sp_HGF_Met_Met,Sp_HGF_Met_HGF_Met,Sp_PPrPPr,Sp_PPrPr,Sp_FFrFFr,Sp_FFrFr,Sp_IIrIr,Sp_IIrIIr];
SCDint = xRP(a+1:a+29); a=a+29;%[EE1E2int,EE1Ev3int,EE1E1int,EE1EE1int,EE1E3int,EE1HE3int,EE1E4int,EE1HE4int,E2HE3int,HE3Ev3int,E1HE3int,HE3E4int,HE3HE4int,E2HE4int,HE4Ev3int,E1HE4int,E3HE4int,HE4E4int,HE4HE4int,HGF_Met_Metint,HGF_Met_HGF_Metint,PPrPPrint,PPrPrint,FFrFFrint,FFrFrint,IIrIrint,IIrIIrint];
pSCDint = xRP(a+1:a+29); a=a+33;%[pEE1E2int,pEE1Ev3int,pEE1E1int,pEE1EE1int,pEE1E3int,pEE1HE3int,pEE1E4int,pEE1HE4int,pE2HE3int,pHE3Ev3int,pE1HE3int,pHE3E4int,pHE3HE4int,pE2HE4int,pHE4Ev3int,pE1HE4int,pE3HE4int,pHE4E4int,pHE4HE4int,pHGF_Met_Metint,pHGF_Met_HGF_Metint,pPPrPPrint,pPPrPrint,pFFrFFrint,pFFrFrint,pIIrIrint,pIIrIIrint];

pSCD_bind = xRP([1:25,30:33]); %[pEE1E2,pEE1Ev3,pEE1E1,pEE1EE1,pEE1E3,pEE1HE3,pEE1E4,pEE1HE4,pE2HE3,pHE3Ev3,pE1HE3,pHE3E4,pHE3HE4,pE2HE4,pHE4Ev3,pE1HE4,pE3HE4,pHE4E4,pHE4HE4,pHGF_Met_Met,pHGF_Met_HGF_Met,pPPrPPr,pPPrPr,pFFrFFr,pFFrFr,pIIrIr_IRS,pIIrIIr_IRS]; %Same as pSCD but with the IRS on the IGF1R. Had to do this because it has an extra step (binding of IRS)
pSCDint_bind = xRP([92:116,121:124]); %[pEE1E2int,pEE1Ev3int,pEE1E1int,pEE1EE1int,pEE1E3int,pEE1HE3int,pEE1E4int,pEE1HE4int,pE2HE3int,pHE3Ev3int,pE1HE3int,pHE3E4int,pHE3HE4int,pE2HE4int,pHE4Ev3int,pE1HE4int,pE3HE4int,pHE4E4int,pHE4HE4int,pHGF_Met_Metint,pHGF_Met_HGF_Metint,pPPrPPrint,pPPrPrint,pFFrFFrint,pFFrFrint,pIIrIrint_IRS,pIIrIIrint_IRS];

pSCD_G2_SOS = xRP(a+1:a+29); a=a+29; %[pEE1E2_G2_SOS,pEE1Ev3_G2_SOS,pEE1E1_G2_SOS,pEE1EE1_G2_SOS,pEE1E3_G2_SOS,pEE1HE3_G2_SOS,pEE1E4_G2_SOS,pEE1HE4_G2_SOS,pE2HE3_G2_SOS,pHE3Ev3_G2_SOS,pE1HE3_G2_SOS,pHE3E4_G2_SOS,pHE3HE4_G2_SOS,pE2HE4_G2_SOS,pHE4Ev3_G2_SOS,pE1HE4_G2_SOS,pE3HE4_G2_SOS,pHE4E4_G2_SOS,pHE4HE4_G2_SOS,pHGF_Met_Met_G2_SOS,pHGF_Met_HGF_Met_G2_SOS,pPPrPPr_G2_SOS,pPPrPr_G2_SOS,pFFrFFr_G2_SOS,pFFrFr_G2_SOS,pIIrIr_IRS_G2_SOS,pIIrIIr_IRS_G2_SOS];
pSCDint_G2_SOS = xRP(a+1:a+29); a=a+29; %[pEE1E2int_G2_SOS,pEE1Ev3int_G2_SOS,pEE1E1int_G2_SOS,pEE1EE1int_G2_SOS,pEE1E3int_G2_SOS,pEE1HE3int_G2_SOS,pEE1E4int_G2_SOS,pEE1HE4int_G2_SOS,pE2HE3int_G2_SOS,pHE3Ev3int_G2_SOS,pE1HE3int_G2_SOS,pHE3E4int_G2_SOS,pHE3HE4int_G2_SOS,pE2HE4int_G2_SOS,pHE4Ev3int_G2_SOS,pE1HE4int_G2_SOS,pE3HE4int_G2_SOS,pHE4E4int_G2_SOS,pHE4HE4int_G2_SOS,pHGF_Met_Metint_G2_SOS,pHGF_Met_HGF_Metint_G2_SOS,pPPrPPrint_G2_SOS,pPPrPrint_G2_SOS,pFFrFFrint_G2_SOS,pFFrFrint_G2_SOS,pIIrIrint_IRS_G2_SOS,pIIrIIrint_IRS_G2_SOS];
pSCD_PLCg = xRP(a+1:a+29); a=a+29; %[pEE1E2_PLCg,pEE1Ev3_PLCg,pEE1E1_PLCg,pEE1EE1_PLCg,pEE1E3_PLCg,pEE1HE3_PLCg,pEE1E4_PLCg,pEE1HE4_PLCg,pE2HE3_PLCg,pHE3Ev3_PLCg,pE1HE3_PLCg,pHE3E4_PLCg,pHE3HE4_PLCg,pE2HE4_PLCg,pHE4Ev3_PLCg,pE1HE4_PLCg,pE3HE4_PLCg,pHE4E4_PLCg,pHE4HE4_PLCg,pHGF_Met_Met_PLCg,pHGF_Met_HGF_Met_PLCg,pPPrPPr_PLCg,pPPrPr_PLCg,pFFrFFr_PLCg,pFFrFr_PLCg,pIIrIr_IRS_PLCg,pIIrIIr_IRS_PLCg];
pSCD_PI3K1 = xRP(a+1:a+29); a=a+29; %[pEE1E2_PI3K1,pEE1Ev3_PI3K1,pEE1E1_PI3K1,pEE1EE1_PI3K1,pEE1E3_PI3K1,pEE1HE3_PI3K1,pEE1E4_PI3K1,pEE1HE4_PI3K1,pE2HE3_PI3K1,pHE3Ev3_PI3K1,pE1HE3_PI3K1,pHE3E4_PI3K1,pHE3HE4_PI3K1,pE2HE4_PI3K1,pHE4Ev3_PI3K1,pE1HE4_PI3K1,pE3HE4_PI3K1,pHE4E4_PI3K1,pHE4HE4_PI3K1,pHGF_Met_Met_PI3K1,pHGF_Met_HGF_Met_PI3K1,pPPrPPr_PI3K1,pPPrPr_PI3K1,pFFrFFr_PI3K1,pFFrFr_PI3K1,pIIrIr_IRS_PI3K1,pIIrIIr_IRS_PI3K1];
pSCD_PI3K2 = xRP(a+1:a+29); a=a+29; %[pEE1E2_PI3K2,pEE1Ev3_PI3K2,pEE1E1_PI3K2,pEE1EE1_PI3K2,pEE1E3_PI3K2,pEE1HE3_PI3K2,pEE1E4_PI3K2,pEE1HE4_PI3K2,pE2HE3_PI3K2,pHE3Ev3_PI3K2,pE1HE3_PI3K2,pHE3E4_PI3K2,pHE3HE4_PI3K2,pE2HE4_PI3K2,pHE4Ev3_PI3K2,pE1HE4_PI3K2,pE3HE4_PI3K2,pHE4E4_PI3K2,pHE4HE4_PI3K2,pHGF_Met_Met_PI3K2,pHGF_Met_HGF_Met_PI3K2,pPPrPPr_PI3K2,pPPrPr_PI3K2,pFFrFFr_PI3K2,pFFrFr_PI3K2,pIIrIr_IRS_PI3K2,pIIrIIr_IRS_PI3K2];

pSCDint_G2_SOS_RasD = xRP(a+1:a+29); a=a+29; %[pEE1E2int_G2_SOS_RasDint,pEE1Ev3int_G2_SOS_RasDint,pEE1E1int_G2_SOS_RasDint,pEE1EE1int_G2_SOS_RasDint,pEE1E3int_G2_SOS_RasDint,pEE1HE3int_G2_SOS_RasDint,pEE1E4int_G2_SOS_RasDint,pEE1HE4int_G2_SOS_RasDint,pE2HE3int_G2_SOS_RasDint,pHE3Ev3int_G2_SOS_RasDint,pE1HE3int_G2_SOS_RasDint,pHE3E4int_G2_SOS_RasDint,pHE3HE4int_G2_SOS_RasDint,pE2HE4int_G2_SOS_RasDint,pHE4Ev3int_G2_SOS_RasDint,pE1HE4int_G2_SOS_RasDint,pE3HE4int_G2_SOS_RasDint,pHE4E4int_G2_SOS_RasDint,pHE4HE4int_G2_SOS_RasDint,pHGF_Met_Metint_G2_SOS_RasDint,pHGF_Met_HGF_Metint_G2_SOS_RasDint,pPPrPPrint_G2_SOS_RasDint,pPPrPrint_G2_SOS_RasDint,pFFrFFrint_G2_SOS_RasDint,pFFrFrint_G2_SOS_RasDint,pIIrIrint_IRS_G2_SOS_RasDint,pIIrIIrint_IRS_G2_SOS_RasDint];
pSCD_G2_SOS_RasD = xRP(a+1:a+29); a=a+29; %[pEE1E2_G2_SOS_RasD,pEE1Ev3_G2_SOS_RasD,pEE1E1_G2_SOS_RasD,pEE1EE1_G2_SOS_RasD,pEE1E3_G2_SOS_RasD,pEE1HE3_G2_SOS_RasD,pEE1E4_G2_SOS_RasD,pEE1HE4_G2_SOS_RasD,pE2HE3_G2_SOS_RasD,pHE3Ev3_G2_SOS_RasD,pE1HE3_G2_SOS_RasD,pHE3E4_G2_SOS_RasD,pHE3HE4_G2_SOS_RasD,pE2HE4_G2_SOS_RasD,pHE4Ev3_G2_SOS_RasD,pE1HE4_G2_SOS_RasD,pE3HE4_G2_SOS_RasD,pHE4E4_G2_SOS_RasD,pHE4HE4_G2_SOS_RasD,pHGF_Met_Met_G2_SOS_RasD,pHGF_Met_HGF_Met_G2_SOS_RasD,pPPrPPr_G2_SOS_RasD,pPPrPr_G2_SOS_RasD,pFFrFFr_G2_SOS_RasD,pFFrFr_G2_SOS_RasD,pIIrIr_IRS_G2_SOS_RasD,pIIrIIr_IRS_G2_SOS_RasD];
pSCD_PLCg_PIP2 = xRP(a+1:a+29); a=a+29; %[pEE1E2_PLCg_PIP2,pEE1Ev3_PLCg_PIP2,pEE1E1_PLCg_PIP2,pEE1EE1_PLCg_PIP2,pEE1E3_PLCg_PIP2,pEE1HE3_PLCg_PIP2,pEE1E4_PLCg_PIP2,pEE1HE4_PLCg_PIP2,pE2HE3_PLCg_PIP2,pHE3Ev3_PLCg_PIP2,pE1HE3_PLCg_PIP2,pHE3E4_PLCg_PIP2,pHE3HE4_PLCg_PIP2,pE2HE4_PLCg_PIP2,pHE4Ev3_PLCg_PIP2,pE1HE4_PLCg_PIP2,pE3HE4_PLCg_PIP2,pHE4E4_PLCg_PIP2,pHE4HE4_PLCg_PIP2,pHGF_Met_Met_PLCg_PIP2,pHGF_Met_HGF_Met_PLCg_PIP2,pPPrPPr_PLCg_PIP2,pPPrPr_PLCg_PIP2,pFFrFFr_PLCg_PIP2,pFFrFr_PLCg_PIP2,pIIrIr_IRS_PLCg_PIP2,pIIrIIr_IRS_PLCg_PIP2];
pSCD_PI3K1_PIP2 = xRP(a+1:a+29); a=a+29; %[pEE1E2_PI3K1_PIP2,pEE1Ev3_PI3K1_PIP2,pEE1E1_PI3K1_PIP2,pEE1EE1_PI3K1_PIP2,pEE1E3_PI3K1_PIP2,pEE1HE3_PI3K1_PIP2,pEE1E4_PI3K1_PIP2,pEE1HE4_PI3K1_PIP2,pE2HE3_PI3K1_PIP2,pHE3Ev3_PI3K1_PIP2,pE1HE3_PI3K1_PIP2,pHE3E4_PI3K1_PIP2,pHE3HE4_PI3K1_PIP2,pE2HE4_PI3K1_PIP2,pHE4Ev3_PI3K1_PIP2,pE1HE4_PI3K1_PIP2,pE3HE4_PI3K1_PIP2,pHE4E4_PI3K1_PIP2,pHE4HE4_PI3K1_PIP2,pHGF_Met_Met_PI3K1_PIP2,pHGF_Met_HGF_Met_PI3K1_PIP2,pPPrPPr_PI3K1_PIP2,pPPrPr_PI3K1_PIP2,pFFrFFr_PI3K1_PIP2,pFFrFr_PI3K1_PIP2,pIIrIr_IRS_PI3K1_PIP2,pIIrIIr_IRS_PI3K1_PIP2];
pSCD_PI3K2_PIP = xRP(a+1:a+29); %[pEE1E2_PI3K2_PIP,pEE1Ev3_PI3K2_PIP,pEE1E1_PI3K2_PIP,pEE1EE1_PI3K2_PIP,pEE1E3_PI3K2_PIP,pEE1HE3_PI3K2_PIP,pEE1E4_PI3K2_PIP,pEE1HE4_PI3K2_PIP,pE2HE3_PI3K2_PIP,pHE3Ev3_PI3K2_PIP,pE1HE3_PI3K2_PIP,pHE3E4_PI3K2_PIP,pHE3HE4_PI3K2_PIP,pE2HE4_PI3K2_PIP,pHE4Ev3_PI3K2_PIP,pE1HE4_PI3K2_PIP,pE3HE4_PI3K2_PIP,pHE4E4_PI3K2_PIP,pHE4HE4_PI3K2_PIP,pHGF_Met_Met_PI3K2_PIP,pHGF_Met_HGF_Met_PI3K2_PIP,pPPrPPr_PI3K2_PIP,pPPrPr_PI3K2_PIP,pFFrFFr_PI3K2_PIP,pFFrFr_PI3K2_PIP,pIIrIr_IRS_PI3K2_PIP,pIIrIIr_IRS_PI3K2_PIP];

% Need these independently because of binding to IRS
pIIrIr=xRP(26);
pIN_Isr_Isr=xRP(27);
pIIrIrI=xRP(28);
pIN_Isr_Isr_IN=xRP(29);
pIIrIr_IRS=xRP(30);
pIN_Isr_Isr_IRS=xRP(31);
pIIrIrI_IRS=xRP(32);
pIN_Isr_Isr_IN_IRS=xRP(33);

pIIrIr_int=xRP(117);
pIN_Isr_Isr_int=xRP(118);
pIIrIrI_int=xRP(119);
pIN_Isr_Isr_IN_int=xRP(120);
pIIrIr_int_IRS=xRP(121);
pIN_Isr_Isr_int_IRS=xRP(122);
pIIrIrI_int_IRS=xRP(123);
pIN_Isr_Isr_IN_int_IRS=xRP(124);



%% xP **
IRS=xP(1);
Sp=xP(2);
Cbl=xP(3);
G2=xP(4);
G2_SOS=xP(5);
G2_pSOS=xP(6);
PLCg=xP(7);
PI3KC1=xP(8);
PI3KR1=xP(9);
PI3K1=xP(10);
pPI3K1=xP(11);
PI3K2=xP(12);
mTORC1=xP(13);
mTORC1active=xP(14);
PIP=xP(15);
PI3P=xP(16);
DAG=xP(17);
GRP=xP(18);
DAG_GRP=xP(19);
RasT=xP(20);
RasD=xP(21);
NF1=xP(22);
pNF1=xP(23);
pCRaf=xP(24);
CRaf=xP(25);
RasT_CRaf=xP(26);
BRaf=xP(27);
RasT_CRaf_BRaf=xP(28);
MEK=xP(29);
pMEK=xP(30);
ppMEK=xP(31);
MKP3=xP(32);
ERKnuc=xP(33);
ppERKnuc=xP(34);
RSK=xP(35);
pRSK=xP(36);
pRSKnuc=xP(37);
MKP1=xP(38);
pMKP1=xP(39);
cFos=xP(40);
pcFos=xP(41);
cJun=xP(42);
pcFos_cJun=xP(43);
cMyc=xP(44);
bCATENINnuc=xP(45);
bCATENIN=xP(46);
pbCATENIN=xP(47);
IP3=xP(48);
PIP2=xP(49);
PIP3=xP(50);
PTEN=xP(51);
PIP3_AKT=xP(52);
AKT=xP(53);
pAKT=xP(54);
ppAKT=xP(55);
PDK1=xP(56);
PIP3_PDK1=xP(57);
PIP3_pAKT=xP(58);
Rictor=xP(59);
mTOR=xP(60);
mTORC2=xP(61);
PIP3_ppAKT=xP(62);
GSK3b=xP(63);
pGSK3b=xP(64);
TSC1=xP(65);
TSC2=xP(66);
pTSC2=xP(67);
TSC=xP(68);
PKC=xP(69);
DAG_PKC=xP(70);
pRKIP=xP(71);
RKIP=xP(72);
RKIP_CRaf=xP(73);
ERK=xP(74);
pERK=xP(75);
ppERK=xP(76);
FOXO=xP(77);
pFOXO=xP(78);
RhebD=xP(79);
RhebT=xP(80);
Raptor=xP(81);
S6K=xP(82);
pS6K=xP(83);
EIF4EBP1=xP(84);
pEIF4EBP1=xP(85);
SOS=xP(86);
G2_SOS_ppERK=xP(87);
CRaf_ppERK=xP(88);
RasD_DAG_GRP=xP(89);
RasT_NF1=xP(90);
NF1_ppERK=xP(91);
MEK_RasT_CRaf_BRaf=xP(92);
pMEK_RasT_CRaf_BRaf=xP(93);
ERK_ppMEK=xP(94);
pERK_ppMEK=xP(95);
RSK_ppERK=xP(96);
pRSKnuc_MKP1=xP(97);
ppERKnuc_MKP1=xP(98);
cFos_pRSKnuc=xP(99);
cFos_ppERKnuc=xP(100);
RKIP_DAG_PKC=xP(101);
PIP3_PTEN=xP(102);
PIP3_AKT_PIP3_PDK1=xP(103);
PIP3_pAKT_mTORC2=xP(104);
GSK3b_ppAKT=xP(105);
bCATENIN_GSK3b=xP(106);
TSC2_ppAKT=xP(107);
TSC2_ppERK=xP(108);
RhebT_TSC=xP(109);
EIF4EBP1_mTORC1active=xP(110);
S6K_mTORC1active=xP(111);
FOXO_ppAKT=xP(112);
PI3K1_mTORC1active=xP(113);
pERK_MKP3=xP(114);
ppERK_MKP3=xP(115);
ppERKnuc_pMKP1=xP(116);
RasT_BRaf=xP(117);
RasT_BRaf_BRaf=xP(118);
MEK_RasT_BRaf_BRaf=xP(119);
pMEK_RasT_BRaf_BRaf=xP(120);
EIF4E=xP(121);
EIF4EBP1_EIF4E=xP(122);
RasT_CRaf_CRaf=xP(123);
MEK_RasT_CRaf_CRaf=xP(124);
pMEK_RasT_CRaf_CRaf=xP(125);
FOXOnuc=xP(126);

MEKi=xP(127);
MEKi_ppMEK=xP(128);
AKTi=xP(129);
AKTi_AKT=xP(130);




%%%%%% RATE CONSTANTS ASSIGNMENT

%%%DNA Damage Rate Constants
bp=kD(1);
ampi=kD(2);
api=kD(3);
bsp=kD(4);
ns=kD(5);
Ts=kD(6);
ampa=kD(7);
awpa=kD(8);
bm=kD(9);
bmi=kD(10);
am=kD(11);
asm=kD(12);
bw=kD(13);
aw=kD(14);
asm2=kD(15);
aws=kD(16);
nw=kD(17);
Tw=kD(18);
as=kD(19);
bs=kD(20);
bs2=kD(21);
tau1=kD(22);
tau2=kD(23);
tau3=kD(24);
tau4=kD(25);
tau5=kD(26);
tau6=kD(27);
tau7=kD(28);
tau8=kD(29);
tau9=kD(30);
tau10=kD(31);
tau11=kD(32);
tau12=kD(33);
tau13=kD(34);
tau14=kD(35);
tau15=kD(36);
tau16=kD(37);
tau17=kD(38);
tau18=kD(39);
tau19=kD(40);
tau20=kD(41);
kD(42)=kD(42);
kD(43)=kD(43);
kD(44)=kD(44);
kD(45)=kD(45);
fixdsb1 =kD(46);
fixmsh =kD(47);
fixmgmt =kD(48);
basalp53act=kD(49);
kDDbasal=kD(50);
kDDE=kD(51);
kDEtop=kD(52);
Etop=kD(53);
kDnSP=kD(54,1);
kDkmSP=kD(55,1);
kDnSS=kD(56,1); 
kDnDS=kD(57,1); 
kDkmSS=kD(58,1); 
kDkmDS=kD(59,1); 


%%% Cell Cycle Rate Constants
aa=kC(1);
ab=kC(2);
ae=kC(3);
cdk1tot=kC(4); 
cdk2tot=kC(5); 
cdk4tot=kC(6);
Chk1tot=kC(7); 
eps=kC(8);
ib=kC(9);
ib1=kC(10);
ib2=kC(11);
ib3=kC(12);
K1=kC(13);
K1a=kC(14);
K1b=kC(15);
K1cdh1=kC(16);
K1chk=kC(17);
K1d=kC(18);
K1e=kC(19);
K1e2f=kC(20);
K1p27=kC(21);
K2=kC(22);
K2a=kC(23);
K2b=kC(24);
K2cdh1=kC(25);
K2chk=kC(26);
K2d=kC(27);
K2e=kC(28);
K2e2f=kC(29);
K2p27=kC(30);
K3=kC(31);
K3b=kC(32);
K4=kC(33);
K4b=kC(34);
K5a=kC(35);
K5b=kC(36);
K5e=kC(37);
K6a=kC(38);
K6b=kC(39);
K6e=kC(40);
K7b=kC(41);
K8b=kC(42);
Kacdc20=kC(43);
kc1=kC(44);
kc10=kC(45);
kc11=kC(46);
kc12=kC(47);
kc13=kC(48);
kc14=kC(49);
kc15=kC(50);
kc16=kC(51);
kc2=kC(52);
kc3=kC(53);
kc4=kC(54);
kc5=kC(55);
kc6=kC(56);
kc7=kC(57);
kc8=kC(58);
kc9=kC(59);
kca=kC(60);
kcd2=kC(61);
Kcdh1=kC(62);
kce=kC(63);
kcom1=kC(64);
kcom2=kC(65);
kcom3=kC(66);
kcom4=kC(67);
Kda=kC(68);
Kdb=kC(69);
Kdbcdc20=kC(70);
Kdbcdh1=kC(71);
kdcdc20a=kC(72);
kdcdc20i=kC(73);
kdcdh1a=kC(74);
kdcdh1i=kC(75);
Kdceskp2=kC(76);
Kdd=kC(77);
kdda=kC(78);
kddb=kC(79);
kddd=kC(80);
kdde=kC(81);
kddp21=kC(82);
kddp27=kC(83);
kddp27p=kC(84);
kddskp2=kC(85);
Kde=kC(86);
kde2f=kC(87);
kde2fp=kC(88);
kdecom1=kC(89);
kdecom2=kC(90);
kdecom3=kC(91);
kdecom4=kC(92);
Kdp27p=kC(93);
Kdp27skp2=kC(94);
kdpa=kC(95);
kdpai=kC(96);
kdpb=kC(97);
kdpbi=kC(98);
kdpe=kC(99);
kdpei=kC(100);
kdprb=kC(101);
kdprbp=kC(102);
kdprbpp=kC(103);
Kdskp2=kC(104);
kdWee1=kC(105);
kdWee1p=kC(106);
Ki10=kC(107);
Ki11=kC(108);
Ki12=kC(109);
Ki13=kC(110);
Ki14=kC(111);
Ki7=kC(112);
Ki8=kC(113);
Ki9=kC(114);
kpc1=kC(115);
kpc2=kC(116);
kpc3=kC(117);
kpc4=kC(118);
V1=kC(119);
V1cdh1=kC(120);
V1chk=kC(121);
V1e2f=kC(122);
V1p27=kC(123);
V2=kC(124);
V2cdh1=kC(125);
V2chk=kC(126);
V2e2f=kC(127);
V2p27=kC(128);
V3=kC(129);
V4=kC(130);
V6a=kC(131);
V6b=kC(132);
V6e=kC(133);
vcb=kC(134);
Vda=kC(135);
Vdb=kC(136);
Vdd=kC(137);
Vde=kC(138);
Vdp27p=kC(139);
Vdskp2=kC(140);
Vm1a=kC(141);
Vm1b=kC(142);
Vm1d=kC(143);
Vm1e=kC(144);
Vm2a=kC(145);
Vm2b=kC(146);
Vm2d=kC(147);
Vm2e=kC(148);
Vm3b=kC(149);
Vm4b=kC(150);
Vm5a=kC(151);
Vm5b=kC(152);
Vm5e=kC(153);
Vm7b=kC(154);
Vm8b=kC(155);
vs1p27=kC(156);
vs2p27=kC(157);
vscdc20i=kC(158);
vscdh1a=kC(159);
vse2f=kC(160);
vspai=kC(161);
vspbi=kC(162);
vspei=kC(163);
vsprb=kC(164);
vsskp2=kC(165);
vswee1=kC(166);
xa1=kC(167);
xa2=kC(168);
xb1=kC(169);
xb2=kC(170);
xe1=kC(171);
xe2=kC(172);
kcd1=kC(173); 

%%%% RATE LAWS

%% EIF4E sequestration by total mRNA in cell

mT=xE(1);
EIF4E_mT=xE(2);

ksynth_mT=kE(1);
kdeg_mT=kE(2);
kdeg_EIF4E_mT=kE(3);

vE(1)=kT1*EIF4E*mT;
vE(2)=kT2*EIF4E_mT;
vE(3)=ksynth_mT;
vE(4)=kdeg_mT*mT;
vE(5)=kdeg_EIF4E_mT*EIF4E_mT;


%% vRIBOSOME **

f1=((pS6K^nR)/(k50R^nR+pS6K^nR));
vbR=kbR0+(kbRi*f1); 
vdR=kdR0*xRibosome;

%% vTL **
rhs=(1/((EIF4Efree/(k50E+EIF4Efree))))*(Vn/Vc);

% % DNA Damage Mods
kTL(1)=bp/mExp_nM(1)*rhs;
kTL(2)=(bmi+bm*Mdm2pro)/mExp_nM(2)*rhs;
kTL(3)=(bw*Wip1pro)/mExp_nM(3)*rhs;
% Cell Cycle Mods
kTL(6)=vsprb/mExp_nM(6)*rhs*eps; 
kTL(7:9)=vse2f/sum(mExp_nM(7:9))*rhs*eps; 
kTL(10:12)=(kcd1+kcd2*E2F*(Ki7/(Ki7+pRB))*(Ki8/(Ki8+pRBp)))/sum(mExp_nM(10:12))*rhs*eps; 
kTL(13:14)=kce*E2F*(Ki9/(Ki9+pRB))*(Ki10/(Ki10+pRBp))/sum(mExp_nM(13:14))*rhs*eps; 
kTL(15)=vsskp2/mExp_nM(15)*rhs*eps; 
kTL(16)=vspei/mExp_nM(16)*rhs*eps;  
kTL(17)=vspai/mExp_nM(17)*rhs*eps;
kTL(18)=vspbi/mExp_nM(18)*rhs*eps;  
kTL(19)=kca*E2F*(Ki11/(Ki11+pRB))*(Ki12/(Ki12+pRBp))/mExp_nM(19)*rhs*eps;
kTL(20)=(vs1p27+(vs2p27*E2F*(Ki13/(Ki13+pRB))*(Ki14/(Ki14+pRBp))))/mExp_nM(20)*rhs*eps; 
kTL(21)=vscdh1a/mExp_nM(21)*rhs*eps; 
kTL(22)=vcb/mExp_nM(22)*rhs*eps;
kTL(23)=vscdc20i/mExp_nM(23)*rhs*eps; 
kTL(24)=vswee1/mExp_nM(24)*rhs*eps; 

kTLcdk1tot=(cdk1tot/mExp_nM(27))*kTLd(27); 
cdk1tot=(kTLcdk1tot*mExp_nM(27))/kTLd(27);

kTLcdk2tot=(cdk2tot/mExp_nM(28))*kTLd(28); 
cdk2tot=(kTLcdk2tot*mExp_nM(28))/kTLd(28);

avgktl=mean(kTLd(29:30));
summExp=sum(mExp_nM(29:30));
kTLcdk4tot=(cdk4tot/summExp)*avgktl; 
cdk4tot=(kTLcdk4tot*summExp)/avgktl;

kTLChk1tot=(kTLd(25)*Chk1tot)/mExp_nM(25); 
Chk1tot=(kTLChk1tot*mExp_nM(25))/kTLd(25);




if flagE
    vTL=kTL.*mMod*(EIF4E/(k50E+EIF4E)); 
else
    vTL=kTL.*mMod*(EIF4Efree/(k50E+EIF4Efree));
end

% Cell cycle proteins should not be affected by EIF4E, nor SGE.
vTL([6:9,13:30])=kTL([6:9,13:30]).*mMod([6:9,13:30])*(EIF4Efree/(k50E+EIF4Efree));

%% vTLCd Protein Conglomerates Degradation Reactions
vTLCd(1)=kTLCd(1)*p53inac;
vTLCd(2)=kTLCd(2)*Mdm2 ;
vTLCd(3)=kTLCd(3)*Wip1 ;
vTLCd(4)=kTLCd(4)*BRCA2;
vTLCd(5)=kTLCd(5)*MSH6;
vTLCd(6)=kTLCd(6)*MGMT;
vTLCd(7)=kTLCd(7)*ARF;
vTLCd(8)=kTLCd(8)*MDM4;
vTLCd(9)=kTLCd(9)*ATMinac;
vTLCd(10)=kTLCd(10)*ATRinac;
vTLCd(11)=0;%kTLCd(11)*pRB;
vTLCd(12)=0;%kTLCd(12)*E2F;
vTLCd(13)=0;%kTLCd(13)*Cd;
vTLCd(14)=0;%kTLCd(14)*Ce;
vTLCd(15)=0;%kTLCd(15)*Skp2;
vTLCd(16)=0;%kTLCd(16)*Pai;
vTLCd(17)=0;%kTLCd(17)*Pei;
vTLCd(18)=0;%kTLCd(18)*Pbi;
vTLCd(19)=0;%kTLCd(19)*Ca;
vTLCd(20)=0;%kTLCd(20)*p27;
vTLCd(21)=0;%kTLCd(21)*Cdh1a;
vTLCd(22)=0;%kTLCd(22)*Cb;
vTLCd(23)=0;%kTLCd(23)*Cdc20i;
vTLCd(24)=0;%kTLCd(24)*Wee1;
vTLCd(25)=0;%kTLCd(25)*Chk1;
vTLCd(26)=0;%kTLCd(26)*p21;
vTLCd(27)=kTLCd(27)*L ;
vTLCd(28)=kTLCd(28)*R ;
vTLCd(29)=kTLCd(29)*flip ;
vTLCd(30)=kTLCd(30)*pC8 ;
vTLCd(31)=kTLCd(31)*Bar ;
vTLCd(32)=kTLCd(32)*pC3 ;
vTLCd(33)=kTLCd(33)*pC6 ;
vTLCd(34)=kTLCd(34)*XIAP ;
vTLCd(35)=kTLCd(35)*PARP ;
vTLCd(36)=kTLCd(36)*Bid ;
vTLCd(37)=kTLCd(37)*Bcl2c ;
vTLCd(38)=kTLCd(38)*Bax ;
vTLCd(39)=kTLCd(39)*CytoCm ;
vTLCd(40)=kTLCd(40)*Smacm ;
vTLCd(41)=kTLCd(41)*Apaf ;
vTLCd(42)=kTLCd(42)*pC9 ;
vTLCd(43)=kTLCd(43)*BAD;
vTLCd(44)=kTLCd(44)*PUMA;
vTLCd(45)=kTLCd(45)*NOXA;
vTLCd(46)=kTLCd(46)*BIM;
vTLCd(47)=kTLCd(47)*E;
vTLCd(48)=kTLCd(48)*H;
vTLCd(49)=kTLCd(49)*HGF;
vTLCd(50)=kTLCd(50)*P;
vTLCd(51)=kTLCd(51)*F;
vTLCd(52)=kTLCd(52)*I;
vTLCd(53)=kTLCd(53)*IN;
vTLCd(54)=kTLCd(54)*E1;
vTLCd(55)=kTLCd(55)*E2;
vTLCd(56)=kTLCd(56)*E3;
vTLCd(57)=kTLCd(57)*E4;
vTLCd(58)=kTLCd(58)*Ev3;
vTLCd(59)=kTLCd(59)*Met;
vTLCd(60)=kTLCd(60)*Pr;
vTLCd(61)=kTLCd(61)*Fr;
vTLCd(62)=kTLCd(62)*Ir;
vTLCd(63)=kTLCd(63)*Isr;
vTLCd(64)=kTLCd(64)*IRS;
vTLCd(65)=kTLCd(65)*Sp;
vTLCd(66)=kTLCd(66)*Cbl;
vTLCd(67)=kTLCd(67)*G2;
vTLCd(68)=kTLCd(68)*PLCg;
vTLCd(69)=kTLCd(69)*PI3KC1;
vTLCd(70)=kTLCd(70)*PI3KR1;
vTLCd(71)=kTLCd(71)*PI3K2;
vTLCd(72)=kTLCd(72)*GRP;
vTLCd(73)=kTLCd(73)*RasD;
vTLCd(74)=kTLCd(74)*NF1;
vTLCd(75)=kTLCd(75)*CRaf;
vTLCd(76)=kTLCd(76)*BRaf;
vTLCd(77)=kTLCd(77)*MEK;
vTLCd(78)=kTLCd(78)*MKP3;
vTLCd(79)=kTLCd(79)*RSK;
vTLCd(80)=kTLCd(80)*MKP1;
vTLCd(81)=kTLCd(81)*cFos;
vTLCd(82)=kTLCd(82)*cJun;
vTLCd(83)=kTLCd(83)*cMyc;
vTLCd(84)=kTLCd(84)*bCATENIN;
vTLCd(85)=kTLCd(85)*PTEN;
vTLCd(86)=kTLCd(86)*AKT;
vTLCd(87)=kTLCd(87)*PDK1;
vTLCd(88)=kTLCd(88)*Rictor;
vTLCd(89)=kTLCd(89)*mTOR;
vTLCd(90)=kTLCd(90)*GSK3b;
vTLCd(91)=kTLCd(91)*TSC1;
vTLCd(92)=kTLCd(92)*TSC2;
vTLCd(93)=kTLCd(93)*PKC;
vTLCd(94)=kTLCd(94)*RKIP;
vTLCd(95)=kTLCd(95)*ERK;
vTLCd(96)=kTLCd(96)*FOXO;
vTLCd(97)=kTLCd(97)*RhebD;
vTLCd(98)=kTLCd(98)*Raptor;
vTLCd(99)=kTLCd(99)*S6K;
vTLCd(100)=kTLCd(100)*EIF4EBP1;
vTLCd(101)=kTLCd(101)*SOS;
vTLCd(102)=kTLCd(102)*EIF4E;

%% vD **
%vD(1)= bp;
vD(2)= ampi*Mdm2*p53inac; 
vD(3)= basalp53act + bsp*p53inac*(ATMP.^ns./(ATMP.^ns+Ts.^ns)+ATRac.^ns./(ATRac.^ns+Ts.^ns)); 
vD(4)= awpa*Wip1*p53ac;
vD(5)= api*p53inac; 
vD(6)= ampa*Mdm2*p53ac; 
%vD(7)= bm*Mdm2pro;
%vD(8)= bmi;
vD(9)= asm*Mdm2*ATMP; 
vD(10)= asm2*Mdm2*ATRac; 
vD(11)= am*Mdm2;
%vD(12)= bw*Wip1pro;
vD(13)= aw*Wip1; 
vD(14)= bs*((damageDSB^kDnDS)/((kDkmDS^kDnDS)+(damageDSB^kDnDS)));
vD(15)= aws*ATMP*Wip1.^nw./(Wip1.^nw+Tw.^nw); 
vD(16)= as*ATMP; 
vD(17)= bs2*((damageSSB^kDnSS)/((kDkmSS^kDnSS)+(damageSSB^kDnSS)));
vD(18)= as*ATRac; 
vD(19)= Mdm2product1/tau1;
vD(20)= p53ac/tau1;
vD(21)= Mdm2product2/tau2;
vD(22)= Mdm2product1/tau2;
vD(23)= Mdm2product3/tau3;
vD(24)= Mdm2product2/tau3;
vD(25)= Mdm2product4/tau4;
vD(26)= Mdm2product3/tau4;
vD(27)= Mdm2product5/tau5;
vD(28)= Mdm2product4/tau5;
vD(29)= Mdm2product6/tau6;
vD(30)= Mdm2product5/tau6;
vD(31)= Mdm2product7/tau7;
vD(32)= Mdm2product6/tau7;
vD(33)= Mdm2product8/tau8;
vD(34)= Mdm2product7/tau8;
vD(35)= Mdm2product9/tau9;
vD(36)= Mdm2product8/tau9;
vD(37)= Mdm2pro/tau10;
vD(38)= Mdm2product9/tau10;
vD(39)= Wip1product1/tau11;
vD(40)= p53ac/tau11;
vD(41)= Wip1product2/tau12;
vD(42)= Wip1product1/tau12;
vD(43)= Wip1product3/tau13;
vD(44)= Wip1product2/tau13;
vD(45)= Wip1product4/tau14;
vD(46)= Wip1product3/tau14;
vD(47)= Wip1product5/tau15;
vD(48)= Wip1product4/tau15;
vD(49)= Wip1product6/tau16;
vD(50)= Wip1product5/tau16;
vD(51)= Wip1product7/tau17;
vD(52)= Wip1product6/tau17;
vD(53)= Wip1product8/tau18;
vD(54)= Wip1product7/tau18;
vD(55)= Wip1product9/tau19;
vD(56)= Wip1product8/tau19;
vD(57)= Wip1pro/tau20;
vD(58)= Wip1product9/tau20;
vD(59)=kD(42)*ARF*Mdm2;
vD(60)=kD(43)*ARF*pMdm2;
vD(61)=kD(44)*MDM4*p53ac;
vD(62)=kD(45)*p53ac_MDM4;
vD(63)=fixdsb1*BRCA2*damageDSB;
vD(64)=fixmsh*MSH6*damageSSB;
vD(65)=fixmgmt*MGMT*damageSSB;
vD(66)=(kDDbasal+kDDE*(Etop/(Etop+kDEtop)))*(((Ma+Me)^kDnSP)/(((Ma+Me)^kDnSP)+(kDkmSP^kDnSP))); 


%% vC **
%vC(1)=vsprb; %pRb(Rb) synthesis 
vC(2)=kpc1*pRB*E2F; %pRb binding to E2F  
vC(3)=kpc2*pRBc1; %dissociation of that pRbE2F complex
vC(4)=V1*(pRB/(K1+pRB))*(Md+Mdp27); %pRb phosphorylation
vC(5)=V2*(pRBp/(K2+pRBp)); %pRb dephosphorylation
vC(6)=kdprb*pRB; %pRb degradation
vC(7)=V3*(pRBp/(K3+pRBp))*Me; %pRbp phosphorylating to pRbpp
vC(8)=V4*(pRBpp/(K4+pRBpp)); %pRbpp dephosphorylation to pRbp
vC(9)=kpc3*pRBp*E2F; %pRbp binding to E2F
vC(10)=kpc4*pRBc2; %pRbp dissociation from E2F
vC(11)=kdprbp*pRBp; %pRbp degradation
vC(12)=kdprbpp*pRBpp; %pRbpp degradation
%vC(13)=vse2f; %E2F synthesis
vC(14)=V1e2f*Ma*(E2F/(E2F+K1e2f)); %E2F phosphorylation by active CYCACDK2
vC(15)=V2e2f*(E2Fp/(E2Fp+K2e2f)); %pE2F dephosphorylation to E2F
vC(16)=kde2f*E2F; %E2F degradation
vC(17)=kde2fp*E2Fp; %pE2F degradation
%vC(18)=kcd1+kcd2*E2F*(Ki7/(Ki7+pRB))*(Ki8/(Ki8+pRBp)); %cycd synthesis due to E2f inhibited by pRb and pRbp
vC(19)=kcom1*Cd*(cdk4tot-(Mdi+Md+Mdp27+Mdp21)); %cycd binding to cdk4/6 
vC(20)=kdecom1*Mdi; %dissociation of CYCD from CYCDCDK4/6(i) complex
vC(21)=Vdd*(Cd/(Kdd+Cd)); %maximum cycd degradation 
vC(22)=kddd*Cd; %non specific cycd degradation
vC(23)=Vm2d*(Md/(K2d+Md)); %inactivation of CYCDCK4/6(a) or Md
vC(24)=Vm1d*(Mdi/(K1d+Mdi)); %activation of CYCDCDK4/6(i) or Mdi
vC(25)=kc1*Md*p27; %CYCDCDK4/6(a) or Md binding to p27
vC(26)=kc2*Mdp27; %dissociation of CYCDCDK4/6p27 coplex or Mdp27 to Md
%vC(27)=kce*E2F*(Ki9/(Ki9+pRB))*(Ki10/(Ki10+pRBp)); %CYCE synthesis %%Keep like this for now.
vC(28)=kcom2*Ce*(cdk2tot-(Mei+Me+Mep27+Mep21+Mai+Ma+Map27+Map21)); %CYCE going into complex with cdk2
vC(29)=kdecom2*Mei; %dissociation of complex CYCECDK2(i)
vC(30)=Vde*(Skp2/(Kdceskp2+Skp2))*(Ce/(Kde+Ce)); %CYCE degradation
vC(31)=Vm2e*(Wee1+ib2)*(Me/(K2e+Me)); %inactivation of CYCECDK2
vC(32)=Vm1e*Pe*(Mei/(K1e+Mei)); %activation of CYCECDK2
vC(33)=kc3*Me*p27; %CYCECDK2 binding to p27
vC(34)=kc4*Mep27; %dissociation of CYCECDK2p27
%vC(35)=vsskp2; %synthesis of SKP2
vC(36)=Vdskp2*(Skp2/(Kdskp2+Skp2))*(Cdh1a/(Kcdh1+Cdh1a)); %max SKP2 degradation
vC(37)=kddskp2*Skp2; %non-specific SKP2 degradation
%vC(38)=vspei; %synthesis of inactivated CDC25(i) or pei
vC(39)=V6e*(xe1+xe2*Chk1)*(Pe/(K6e+Pe)); %Pei inactivation
vC(40)=Vm5e*(Me+ae)*(Pei/(K5e+Pei)); %Pei activation
vC(41)=kdpei*Pei; %Pei degradation
vC(42)=kdpe*Pe; %Pe degradation
%vC(43)=kca*E2F*(Ki11/(Ki11+pRB))*(Ki12/(Ki12+pRBp)); %CYCA synthesis induced by E2F inhibited by pRb and pRbp %%Keep like this for now.
vC(44)=kcom3*Ca*(cdk2tot-(Mei+Me+Mep27+Mep21+Mai+Ma+Map27+Map21)); %CYCA binding to CDK2 making CYCACDK2(i) Mai
vC(45)=kdecom3*Mai; %dissociation of Mai complex
vC(46)=Vda*(Ca/(Kda+Ca))*(Cdc20a/(Kacdc20+Cdc20a)); %CYCA degradation induced by Cdc20(a)
vC(47)=kdda*Ca; %non specific CYCA degradation
vC(48)=Vm2a*(Wee1+ib3)*(Ma/(K2a+Ma)); %inactivation of CYCACDK2(a) or Ma complex
vC(49)=Vm1a*Pa*(Mai/(K1a+Mai)); %activation of CYCACDK2(i) or Mai complex
vC(50)=kc5*Ma*p27; %Ma binding to P27 making CYCACDK2p27
vC(51)=kc6*Map27; %dissociation of Map27 to Map and p27
vC(52)=V1p27*Me*(p27/(p27+K1p27)); %inactivation of p27 by CYCECDK2(a)
vC(53)=V2p27*(p27p/(p27p+K2p27)); %activation of p27
vC(54)=Vdp27p*(Skp2/(Skp2+Kdp27skp2))*(p27p/(p27p+Kdp27p)); %max p27p degradation
vC(55)=kddp27p*p27p; %non specific p27p degradation
vC(56)=V2cdh1*(Cdh1a/(K2cdh1+Cdh1a))*(Ma+Mb); %inactivation of Cdh1(a)
vC(57)=V1cdh1*(Cdh1i/(K1cdh1+Cdh1i)); %activation of Cdh1(1)
vC(58)=kdcdh1i*Cdh1i; %inactive Cdh1(i) degradation
%vC(59)=vscdh1a; %synthesis of Cdh1(a)
vC(60)=kdcdh1a*Cdh1a; %Cdh1(a) degradation
%vC(61)=vspai; %CDC25(i) Pai synthesis
vC(62)=V6a*(xa1+xa2*Chk1)*(Pa/(K6a+Pa)); %Pa Cdc25(a) inactivation to Pai
vC(63)=Vm5a*(Ma+aa)*(Pai/(K5a+Pai)); %Pai activation to Pa
vC(64)=kdpai*Pai; %Pai degradation
vC(65)=kdpa*Pa; %Pa degradation
%vC(66)=vcb; %CYCB synthesis
vC(67)=kcom4*Cb*(cdk1tot-(Mbi+Mb+Mbp27+Mbp21)); %CYCB forming CYCBCDK1(i) complex
vC(68)=kdecom4*Mbi; %dissociation of CYCBCDK1(i) complex
vC(69)=Vdb*(Cb/(Kdb+Cb))*((Cdc20a/(Kdbcdc20+Cdc20a))+(Cdh1a/(Kdbcdh1+Cdh1a))); %max CYCB degradation induced by Cdc20(a)
vC(70)=kddb*Cb; %non specific CYCB degradation
vC(71)=Vm2b*(Wee1+ib1)*(Mb/(K2b+Mb)); %Mb inactivation to Mbi
vC(72)=Vm1b*Pb*(Mbi/(K1b+Mbi)); %Mbi becoming active
vC(73)=kc7*Mb*p27; %Mb (active CYCBCDK1) binding to p27
vC(74)=kc8*Mbp27; %dissociation of Mb and p27
%vC(75)=vscdc20i; %cdc20(i) synthesis
vC(76)=Vm3b*Mb*(Cdc20i/(K3b+Cdc20i)); %activation of cdc20(i) through phosphorylation by CYCBCDK1
vC(77)=Vm4b*(Cdc20a/(K4b+Cdc20a)); %cdc20(a) inactivation to cdc20(i)
vC(78)=kdcdc20i*Cdc20i; %cdc20i degradation
vC(79)=kdcdc20a*Cdc20a; %degradation of cdc20a
%vC(80)=vspbi; %synthesis of inactive CDC25 Pbi
vC(81)=V6b*(xb1+xb2*Chk1)*(Pb/(K6b+Pb)); %inactivation of CDC25 (Pb)
vC(82)=Vm5b*(Mb+ab)*(Pbi/(K5b+Pbi)); %activation of CDC25(i) (Pbi)
vC(83)=kdpbi*Pbi; %degradation of inactive CDC25 (Pbi)
vC(84)=kdpb*Pb; %degradation of active CDC25 (Pb)
%vC(85)=vswee1; %Wee1 synthesis
vC(86)=Vm7b*(Mb+ib)*(Wee1/(K7b+Wee1)); %wee1 inactivation through phosphorylation due to CYCBCDK1
vC(87)=Vm8b*(Wee1p/(K8b+Wee1p)); %wee1p(wee1(i)) activation through dephosphorylation
vC(88)=kdWee1*Wee1; %wee1 degradation
vC(89)=kdWee1p*Wee1p; %wee1p degradation
vC(90)=V1chk*ATRac*((Chk1tot-Chk1)/(K1chk+(Chk1tot-Chk1))); %activation of Chk1 through phosphorylation by kinase ATR 
vC(91)=V2chk*(Chk1/(K2chk+Chk1)); %Chk1 inactivation through dephosphorylation
vC(92)=kdde*Ce;%non-specific cycE degradation
%vC(93)=vs1p27; %+ basal p27 synthesis
%vC(94)=vs2p27*E2F*(Ki13/(Ki13+pRB))*(Ki14/(Ki14+pRBp)); %+ synthesis of p27 induced by E2F, inhibited by pRB and pRBp %%Keep like this for now.
vC(95)=kddp27*p27; %- non-specific p27 degradation.
vC(96)=kc9*Md*p21; %CYCDCDK4/6(a) or Md binding to p21
vC(97)=kc10*Mdp21; %dissociation of CYCDCDK4/6p21 coplex or Mdp21 to Md
vC(98)=kc11*Me*p21; %CYCECDK2 binding to p21
vC(99)=kc12*Mep21; %dissociation of CYCECDK2p21
vC(100)=kc13*Ma*p21; %Ma binding to P21 making CYCACDK2p21
vC(101)=kc14*Map21; %dissociation of Map21 to Map and p21
vC(102)=kc15*Mb*p21; %Mb (active CYCBCDK1) binding to p21
vC(103)=kc16*Mbp21; %dissociation of Mb and p21
vC(104)=kddp21*p21; %-non-specific p21 degradation
vC(1:104)=vC(1:104)*eps;


%% vA **
vA(1)=kA(1)*L*R *Vc/Ve;
vA(2)=kA(2)*L_R;
vA(3)=kA(3)*L_R;
vA(4)=kA(4)*Ractive*flip;
vA(5)=kA(5)*Ractive_flip;
vA(6)=kA(6)*Ractive*pC8;
vA(7)=kA(7)*Ractive_pC8;
vA(8)=kA(8)*Ractive_pC8;
vA(9)=kA(9)*C8*Bar;
vA(10)=kA(10)*C8_Bar;
vA(11)=kA(11)*C8*pC3;
vA(12)=kA(12)*C8_pC3;
vA(13)=kA(13)*C8_pC3;
vA(14)=kA(14)*C3*pC6;
vA(15)=kA(15)*C3_pC6;
vA(16)=kA(16)*C3_pC6;
vA(17)=kA(17)*pC8*C6;
vA(18)=kA(18)*C6_C8;
vA(19)=kA(19)*C6_C8;
vA(20)=kA(20)*C3*XIAP;
vA(21)=kA(21)*C3_XIAP;
vA(22)=kA(22)*C3_XIAP;
vA(23)=kA(23)*C3*PARP;
vA(24)=kA(24)*C3_PARP;
vA(25)=kA(25)*C3_PARP;
vA(26)=kA(26)*C8*Bid;
vA(27)=kA(27)*C8_Bid;
vA(28)=kA(28)*C8_Bid;
vA(29)=kA(29)*tBid*Bcl2c;
vA(30)=kA(30)*tBid_Bcl2c;
vA(31)=kA(31)*tBid*Bax; 
vA(32)=kA(32)*tBid_Bax; 
vA(33)=kA(33)*tBid_Bax; 
vA(34)=kA(34)*Baxactive;
vA(35)=kA(35)*Baxm;
vA(36)=kA(36)*Baxm*Bcl2;
vA(37)=kA(37)*Baxm_Bcl2;
vA(38)=kA(38)*Baxm^2;
vA(39)=kA(39)*Bax2;
vA(40)=kA(40)*Bcl2*Bax2;
vA(41)=kA(41)*Bax2_Bcl2;
vA(42)=kA(42)*Bax2^2;
vA(43)=kA(43)*Bax4;
vA(44)=kA(44)*Bcl2*Bax4;
vA(45)=kA(45)*Bax4_Bcl2;
vA(46)=kA(46)*M*Bax4;
vA(47)=kA(47)*Bax4_M;
vA(48)=kA(48)*Bax4_M;
vA(49)=kA(49)*Mactive*CytoCm;
vA(50)=kA(50)*Mactive_CytoCm;
vA(51)=kA(51)*Mactive_CytoCm;
vA(52)=kA(52)*Mactive*Smacm;
vA(53)=kA(53)*Mactive_Smacm;
vA(54)=kA(54)*Mactive_Smacm;
vA(55)=kA(55)*CytoCr;
vA(56)=kA(56)*CytoC;
vA(57)=kA(57)*CytoC*Apaf;
vA(58)=kA(58)*CytoC_Apaf;
vA(59)=kA(59)*CytoC_Apaf;
vA(60)=kA(60)*Apafactive*pC9;
vA(61)=kA(61)*Apop;
vA(62)=kA(62)*pC3*Apop;
vA(63)=kA(63)*Apop_C3;
vA(64)=kA(64)*Apop_C3;
vA(65)=kA(65)*Smacr;
vA(66)=kA(66)*Smac;
vA(67)=kA(67)*XIAP*Apop;
vA(68)=kA(68)*Apop_XIAP;
vA(69)=kA(69)*XIAP*Smac;
vA(70)=kA(70)*Smac_XIAP;
vA(71)=kA(71)*Bcl2c*BAD; 
vA(72)=kA(72)*Bcl2c_BAD;
vA(73)=kA(73)*Bcl2c*PUMA;
vA(74)=kA(74)*Bcl2c_PUMA;
vA(75)=kA(75)*Bcl2c*NOXA;
vA(76)=kA(76)*Bcl2c_NOXA;
vA(77)=kA(77)*BIM*Bax; 
vA(78)=kA(78)*BIM_Bax; 
vA(79)=kA(79)*BIM_Bax; 
vA(80)=kA(80)*Bcl2c*BIM; 
vA(81)=kA(81)*Bcl2c_BIM;
vA(82)=kA(82)*Mactive; 
vA(83)=kA(83)*Smacr;
vA(84)=kA(84)*CytoCr;
vA(85)=kA(85)*Bcl2c; 
vA(86)=kA(86)*Bcl2; 
vA(87)=kA(87)*pC8; 

%% vR **
vR(1)=kR(1)*E1*E *Vc/Ve;
vR(2)=kR(2)*EE1 ;
vR(3)=kR(3)*EE1*E2 ;
vR(4)=kR(4)*EE1E2 ;
vR(5)=kR(5)*EE1*E1 ;
vR(6)=kR(6)*EE1E1 ;
vR(7)=kR(7)*EE1E1*E *Vc/Ve;
vR(8)=2*kR(8)*EE1EE1 ;
vR(9)=kR(9)*EE1EE1 ;
vR(10)=kR(10)*EE1*EE1 ;
vR(11)=kR(11)*EE1*E3 ; 
vR(12)=kR(12)*EE1E3 ;
vR(13)=kR(13)*EE1E3*H *Vc/Ve;
vR(14)=kR(14)*EE1HE3 ;
vR(15)=kR(15)*EE1HE3 ;
vR(16)=kR(16)*EE1*HE3 ;
vR(17)=kR(17)*EE1*Ev3 ;
vR(18)=kR(18)*EE1Ev3 ;
vR(19)=kR(19)*EE1*E4 ; 
vR(20)=kR(20)*EE1E4 ;
vR(21)=kR(21)*EE1E4*H *Vc/Ve;
vR(22)=kR(22)*EE1HE4 ;
vR(23)=kR(23)*EE1HE4 ;
vR(24)=kR(24)*EE1*HE4 ;
vR(25)=kR(25)*H*E3 *Vc/Ve;
vR(26)=kR(26)*HE3 ;
vR(27)=kR(27)*HE3*E2 ;
vR(28)=kR(28)*E2HE3 ;
vR(29)=kR(29)*HE3*E1 ;
vR(30)=kR(30)*E1HE3 ;
vR(31)=kR(31)*E1HE3*E *Vc/Ve;
vR(32)=kR(32)*EE1HE3 ;
vR(33)=kR(33)*HE3*E3 ;
vR(34)=kR(34)*HE3E3 ;
vR(35)=kR(35)*HE3E3*H *Vc/Ve;
vR(36)=2*kR(36)*HE3HE3 ;
vR(37)=kR(37)*HE3HE3 ;
vR(38)=kR(38)*HE3*HE3 ;
vR(39)=kR(39)*HE3*Ev3 ;
vR(40)=kR(40)*HE3Ev3 ;
vR(41)=kR(41)*HE3*E4 ;
vR(42)=kR(42)*HE3E4 ;
vR(43)=kR(43)*HE3E4*H *Vc/Ve; 
vR(44)=kR(44)*HE3HE4 ;
vR(45)=kR(45)*HE3HE4 ;
vR(46)=kR(46)*HE3*HE4 ;
vR(47)=kR(47)*H*E4 *Vc/Ve;
vR(48)=kR(48)*HE4 ;
vR(49)=kR(49)*HE4*E2 ;
vR(50)=kR(50)*E2HE4 ;
vR(51)=kR(51)*HE4*E1 ;
vR(52)=kR(52)*E1HE4 ;
vR(53)=kR(53)*E1HE4*E *Vc/Ve;
vR(54)=kR(54)*EE1HE4 ;
vR(55)=kR(55)*HE4*E3 ;
vR(56)=kR(56)*E3HE4 ;
vR(57)=kR(57)*E3HE4*H *Vc/Ve;
vR(58)=kR(58)*HE3HE4 ;
vR(59)=kR(59)*HE4*Ev3 ;
vR(60)=kR(60)*HE4Ev3 ;
vR(61)=kR(61)*HE4*E4 ;
vR(62)=kR(62)*HE4E4 ;
vR(63)=kR(63)*HE4E4*H *Vc/Ve;
vR(64)=2*kR(64)*HE4HE4 ;
vR(65)=kR(65)*HE4HE4 ;
vR(66)=kR(66)*HE4*HE4 ;
vR(67)=kR(67)*E1*ppERK ; 
vR(68)=kR(68)*pE1 ;
vR(69)=kR(69)*E2*ppERK ; 
vR(70)=kR(70)*pE2 ;
vR(71)=kR(71)*E4*ppERK ; 
vR(72)=kR(72)*pE4 ;
vR(73)=kR(73)*Met*HGF *Vc/Ve;
vR(74)=kR(74)*HGF_Met ;
vR(75)=kR(75)*HGF_Met*Met ;
vR(76)=kR(76)*HGF_Met_Met ;
vR(77)=kR(77)*HGF_Met_Met*HGF *Vc/Ve;
vR(78)=2*kR(78)*HGF_Met_HGF_Met ;
vR(79)=kR(79)*HGF_Met_HGF_Met ;
vR(80)=kR(80)*HGF_Met*HGF_Met ;
vR(81)=kR(81)*Pr*P *Vc/Ve;
vR(82)=kR(82)*PPr ;
vR(83)=kR(83)*PPr*PPr ;
vR(84)=kR(84)*PPrPPr ;
vR(85)=2*kR(85)*PPrPPr ;
vR(86)=kR(86)*PPrPr ;
vR(87)=kR(87)*Fr*F *Vc/Ve;
vR(88)=kR(88)*FFr;
vR(89)=kR(89)*FFr*FFr;
vR(90)=kR(90)*FFrFFr;
vR(91)=2*kR(91)*FFrFFr;
vR(92)=kR(92)*F*FFrFr *Vc/Ve;
vR(93)=kR(93)*FFrFr;
vR(94)=kR(94)*Fr*FFr;
vR(95)=kR(95)*Ir*Ir;
vR(96)=kR(96)*IrIr;
vR(97)=kR(97)*Isr*Isr;
vR(98)=kR(98)*Isr_Isr;
vR(99)=kR(99)*E1_ppERK;
vR(100)=kR(100)*E1_ppERK;
vR(101)=kR(101)*E2_ppERK;
vR(102)=kR(102)*E2_ppERK;
vR(103)=kR(103)*E4_ppERK;
vR(104)=kR(104)*E4_ppERK;
vR(105)=kR(105)*E1*E1;
vR(106)=kR(106)*E1E1;
vR(107)=kR(107)*E1*E2;
vR(108)=kR(108)*E1E2;
vR(109)=kR(109)*E1*E3;
vR(110)=kR(110)*E1E3;
vR(111)=kR(111)*E1*E4;
vR(112)=kR(112)*E1E4;
vR(113)=kR(113)*E2*E2;
vR(114)=kR(114)*E2E2;
vR(115)=kR(115)*E2*E3;
vR(116)=kR(116)*E2E3;
vR(117)=kR(117)*E2*E4;
vR(118)=kR(118)*E2E4;
vR(119)=kR(119)*E3*E4;
vR(120)=kR(120)*E3E4;
vR(121)=kR(121)*E4*E4;
vR(122)=kR(122)*E4E4;
vR(123)=2*kR(123)*E*E1E1 *Vc/Ve;
vR(124)=kR(124)*EE1E1;
vR(125)=kR(125)*E*E1E2 *Vc/Ve;
vR(126)=kR(126)*EE1E2;
vR(127)=kR(127)*E*E1E3 *Vc/Ve;
vR(128)=kR(128)*EE1E3;
vR(129)=kR(129)*H*E1E3 *Vc/Ve;
vR(130)=kR(130)*E1HE3;
vR(131)=kR(131)*E*E1E4 *Vc/Ve;
vR(132)=kR(132)*EE1E4;
vR(133)=kR(133)*H*E1E4 *Vc/Ve;
vR(134)=kR(134)*E1HE4;
vR(135)=kR(135)*H*E2E3 *Vc/Ve;
vR(136)=kR(136)*E2HE3;
vR(137)=kR(137)*H*E2E4 *Vc/Ve;
vR(138)=kR(138)*E2HE4;
vR(139)=2*kR(139)*H*E3E4 *Vc/Ve;
vR(140)=kR(140)*E3HE4;
vR(141)=2*kR(141)*H*E4E4 *Vc/Ve;
vR(142)=kR(142)*HE4E4;
vR(143)=kR(143)*Met*Met;
vR(144)=kR(144)*Met_Met;
vR(145)=2*kR(145)*HGF*Met_Met;
vR(146)=kR(146)*HGF_Met_Met;
vR(147)=kR(147)*Fr*Fr;
vR(148)=kR(148)*FrFr;
vR(149)=2*kR(149)*F*FrFr;
vR(150)=kR(150)*FFrFr;
vR(151)=kR(151)*IrIr*I *Vc/Ve;
vR(152)=kR(152)*IIrIr;
vR(153)=kR(153)*IIrIr*I *Vc/Ve;
vR(154)=kR(154)*IIrIrI;  
vR(155)=kR(155)*Isr_Isr*IN *Vc/Ve;
vR(156)=kR(156)*IN_Isr_Isr; 
vR(157)=kR(157)*IN_Isr_Isr*IN *Vc/Ve; 
vR(158)=kR(158)*IN_Isr_Isr_IN; 


%% vRP **
fac=1; 

ps_G2_SOS=[9;10;10;10;8;8;11;11;7;8;8;9;9;10;11;11;9;12;12;2;2;1;1;8;8;5;7;5;7]*fac;
ps_PI3K1=[4;2;2;2;10;10;4;4;12;10;10;12;12;6;4;4;12;6;6;6;6;6;6;14;16;10.5;14.5;10.5;14.5]*fac;
ps_PI3K2=ps_PI3K1;
ps_PLCg=[4;4;4;4;3;3;3;3;3;3;3;2;2;3;3;3;2;2;2;4;4;5;5;2;6;1.5;3.5;1.5;3.5]*fac;

% Equations
vRP1 = kRP1.*SCD; % Phosphorylation
vRP2 = kRP2.*pSCD; % Dephosphorylation

vRP3 = kRP3.*SCD*Sp; % Binding Spouty
vRP4 = kRP4.*Sp_SCD; % Dissociation Sprouty

kmCblint=2;
vRP5 = kRP5.*(pSCD./(kmCblint+pSCD))*Cbl + kRP5(9).*pSCD; % Internalization

vRP6 = kRP6.*pSCDint; % Dephosphorylation
vRP7 = kRP7.*SCDint; % Phosphorylation

vRP8 = kRP8.*SCDint; % Recycling

vRP9 = ps_G2_SOS.*kRP9.*pSCDint_bind*G2_SOS; % Int binding G2-SOS
vRP10 = kRP10.*pSCDint_G2_SOS; % Int dissociating from G2-SOS

vRP11 = ps_G2_SOS.*kRP11.*pSCD_bind.*G2_SOS; %pSCD binding G2-SOS
vRP12 = kRP12.*pSCD_G2_SOS; %pSCD dissociating G2-SOS

vRP13 = ps_PLCg.*kRP13.*pSCD_bind.*PLCg; 
vRP14 = kRP14.*pSCD_PLCg;

vRP15 = ps_PI3K1.*kRP15.*pSCD_bind.*PI3K1;
vRP16 = kRP16.*pSCD_PI3K1;

vRP17 = ps_PI3K2.*kRP17.*pSCD_bind.*PI3K2;
vRP18 = kRP18.*pSCD_PI3K2;

vRP19 = kRP19.*pSCDint_G2_SOS.*RasD; %Binding to RasDint
vRP20 = kRP20.*pSCDint_G2_SOS_RasD; %Dissociation
vRP21 = kRP21.*pSCDint_G2_SOS_RasD; %Making product

vRP22 = kRP22.*pSCD_G2_SOS.*RasD; %Binding to RasD
vRP23 = kRP23.*pSCD_G2_SOS_RasD; %Dissociation
vRP24 = kRP24.*pSCD_G2_SOS_RasD; %Making product

vRP25 = kRP25.*pSCD_PLCg.*PIP2; %Binding to PIP2
vRP26 = kRP26.*pSCD_PLCg_PIP2; %Dissociation
vRP27 = kRP27.*pSCD_PLCg_PIP2; %Making product

vRP28 = kRP28.*pSCD_PI3K1.*PIP2; %Binding to PIP2
vRP29 = kRP29.*pSCD_PI3K1_PIP2; %Dissociation
vRP30 = kRP30.*pSCD_PI3K1_PIP2; %Making product

vRP31 = kRP31.*pSCD_PI3K2.*PIP; %Binding to PIP2
vRP32 = kRP32.*pSCD_PI3K2_PIP; %Dissociation
vRP33 = kRP33.*pSCD_PI3K2_PIP; %Making product

% IGFR and INSR binding IRS
ps_IRS=2*4;
vRP34(1)=ps_IRS*fac.*kRP34(1).*pIIrIr.*IRS; 
vRP34(2)=kRP34(2).*pIIrIr_IRS;
vRP34(3)=ps_IRS*fac.*kRP34(3).*pIN_Isr_Isr.*IRS;
vRP34(4)=kRP34(4).*pIN_Isr_Isr_IRS;
vRP34(5)=ps_IRS*fac.*kRP34(5).*pIIrIrI.*IRS; 
vRP34(6)=kRP34(6).*pIIrIrI_IRS;
vRP34(7)=ps_IRS*fac.*kRP34(7).*pIN_Isr_Isr_IN.*IRS;
vRP34(8)=kRP34(8).*pIN_Isr_Isr_IN_IRS;
vRP34(9)=ps_IRS*fac.*kRP34(9).*pIIrIr_int.*IRS; 
vRP34(10)=kRP34(10).*pIIrIr_int_IRS;
vRP34(11)=ps_IRS*fac.*kRP34(11).*pIN_Isr_Isr_int.*IRS;
vRP34(12)=kRP34(12).*pIN_Isr_Isr_int_IRS;
vRP34(13)=ps_IRS*fac.*kRP34(13).*pIIrIrI_int.*IRS; 
vRP34(14)=kRP34(14).*pIIrIrI_int_IRS;
vRP34(15)=ps_IRS*fac.*kRP34(15).*pIN_Isr_Isr_IN_int.*IRS;
vRP34(16)=kRP34(16).*pIN_Isr_Isr_IN_int_IRS;
vRP34=vRP34';



%% vP **
vP(1)=kP(1)*RasT*CRaf;
vP(2)=kP(2)*RasT_CRaf;
vP(3)=kP(3)*RasT*BRaf;
vP(4)=kP(4)*RasT_BRaf;
vP(5)=kP(5)*RasT_CRaf*RasT_CRaf;
vP(6)=kP(6)*RasT_CRaf_CRaf;
vP(7)=kP(7)*RasT_BRaf*RasT_BRaf;
vP(8)=kP(8)*RasT_BRaf_BRaf;
vP(9)=kP(9)*RasT_CRaf*RasT_BRaf;
vP(10)=kP(10)*RasT_CRaf_BRaf;
vP(11)=0;%kP(11)*RasT_BRaf*CRaf;
vP(12)=0;%kP(12)*RasT_CRaf_BRaf;
vP(13)=kP(13)*MEK*RasT_CRaf_CRaf;
vP(14)=kP(14)*MEK_RasT_CRaf_CRaf; 
vP(15)=kP(15)*MEK_RasT_CRaf_CRaf; 
vP(16)=kP(16)*pMEK*RasT_CRaf_CRaf;
vP(17)=kP(17)*pMEK_RasT_CRaf_CRaf;
vP(18)=kP(18)*pMEK_RasT_CRaf_CRaf; 
vP(19)=kP(19)*MEK*RasT_BRaf_BRaf;
vP(20)=kP(20)*MEK_RasT_BRaf_BRaf; 
vP(21)=kP(21)*MEK_RasT_BRaf_BRaf; 
vP(22)=kP(22)*pMEK*RasT_BRaf_BRaf;
vP(23)=kP(23)*pMEK_RasT_BRaf_BRaf; 
vP(24)=kP(24)*pMEK_RasT_BRaf_BRaf; 
vP(25)=kP(25)*MEK*RasT_CRaf_BRaf;
vP(26)=kP(26)*MEK_RasT_CRaf_BRaf;
vP(27)=kP(27)*MEK_RasT_CRaf_BRaf;
vP(28)=kP(28)*pMEK*RasT_CRaf_BRaf;
vP(29)=kP(29)*pMEK_RasT_CRaf_BRaf;
vP(30)=kP(30)*pMEK_RasT_CRaf_BRaf;
vP(31)=kP(31)*pMEK;
vP(32)=kP(32)*ppMEK;
vP(33)=kP(33)*ERK*ppMEK;
vP(34)=kP(34)*ERK_ppMEK; 
vP(35)=kP(35)*ERK_ppMEK; 
vP(36)=kP(36)*pERK*ppMEK;
vP(37)=kP(37)*pERK_ppMEK;
vP(38)=kP(38)*pERK_ppMEK;
vP(39)=kP(39)*ppERKnuc*MKP1; 
vP(40)=kP(40)*ppERKnuc_MKP1;
vP(41)=kP(41)*ppERKnuc_MKP1;
vP(42)=kP(42)*pMKP1; 
vP(43)=kP(43)*pcFos; 
vP(44)=kP(44)*pTSC2; 
vP(45)=kP(45)*pRSKnuc; 
vP(46)=kP(46)*G2_SOS*ppERK; 
vP(47)=kP(47)*G2_SOS_ppERK; 
vP(48)=kP(48)*G2_SOS_ppERK; 
vP(49)=kP(49)*G2_pSOS;
vP(50)=kP(50)*CRaf*ppERK;
vP(51)=kP(51)*CRaf_ppERK;
vP(52)=kP(52)*CRaf_ppERK;
vP(53)=kP(53)*pCRaf;
vP(54)=kP(54)*RasD*DAG_GRP;
vP(55)=kP(55)*RasD_DAG_GRP;
vP(56)=kP(56)*RasD_DAG_GRP;
vP(57)=kP(57)*RasT*NF1;
vP(58)=kP(58)*RasT_NF1;
vP(59)=kP(59)*RasT_NF1;
vP(60)=kP(60)*NF1*ppERK;
vP(61)=kP(61)*NF1_ppERK;
vP(62)=kP(62)*NF1_ppERK;
vP(63)=kP(63)*pNF1;
vP(64)=kP(64)*pERK*MKP3; 
vP(65)=kP(65)*pERK_MKP3;
vP(66)=kP(66)*pERK_MKP3;
vP(67)=kP(67)*ppERK*MKP3; 
vP(68)=kP(68)*ppERK_MKP3;
vP(69)=kP(69)*ppERK_MKP3;
vP(70)=kP(70)*RSK*ppERK;
vP(71)=kP(71)*RSK_ppERK;
vP(72)=kP(72)*RSK_ppERK;
vP(73)=kP(73)*pRSK; 
vP(74)=kP(74)*ppERK; 
vP(75)=kP(75)*ERKnuc; 
vP(76)=kP(76)*ppERKnuc*pMKP1; 
vP(77)=kP(77)*ppERKnuc_pMKP1 ;
vP(78)=kP(78)*ppERKnuc_pMKP1 ;
vP(79)=kP(79)*ERKnuc; 
vP(80)=kP(80)*pRSK; 
vP(81)=kP(81)*MKP1*pRSKnuc; 
vP(82)=kP(82)*pRSKnuc_MKP1;
vP(83)=kP(83)*pRSKnuc_MKP1;
vP(84)=kP(84)*MKP1*ppERKnuc; 
vP(85)=kP(85)*ppERKnuc_MKP1; 
vP(86)=kP(86)*ppERKnuc_MKP1;
vP(87)=kP(87)*cFos*pRSKnuc; 
vP(88)=kP(88)*cFos_pRSKnuc;
vP(89)=kP(89)*cFos_pRSKnuc;
vP(90)=kP(90)*cFos*ppERKnuc; 
vP(91)=kP(91)*cFos_ppERKnuc;
vP(92)=kP(92)*cFos_ppERKnuc;
vP(93)=kP(93)*pcFos*cJun; 
vP(94)=kP(94)*pcFos_cJun;
vP(95)=kP(95)*GRP*DAG;
vP(96)=kP(96)*DAG_GRP;
vP(97)=kP(97)*DAG*PKC;
vP(98)=kP(98)*DAG_PKC;
vP(99)=kP(99)*RKIP*DAG_PKC;
vP(100)=kP(100)*RKIP_DAG_PKC;
vP(101)=kP(101)*RKIP_DAG_PKC;
vP(102)=kP(102)*pRKIP;
vP(103)=kP(103)*RKIP*CRaf;
vP(104)=kP(104)*RKIP_CRaf;
vP(105)=kP(105)*PIP3*PTEN;
vP(106)=kP(106)*PIP3_PTEN;
vP(107)=kP(107)*PIP3_PTEN;
vP(108)=kP(108)*PIP3*AKT;
vP(109)=kP(109)*PIP3_AKT;
vP(110)=kP(110)*PIP3*PDK1;
vP(111)=kP(111)*PIP3_PDK1;
vP(112)=kP(112)*Rictor*mTOR;
vP(113)=kP(113)*mTORC2;
vP(114)=kP(114)*PIP3_AKT*PIP3_PDK1;
vP(115)=kP(115)*PIP3_AKT_PIP3_PDK1;
vP(116)=kP(116)*PIP3_AKT_PIP3_PDK1;
vP(117)=kP(117)*PIP3_pAKT;
vP(118)=kP(118)*PIP3_pAKT*mTORC2;
vP(119)=kP(119)*PIP3_pAKT_mTORC2;
vP(120)=kP(120)*PIP3_pAKT_mTORC2;
vP(121)=kP(121)*PIP3_ppAKT; 
vP(122)=kP(122)*PIP3_ppAKT; 
vP(123)=kP(123)*PIP3*ppAKT;
vP(124)=kP(124)*GSK3b*ppAKT;
vP(125)=kP(125)*GSK3b_ppAKT;
vP(126)=kP(126)*GSK3b_ppAKT;
vP(127)=kP(127)*pGSK3b;
vP(128)=(kP(128)*bCATENIN*GSK3b^4)/(200^4+GSK3b^4);
vP(131)=kP(131)*pbCATENIN;
vP(132)=kP(132)*bCATENIN;
vP(133)=kP(133)*bCATENINnuc;
vP(134)=kP(134)*TSC2*ppAKT;
vP(135)=kP(135)*TSC2_ppAKT;
vP(136)=kP(136)*TSC2_ppAKT;
vP(137)=kP(137)*TSC2*ppERK;
vP(138)=kP(138)*TSC2_ppERK;
vP(139)=kP(139)*TSC2_ppERK;
vP(140)=kP(140)*TSC1*TSC2;
vP(141)=kP(141)*TSC;
vP(142)=0;
vP(143)=0;
vP(144)=0;
vP(145)=0;
vP(146)=kP(146)*mTORC1;
vP(147)=(kP(147)*mTORC1active*TSC^4)/(.1+TSC^4);
vP(148)=kP(148)*Raptor*mTOR;
vP(149)=kP(149)*mTORC1;
vP(150)=kP(150)*EIF4EBP1*mTORC1active;
vP(151)=kP(151)*EIF4EBP1_mTORC1active;
vP(152)=kP(152)*EIF4EBP1_mTORC1active;
vP(153)=kP(153)*pEIF4EBP1;
vP(154)=kP(154)*S6K*mTORC1active;
vP(155)=kP(155)*S6K_mTORC1active;
vP(156)=kP(156)*S6K_mTORC1active;
vP(157)=kP(157)*pS6K;
vP(158)=kP(158)*FOXO*ppAKT;
vP(159)=kP(159)*FOXO_ppAKT;
vP(160)=kP(160)*FOXO_ppAKT;
vP(161)=kP(161)*pFOXO;
vP(162)=kP(162)*PI3K1*mTORC1active;
vP(163)=kP(163)*PI3K1_mTORC1active;
vP(164)=kP(164)*PI3K1_mTORC1active;
vP(165)=kP(165)*pPI3K1;
vP(166)=kP(166)*PI3KC1*PI3KR1;
vP(167)=kP(167)*PI3K1;
vP(168)=kP(168)*PI3P;
vP(169)=kP(169)*SOS*G2;
vP(170)=kP(170)*G2_SOS;
vP(171)=kP(171)*EIF4E*EIF4EBP1;
vP(172)=kP(172)*EIF4EBP1_EIF4E;
vP(173)=kP(173); 
vP(174)=kP(174)*pbCATENIN*TSC; 
vP(175)=kP(175)*pERK; 
vP(176)=kP(176)*ppERK; 
vP(177)=kP(177)*ppERKnuc; 
vP(178)=kP(178)*ppAKT;
vP(179)=kP(179)*pAKT;
vP(180)=kP(180)*RhebT;
vP(181)=kP(181)*FOXO;
vP(182)=kP(182)*FOXOnuc; 
vP(183)=kP(183)*PIP2; 
vP(184)=kP(184)*PIP3; 
vP(185)=kP(185)*RasD; 
vP(186)=kP(186)*RasT; 
vP(187)=kP(187)*MEKi*ppMEK;
vP(188)=kP(188)*MEKi_ppMEK;
vP(189)=kP(189)*AKTi*AKT; 
vP(190)=kP(190)*AKTi_AKT;



%% vDA **
% NONE as these are all transcriptional control mechanisms.

%% vDP **
vDP(1)=kDP(1)*ppAKT*Mdm2;
vDP(2)=kDP(2)*ppAKT_Mdm2;
vDP(3)=kDP(3)*ppAKT_Mdm2;
vDP(4)=kDP(4)*pMdm2;

%% vPA **
vPA(1)=kPA(1)*ppERK*BIM;
vPA(2)=kPA(2)*ppERK_BIM;
vPA(3)=kPA(3)*ppERK_BIM;
vPA(4)=kPA(4)*pBIM;
vPA(5)=kPA(5)*ppAKT*BAD;
vPA(6)=kPA(6)*ppAKT_BAD;
vPA(7)=kPA(7)*ppAKT_BAD;
vPA(8)=kPA(8)*ppERK*BAD;
vPA(9)=kPA(9)*ppERK_BAD;
vPA(10)=kPA(10)*ppERK_BAD;
vPA(11)=kPA(11)*pBAD;

%% vPC **
% NONE as these are all transcriptional control mechanisms.

%% vDC **
% NONE as these are all transcriptional control mechanisms.


%% vDd (none)
% Removed all Dd and Cd degradation reactions. These are controlled from
% within these submodels, respectively.

%% vCd (none) 
% % removed all degradation reactions because cell cycle model handles that
% % internally. 

%% vAd 
Ads=[L_R;
Ractive;
Ractive_flip;
Ractive_pC8;
C8;
C8_Bar;
C8_pC3;
C3;
C3_pC6;
C6;
C6_C8;
C3_XIAP;
C3_PARP;
cPARP;
C8_Bid;
tBid;
tBid_Bcl2c;
tBid_Bax;
Baxactive;
Baxm;
Bcl2;
Baxm_Bcl2;
Bax2;
Bax2_Bcl2;
Bax4;
Bax4_Bcl2;
CytoCr;
Smacr;
CytoC;
CytoC_Apaf;
Apafactive;
Apop;
Apop_C3;
Smac;
Apop_XIAP;
Smac_XIAP;
C3_Ub;
Bcl2c_BAD;
Bcl2c_PUMA;
Bcl2c_NOXA;
BIM_Bax;
Bcl2c_BIM;
ppERK_BIM;
pBIM;
ppAKT_BAD;
pBAD;
ppERK_BAD];

%
Rds=[pE1;
pE2;
pE4;
E1E1;
E1E2;
E1E3;
E1E4;
E2E2;
E2E3;
E2E4;
E3E4;
E4E4;
Met_Met;
FrFr;
IrIr;
Isr_Isr;
EE1;
HE3;
HE4;
HGF_Met;
PPr;
FFr;
EE1E2;
EE1Ev3;
EE1E1;
EE1E3;
EE1E4;
E2HE3;
E1HE3;
HE3E3;
HE3Ev3;
HE3E4;
E2HE4;
HE4Ev3;
E1HE4;
E3HE4;
HE4E4;
HGF_Met_Met;
PPrPr;
FFrFr;
IIrIr;
IN_Isr_Isr;
EE1EE1;
EE1HE3;
EE1HE4;
HE3HE3;
HE3HE4;
HE4HE4;
HGF_Met_HGF_Met;
PPrPPr;
FFrFFr;
IIrIrI;
IN_Isr_Isr_IN;
E1_ppERK;
E2_ppERK;
E4_ppERK];

%
RPds=xRP;

%
Pds=[G2_SOS;
G2_pSOS;
PI3K1;
pPI3K1;
mTORC1;
mTORC1active;
PIP;
PI3P;
DAG;
DAG_GRP;
RasT;
pNF1;
pCRaf;
RasT_CRaf;
RasT_CRaf_BRaf;
pMEK;
ppMEK;
ERKnuc;
ppERKnuc;
pRSK;
pRSKnuc;
pMKP1;
pcFos;
pcFos_cJun;
bCATENINnuc;
pbCATENIN;
IP3;
PIP2;
PIP3;
PIP3_AKT;
pAKT;
ppAKT;
PIP3_PDK1;
PIP3_pAKT;
mTORC2;
PIP3_ppAKT;
pGSK3b;
pTSC2;
TSC;
DAG_PKC;
pRKIP;
RKIP_CRaf;
pERK;
ppERK;
pFOXO;
RhebT;
pS6K;
pEIF4EBP1;
G2_SOS_ppERK;
CRaf_ppERK;
RasD_DAG_GRP;
RasT_NF1;
NF1_ppERK;
MEK_RasT_CRaf_BRaf;
pMEK_RasT_CRaf_BRaf;
ERK_ppMEK;
pERK_ppMEK;
RSK_ppERK;
pRSKnuc_MKP1;
ppERKnuc_MKP1;
cFos_pRSKnuc;
cFos_ppERKnuc;
RKIP_DAG_PKC;
PIP3_PTEN;
PIP3_AKT_PIP3_PDK1;
PIP3_pAKT_mTORC2;
GSK3b_ppAKT;
bCATENIN_GSK3b;
TSC2_ppAKT;
TSC2_ppERK;
RhebT_TSC;
EIF4EBP1_mTORC1active;
S6K_mTORC1active;
FOXO_ppAKT;
PI3K1_mTORC1active;
pERK_MKP3;
ppERK_MKP3;
ppERKnuc_pMKP1;
RasT_BRaf;
RasT_BRaf_BRaf;
MEK_RasT_BRaf_BRaf;
pMEK_RasT_BRaf_BRaf;
EIF4EBP1_EIF4E;
RasT_CRaf_CRaf;
MEK_RasT_CRaf_CRaf;
pMEK_RasT_CRaf_CRaf;
FOXOnuc
MEKi_ppMEK
AKTi_AKT];

Xds = [Ads;Rds;RPds;Pds];

vXd=kXd.*Xds;


%% PUTTING IT TOGETHER  %%

vRP=[vRP1;vRP2;vRP3;vRP4;vRP5;vRP6;vRP7;vRP8;vRP9;vRP10;vRP11;vRP12;vRP13;vRP14;vRP15;vRP16;vRP17;vRP18;vRP19;vRP20;vRP21;vRP22;vRP23;vRP24;vRP25;vRP26;vRP27;vRP28;vRP29;vRP30;vRP31;vRP32;vRP33;vRP34];
v=[vbR;vdR;vTL;vTLCd';vE';vD';vC';vA';vR';vRP;vP';vDP';vPA';vXd];


v(isinf(v))=0;
v(isnan(v))=0;

ndot=S_PARCDL*(v.*(VvPARCDL*1E12));
ydot=ndot./(VxPARCDL*1E12);


flag=0;
new_data = [];

end