function [kTLCd,kTL,kXd,xp_mpc] = PARCDL_rateconstants_sd(ke50,mExp_mpc,kTLnat,kTLd,EIF4Efree,S_TL,S_d)

ObsMat=csvread('observables_mat_18.csv',1,1);
ObsMat=ObsMat(:,1:end-1);

xpi=(kTLnat.*mExp_mpc)./kTLd;
xpi(isnan(xpi))=0;


%% kTL **
temp=1/((EIF4Efree/(ke50+EIF4Efree)));
kTL=kTLnat.*temp;

%% kTLCd

xp_mpc=(kTLnat.*mExp_mpc)./kTLd;
xp_mpc(isnan(xp_mpc))=0;

inds=find(sum(S_TL,2));
for i=1:length(inds)
    genes2include=find(S_TL(inds(i),:));
    kTLCd(i)=sum(kTLd(genes2include).*(xp_mpc(genes2include)/sum(xp_mpc(genes2include)))); %Sum of the fraction of each gene
end

kTLCd(isnan(kTLCd))=0;
kTLCd=kTLCd';



%% kXd

for i=1:size(S_d,2)
    ProtInd=find(S_d(:,i)==-1);
    Obs2Include=find(ObsMat(ProtInd,:));
    if Obs2Include
        kXd(i)=max(kTLCd(Obs2Include));
    else kXd(i)=0;
    end
end


a=0;
kAd=kXd(a+1:a+47);a=a+47;
kRd=kXd(a+1:a+56);a=a+56;
kRPd1=kXd(a+1:a+33);a=a+33;
kRPd2=kXd(a+1:a+29);a=a+29;
kRPd3=kXd(a+1:a+29);a=a+29;
kRPd4=kXd(a+1:a+33);a=a+33;
kRPd5=kXd(a+1:a+29);a=a+29;
kRPd6=kXd(a+1:a+29);a=a+29;
kRPd7=kXd(a+1:a+29);a=a+29;
kRPd8=kXd(a+1:a+29);a=a+29;
kRPd9=kXd(a+1:a+29);a=a+29;
kRPd10=kXd(a+1:a+29);a=a+29;
kRPd11=kXd(a+1:a+29);a=a+29;
kRPd12=kXd(a+1:a+29);a=a+29;
kRPd13=kXd(a+1:a+29);a=a+29;
kRPd14=kXd(a+1:a+29);a=a+29;
kPd=kXd(a+1:a+89);


% Rate constant modified for internalized receptors
kRPd3(3:4)=8.3711E-4; kRPd3([1,2,5:8,11,16])=2.1301E-4; kRPd3([9,10,12:15,17:29])=8.7106E-5; 
kRPd4(3:4)=8.3711E-4; kRPd4([1,2,5:8,11,16])=2.1301E-4; kRPd4([9,10,12:15,17:33])=8.7106E-5; 
kRPd6(3:4)=8.3711E-4; kRPd6([1,2,5:8,11,16])=2.1301E-4; kRPd6([9,10,12:15,17:29])=8.7106E-5; 
kRPd10(3:4)=8.3711E-4; kRPd10([1,2,5:8,11,16])=2.1301E-4; kRPd10([9,10,12:15,17:29])=8.7106E-5; 

% Rate constants modified for degradation of lipids
kPd(7)=0.00001; %PIP
kPd(8)=0.00001; %PI3P
kPd(9)=0.00108*10; %DAG
kPd(27)=0.081; %IP3
kPd(28)=4.83E-5; %PIP2
kPd(29)=4.83E-5; %PIP3

%Other rate constant modifications
kPd(22:24)=log(2)/4/3600; %pMPK1 (DUSP), pcFos, pcFos_cJun


kXd=[kAd';kRd';kRPd1';kRPd2';kRPd3';kRPd4';kRPd5';kRPd6';kRPd7';kRPd8';kRPd9';kRPd10';kRPd11';kRPd12';kRPd13';kRPd14';kPd'];


