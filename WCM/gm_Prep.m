function [xgac_mpc,xgin_mpc,xgac_mpc_D,xgin_mpc_D,kTCleak,AllGenesVec,GenePositionMatrix] = gm_Prep(mExp_mpc,gExp_mpc,kTCd,kGac,kGin,numberofgenes)


%% Make Gene Position Matrix 
GenePositionMatrix=zeros(numberofgenes,sum(gExp_mpc));
ind=1;
for i=1:numberofgenes
    GenePositionMatrix(i,ind:ind+gExp_mpc(i)-1)=1;
    ind=ind+gExp_mpc(i);
end

%% xg Deterministic
xgac_mpc_D=(kGac.*gExp_mpc)./(kGin+kGac); %active genes initial condition
xgin_mpc_D=gExp_mpc-xgac_mpc_D; %inactive genes initial condition

%% xg Stochastic
AllGenesVec=zeros(sum(gExp_mpc),1); 
IndsGenesOn=randsample(sum(gExp_mpc),round(sum(gExp_mpc)*kGac(1)/kGin(1))); 
AllGenesVec(IndsGenesOn)=1;
% Calculate Concentration of Active and Inactive Genes for each gene
xgac_mpc=GenePositionMatrix*AllGenesVec;
xgin_mpc=gExp_mpc-xgac_mpc;

%% kTCleak (deterministic)
kTCleak=((kTCd.*mExp_mpc)./xgac_mpc_D);
kTCleak(isnan(kTCleak))=0;
kTCleak(isinf(kTCleak))=0;

