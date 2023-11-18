function [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,xoutS,xoutG,dataS,dataG,kTCleak,kTCmaxs)

% This function runs the model and outputs timecourse simulation results.

% Required Inputs:
% flagD: 1 for deterministic simulations, 0 for stochastic simulations.
% th: simulation time (hours)
% STIM: stimulus vector

% Outputs:
% tout_all: n-by-1 vector of time values (seconds)
% xoutG_all: n-by-g matrix of species (g) through time (n) (g indices lines up to gm tab in Names.xls sheet) 
% xoutS_all: n-by-p matrix of speices (p) through time (n) (p indices lines up to PARCDL tab in Names.xls sheet) 


%% PREP
if isempty(dataS)
    [dataS,~]=RunPrep;
end
if isempty(dataG)
    [~,dataG]=RunPrep;
end

%% RUN

ts=dataS.ts;
ts_up=ts;
N_STEPS=th*3600/ts;

% IMPORT INITIALIZED PARAMETERS
pathi='initialized/'; 

% for PARCDL
kbR0=dlmread(strcat(pathi,'i_kbR0.txt'));
kTL=dlmread(strcat(pathi,'i_kTLF.txt'));
kC173=dlmread(strcat(pathi,'i_kC173.txt'));
kC82=dlmread(strcat(pathi,'i_kC82.txt'));
kA77=dlmread(strcat(pathi,'i_kA77.txt'))*5; 
kA87=dlmread(strcat(pathi,'i_kA87.txt'));
Rt=dlmread(strcat(pathi,'i_Rt.txt'));
EIF4Efree=dlmread(strcat(pathi,'i_EIF4Efree.txt'));
kDDbasal=dlmread(strcat(pathi,'i_kDDbasal.txt'));
Vc=dataS.kS(3);
% for gm
if ~exist('kTCleak','var')
    kTCleak=dlmread(strcat(pathi,'i_kTCleakF.txt'));
end
if ~exist('kTCmaxs','var')
    kTCmaxs=dlmread(strcat(pathi,'i_kTCmaxsF.txt'));
end
% modifying data.S structure
dataS.kS(1)=Rt;
dataS.kS(2)=EIF4Efree;
dataS.kS(12)=kbR0;
dataS.kS(17:157)=kTL;
dataS.kS(632)=kC173;
dataS.kS(541)=kC82;
dataS.kS(709)=kA77;
dataS.kS(719)=kA87;
dataS.kS(450)=kDDbasal;
% modifying data.G structure
dataG.kTCleak=kTCleak;
dataG.kTCmaxs=kTCmaxs;

%species
if isempty(xoutS)
    xoutS=dlmread(strcat(pathi,'i_xoutF.txt')); xoutS=xoutS(end,:);
end

if isempty(xoutG)
if flagD
    xoutG=dataG.x0gm_mpc_D;
else
    xoutG=dataG.x0gm_mpc;
    indsD=dataG.indsD;
    xoutG([indsD,indsD+141,indsD+141*2])=dataG.x0gm_mpc_D([indsD,indsD+141,indsD+141*2]);
end
end


% Apply STIM
Etop=STIM(end);
STIM=STIM(1:end-1);
xoutS(logical(STIM))=STIM(logical(STIM));
dataS.kS(453)=Etop;


% Instantiation
t0 = 0;
optionscvodes = CVodeSetOptions('UserData', dataS,...
                          'RelTol',1.e-3,...
                          'LinearSolver','Dense',...
                          'JacobianFn',@Jeval774);                      
CVodeInit(@createODEs, 'BDF', 'Newton', t0, xoutS', optionscvodes);

%ODE15s options
%optionsode15s=odeset('RelTol',1e-3,'Jacobian',@Jeval774ode15s);

tout_all=zeros(N_STEPS+1,1);
xoutG_all=zeros(N_STEPS+1,length(xoutG));
xoutS_all=zeros(N_STEPS+1,length(xoutS));
tout_all(1)=0;
xoutG_all(1,:)=xoutG';
xoutS_all(1,:)=xoutS';

% Starting simulations
disp('... Starting Sims')
tic;
for i=1:N_STEPS
    
    %gm
    [xginN,xgacN,AllGenesVecN,xmN,~] = gm(flagD,dataG,ts,xoutG,xoutS);
    xoutG=[xgacN;xginN;xmN];
    dataS.mMod=xmN*(1E9/(Vc*6.023E+23)); %convert mRNAs from mpc to nM
    dataG.AllGenesVec=AllGenesVecN;
        
    %CVODEs PARCDL
    CVodeSet('UserData',dataS);
    [status,tout,xoutS] = CVode(ts_up,'Normal');
    %collect data
    tout_all(i+1)=tout;
    xoutG_all(i+1,:)=xoutG';
    xoutS_all(i+1,:)=xoutS';
    
    if xoutS(104)<xoutS(106)
        disp('Apoptosis happened'); 
        tout_all=tout_all(1:i+1);
        xoutG_all=xoutG_all(1:i+1,:);
        xoutS_all=xoutS_all(1:i+1,:);
        return
    end
    
%     %ODE15s
%     [tout,xoutS]=ode15s(@createODEs,[ts_up-ts ts_up],xoutS_all(i,:),optionsode15s,dataS);
%     tout_all(i+1)=tout(end);
%     xoutG_all(i+1,:)=xoutG';
%     xoutS_all(i+1,:)=xoutS(end,:);    

    
    ts_up=ts_up+ts;
    
    if rem(i,1000)==0; disp(strcat(num2str(i),'...')); end
    
end

CVodeFree;
toc;







