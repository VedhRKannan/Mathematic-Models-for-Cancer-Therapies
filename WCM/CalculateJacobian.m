function CalculateJacobian

% This function calculates the Jacobian Matrix for the model. This must be
% re-run if the model structure is altered.


[dataS,dataG]=RunPrep;

x = sym('x',[length(dataS.x0PARCDL),1]);

k = sym('k',[length(dataS.kS),1]);
dataS.kS=k;

Vv = sym('Vv',[length(dataS.VvPARCDL),1]);
dataS.VvPARCDL=Vv;

Vx = sym('Vx',[length(dataS.VxPARCDL),1]);
dataS.VxPARCDL=Vx;

mExp = sym('mExp',[length(dataS.mExp_nM),1]);
dataS.mExp_nM=mExp;

mMod = sym('mMod',[length(dataS.mMod),1]);
dataS.mMod=mMod;

t=0;

[ydot,flag,new_data,v] = createODEs(t,x,dataS);

J=jacobian(ydot,x);
[SPr,SPc]=find(J);

%Writing Jeval.m file
disp('... Writing Jeval.m File')
ExportJactoTxt(J,SPr,SPc,dataS); %Write Jeval.m file; this file will be the input to the ode solver.

end


function ExportJactoTxt(J,SPr,SPc,data)

% This script creates the Jeval.m matlab script using text functions. It
% overwrites an existing file if it exists. The created file, Jeval.m, will
% be used as input to the options for the ode15s solver. It only writes out
% parts of the Jacobian that are non-zero, as indicated by SPr and SPc
% inputs.

n=size(J,1);
filename=strcat('Jeval',num2str(n));

FID = fopen(strcat(filename,'.m'),'w');
formatSpec = strcat('function [Je,flag,new_data]=',filename,'(t,x,fy,data)\n\n');
fprintf(FID,formatSpec);

FID = fopen(strcat(filename,'.m'),'a');
formatSpec = 'flag=0; new_data=[];\n';
fprintf(FID,formatSpec);

FID = fopen(strcat(filename,'.m'),'a');
formatSpec = 'Je=zeros(length(x));\n';
fprintf(FID,formatSpec);

FID = fopen(strcat(filename,'.m'),'a');
formatSpec = 'k=data.kS;\n';
fprintf(FID,formatSpec);

%added
FID = fopen(strcat(filename,'.m'),'a');
formatSpec = 'Vv=data.VvPARCDL;\n';
fprintf(FID,formatSpec);

FID = fopen(strcat(filename,'.m'),'a');
formatSpec = 'Vx=data.VxPARCDL;\n';
fprintf(FID,formatSpec);

FID = fopen(strcat(filename,'.m'),'a');
formatSpec = 'mExp=data.mExp_nM;\n';
fprintf(FID,formatSpec);

FID = fopen(strcat(filename,'.m'),'a');
formatSpec = 'mMod=data.mMod;\n';
fprintf(FID,formatSpec);
%done adding


%print x1=x(1) type values
disp('... Writing x1=x(1) type values')
formatSpec = 'x%i=x(%i);\n';
for i=1:size(J,1)
    fprintf(FID,formatSpec,[i,i]);
end

%print k1=k(1) type values
disp('... Writing k1=k(1) type values')
formatSpec = 'k%i=k(%i);\n';
for i=1:length(data.kS)
    fprintf(FID,formatSpec,[i,i]);
end

%adding
disp('... Writing Vv1=Vv(1) type values')
formatSpec = 'Vv%i=Vv(%i);\n';
for i=1:length(data.VvPARCDL)
    fprintf(FID,formatSpec,[i,i]);
end
disp('... Writing Vx1=Vx(1) type values')
formatSpec = 'Vx%i=Vx(%i);\n';
for i=1:length(data.VxPARCDL)
    fprintf(FID,formatSpec,[i,i]);
end
disp('... Writing mExp1=mExp(1) type values')
formatSpec = 'mExp%i=mExp(%i);\n';
for i=1:length(data.mExp_nM)
    fprintf(FID,formatSpec,[i,i]);
end
disp('... Writing mMod1=mMod(1) type values')
formatSpec = 'mMod%i=mMod(%i);\n';
for i=1:length(data.mMod)
    fprintf(FID,formatSpec,[i,i]);
end
%done adding


%print Jacobian equations
disp('... Writing Jacobian Equations')
FID = fopen(strcat(filename,'.m'),'a');
formatSpec = 'Je(%i,%i)=%s;\n';

for i=1:length(SPr)
    str=char(J(SPr(i),SPc(i)));
    fprintf(FID,formatSpec,[SPr(i),SPc(i),str]);
    %disp(i/length(SPr))
end

end


