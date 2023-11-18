function Obs = GetObservables_matrix(xout,VxPARCDL,Vc)

% Input all species and this will return a sum across of all the species
% that comprise by each observable.

ObsMat=csvread('observables_mat_18.csv',1,1); ObsMat=ObsMat(:,1:end-1);

xout_corr=xout.*repmat((VxPARCDL'/Vc),size(xout,1),1);

Obs=xout_corr*ObsMat;