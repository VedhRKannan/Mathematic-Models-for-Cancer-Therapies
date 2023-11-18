function Obs = GetObservables(xout,VxPARCDL,Vc)

% Input all species and this will return a sum across of all the species
% that comprise by each observable.

ObsMat=csvread('observables_mat_18.csv',1,1);

for i=1:size(ObsMat,2)
    Obs(i)=sum(ObsMat(:,i).*xout'.*(VxPARCDL/Vc));
end

end