function Run_SimCells2(STIM1,STIM2,th1,th2,filename,dataS)

if ~exist('dataS','var')
    dataS=[];
end

matpath='matfiles/';

flagD=0;
St=load([matpath,'RandomPopCells.mat']);
cells0=St.cells0;

for i=1:length(cells0)
    
    xoutG=cells0{i}.xoutG_all;
    xoutS=cells0{i}.xoutS_all;
    [tout_all1,xoutG_all1,xoutS_all1]=RunModel(flagD,th1,STIM1,xoutS,xoutG',dataS,[]);
    
    xoutG=xoutG_all1(end,:);
    xoutS=xoutS_all1(end,:);
    [tout_all2,xoutG_all2,xoutS_all2]=RunModel(flagD,th2,STIM2,xoutS,xoutG',dataS,[]);
    
    tout_all=[tout_all1;(tout_all2(2:end)+tout_all1(end))];
    xoutG_all=[xoutG_all1;xoutG_all2(2:end,:)];
    xoutS_all=[xoutS_all1;xoutS_all2(2:end,:)];
    
    inds=1:30:length(tout_all);
    C.tout_all=tout_all(inds);
    C.xoutG_all=xoutG_all(inds,:);
    C.xoutS_all=xoutS_all(inds,:);

    cells{i,1}=C;
    
    disp(strcat('cell number =',num2str(i)))

end

txt=strcat(matpath,filename,'.mat');
save(txt,'-v7.3','cells');


