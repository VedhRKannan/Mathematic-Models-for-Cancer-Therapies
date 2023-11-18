function Run_SimCells(STIM,th,filename,dataS)

matpath='matfiles/';
[~,dataG]=RunPrep;

flagD=0;
St=load([matpath,'RandomPopCells.mat']);
cells0=St.cells0;

for i=1:length(cells0)
    
    xoutG=cells0{i}.xoutG_all;
    xoutS=cells0{i}.xoutS_all;
    
    if exist('dataS','var')
    [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,xoutS,xoutG',dataS,[]);    
    else
    [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,xoutS,xoutG',[],[]);
    end
    
    inds=1:30:length(tout_all);

    C.tout_all=tout_all(inds);
    C.xoutG_all=xoutG_all(inds,:);
    C.xoutS_all=xoutS_all(inds,:);

    cells{i,1}=C;
    
    disp(strcat('cell number =',num2str(i)))

end

txt=strcat(matpath,filename,'.mat');
save(txt,'-v7.3','cells');


