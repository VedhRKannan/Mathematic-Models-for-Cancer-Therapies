function RunRandomPopCells(numcells,filename,STIM)

cells0=[];
flagD=0;
th=24;
if ~exist('STIM','var')
    STIM=zeros(775,1);
end

for i=1:numcells
    [tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],[],[]);
    C.tout_all=tout_all;
    C.xoutG_all=xoutG_all(end,:);
    C.xoutS_all=xoutS_all(end,:);
    cells0{i}=C;
    disp(i)
end

if exist('filename','var')
    save(filename,'-v7.3','cells0');
else
    save('RandomPopCells.mat','-v7.3','cells0');
end
