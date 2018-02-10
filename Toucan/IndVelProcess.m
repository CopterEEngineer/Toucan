% clear;
% clc;
%%
currentFolder = pwd;
content2={'lambdi'};
content3={'bladedeform','tipgeo','tipgeo_athub'};
suffix='.output';
prefix2='DEBUG';
prefix3={'id0','id1','id2'};
targetfolder = currentFolder;
%%
%两维数据
ni = 72; nj = 40;
for i=1:length(content2)
    var=strcat(prefix2,'_',content2{i});
    targetfile=strcat(var,suffix);
    targetfile=strcat(targetfolder,'\',targetfile);
      
    data=importdata(targetfile,'\t',0);
    eval([var,'=data(1:',num2str(ni),',1:',num2str(nj),');']);
end

%%
%三维数据
nk=1440;
for i=1:length(content3)
    var=strcat(prefix2,'_',content3{i});
    for j=1:length(prefix3)
        targetfile=strcat(prefix3{i},'-',var,suffix);
        targetfile=strcat(targetfolder,'\',targetfile);
        
        data=importdata(targetfile,'\t',0);
    end
end