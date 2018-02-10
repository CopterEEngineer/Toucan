function data = distilldata(infile)
%功能说明：
%将保存数据的原始文件中的数值数据读入到一个data变量中
%使用说明：
% infile——原始数据文件名;
% data=数据变量

tmpfile='tmp2.mat';

fidin=fopen(infile,'r'); % 打开原始数据文件（.list）

fidtmp=fopen(tmpfile,'w'); % 创建保存数据文件（不含说明文字）

while ~feof(fidin) % 判断是否为文件末尾
    tline=fgetl(fidin); % 从文件读入一行文本（不含回车键）
    if ~isempty(tline) % 判断是否空行
        flag=0;
        if (double(tline(1))>=48&&double(tline(1))<=57) || ...
                (tline(1)=='-') || (tline(1)=='+') || (tline(1)=='.')
            flag=1;
        end
        if flag==1 % 如果是数字行，把此行数据写入文件
            fprintf(fidtmp,'%s\n',tline);
        end
    end
end

fclose(fidin);

fclose(fidtmp);

data=textread(tmpfile);

delete(tmpfile);

end