function data = distilldata(infile)
%����˵����
%���������ݵ�ԭʼ�ļ��е���ֵ���ݶ��뵽һ��data������
%ʹ��˵����
% infile����ԭʼ�����ļ���;
% data=���ݱ���

tmpfile='tmp2.mat';

fidin=fopen(infile,'r'); % ��ԭʼ�����ļ���.list��

fidtmp=fopen(tmpfile,'w'); % �������������ļ�������˵�����֣�

while ~feof(fidin) % �ж��Ƿ�Ϊ�ļ�ĩβ
    tline=fgetl(fidin); % ���ļ�����һ���ı��������س�����
    if ~isempty(tline) % �ж��Ƿ����
        flag=0;
        if (double(tline(1))>=48&&double(tline(1))<=57) || ...
                (tline(1)=='-') || (tline(1)=='+') || (tline(1)=='.')
            flag=1;
        end
        if flag==1 % ����������У��Ѵ�������д���ļ�
            fprintf(fidtmp,'%s\n',tline);
        end
    end
end

fclose(fidin);

fclose(fidtmp);

data=textread(tmpfile);

delete(tmpfile);

end