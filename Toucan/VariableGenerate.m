function variable = VariableGenerate(prefix, content, postfix)

prefix = strrep(prefix,'-','_');
content = strrep(content,'-','_');
postfix = strrep(postfix,'-','_');
variable = strcat(prefix, '_', content,'_',postfix);

end
