function file = TargetFile(prefix, content, suffix)
content = strcat(prefix, '-', content);
file = strcat(content, suffix);
end