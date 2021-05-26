function d = dir_ls(path)
a = dir(path);
for dd = 1:numel(a)
   ok(dd) = 1;
   if strcmp(a(dd).name,'.') || strcmp(a(dd).name,'..') || strcmp(a(dd).name,'.DS_Store')
       ok(dd) = 0;
   end
end
d = a(ok==1);