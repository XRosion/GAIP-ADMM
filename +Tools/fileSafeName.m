function safeName = fileSafeName(name)
%FILESAFENAME Convert display text into a conservative file name stem.

safeName = regexprep(lower(char(name)), '[^a-z0-9]+', '_');
safeName = regexprep(safeName, '^_+|_+$', '');
if isempty(safeName)
    safeName = 'item';
end
end
