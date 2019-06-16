%=====sub funtion ==================
function x = setparam(opts,name,default)
if isfield(opts,name);
    x = opts.(name);
else
    x = default;
end
end
