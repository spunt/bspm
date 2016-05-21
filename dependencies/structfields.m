function allnames = structfields(thestruct, name, allnames)
if ~exist('name','var'), name = 'Str'; end
if ~exist('allnames','var'), allnames = []; end

if ischar(thestruct) 
    allnames = [allnames; {sprintf( '%s = %s;\n',name,thestruct)}];
    return;
end

if isstruct(thestruct)
    if length( thestruct ) == 1
        fns = fieldnames(thestruct);
        for n = 1:length(fns)
            fn = getfield(thestruct,fns{n});
            allnames = structfields( fn, [name '.' fns{n}], allnames);
        end
    else
        for k = 1:length( thestruct )
            fns = fieldnames(thestruct(k));
            for n = 1:length(fns)
                fn = getfield(thestruct(k),fns{n});
                allnames = structfields( fn, [name '(' int2str(k) ').' fns{n}], allnames );
            end
        end
    end
    return;
end

if iscell(thestruct)
    for n = 1:length(thestruct)
        allnames = structfields( thestruct{n}, [name '{' int2str(n) '}'], allnames );
    end
    return;
end

if isnumeric(thestruct)
    if length(thestruct) == 1
        allnames = [allnames; {sprintf('%s = %s;\n',name,num2str(thestruct))}];
    else
        for n = 1:length(thestruct)
            allnames = [allnames; {sprintf( '%s[%d] = %s;\n',name,n,num2str(thestruct(n)))}];
        end
    end
    return;
end





