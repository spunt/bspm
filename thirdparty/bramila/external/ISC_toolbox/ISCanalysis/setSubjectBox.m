function handles = setSubjectBox(handles)

disp(' ')
disp('Currently selected subject source files:')
disp(' ')
L = get(handles.editSubj,'String');
%disp(L)

for h = 1:size(handles.Pub.subjectSource,2)
    handles.Pub.subjectSource{get(handles.popupmenuSession,'Value'),h} = [];
end

for h = 1:size(L,1)
    if ischar(L)
        handles.Pub.subjectSource{get(handles.popupmenuSession,'Value'),h} = L(h,:);
    else
        handles.Pub.subjectSource{get(handles.popupmenuSession,'Value'),h} = L{h};
    end
end


flagEmpty = true*ones(1,size(handles.Pub.subjectSource,1));
for s = 1:size(handles.Pub.subjectSource,1)
    for k = 1:size(handles.Pub.subjectSource,2)
        if ~isempty(handles.Pub.subjectSource{s,k})
            flagEmpty(s) = false;
            break
        end
    end
end

F = find(flagEmpty == false);
if ~isempty(F)
    for h = 1:(F(end)+1)
        D{h} = ['Session' num2str(h)];
    end
    set(handles.popupmenuSession,'Value',max(length(D)-1,1))
    for h = 1:size(handles.Pub.subjectSource,2)
        D2{h} = handles.Pub.subjectSource{max(length(D)-1,1),h};
    end
    set(handles.editSubj,'String',D2)    
    set(handles.popupmenuSession,'String',D)
    SubjSour{1} = [];
    for s = 1:F(end)
        for k = 1:size(handles.Pub.subjectSource,2)
            SubjSour{s,k} = handles.Pub.subjectSource{s,k};
        end
    end
    handles.Pub.subjectSource = SubjSour;    
else
    set(handles.popupmenuSession,'Value',1)
    set(handles.popupmenuSession,'String','Session 1')
    handles.Pub.subjectSource = cell(1);
end


for a1 = 1:size(handles.Pub.subjectSource,1)
    disp(['Session ' num2str(a1) ':'])
    disp(' ')
    for a2 = 1:size(handles.Pub.subjectSource,2)
        disp(handles.Pub.subjectSource{a1,a2})
    end
    disp(' ')
end
