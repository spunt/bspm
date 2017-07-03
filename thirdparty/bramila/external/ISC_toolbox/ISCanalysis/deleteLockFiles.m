function deleteLockFiles(Params)

qq = ls(Params.PublicParams.dataDestination);
for m = 1:size(qq,1)
    if length(qq(m,:)) > 4
        if strcmp(qq(m,end-3:end),'lock')
            delete([Params.PublicParams.dataDestination qq(m,:)]);
            disp(['Lock-file ' qq(m,:) ' deleted...'])
        end
    end
end
