function adj = read_svml(filepath)
    svml = table2array(readtable(filepath,"FileType",'text',"Delimiter",'\t','NumHeaderLines',0));
    r = [];
    c = [];
    for i = 1:size(svml,1)
        neighbors = svml(i, ~isnan(svml(i,:)));
        neighbors = neighbors(2:end) + 1;
        r = [r repelem(i,size(neighbors,2))];
        c = [c neighbors];
    end
    adj = sparse(r,c,repelem(1,length(r)));
end