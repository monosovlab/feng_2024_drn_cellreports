function C = cateachcell(dim,A,B)
    assert(all(size(A)==size(B)),'A and B must be same size');
    C = cell(size(A));
    for ii=1:numel(C)
        C{ii} = cat(dim,A{ii},B{ii});
    end
end