function B = set_up_sub_matrix(G, M)
%set_up_matrix Set up matrix for the finite difference equations

B = zeros(G.N, G.N);

for i=1:G.N
    for j=find(G.contacts(i,:))
        if i ~= j
            B(i,i) = B(i,i) + M(i,j);
            B(i,j) = B(i,j) - M(i,j);
        end
    end
end

B = sparse(B);

end

