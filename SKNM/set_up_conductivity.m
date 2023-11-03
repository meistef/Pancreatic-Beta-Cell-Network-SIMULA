function Mi = set_up_conductivity(G)
%Mi = set_up_conductivities(G)
Mi = zeros(G.N, G.N);

for i=1:G.N
    for j=find(G.contacts(i,:))
        if i ~= j
            Mi(i, j) = G.Gg(i,j);
        end
    end
end

end