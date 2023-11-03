function GJ_areas = set_up_GJ_areas(G)
%GJ_areas = set_up_GJ_areas(G)
GJ_areas = zeros(size(G.contacts));

for i=1:G.N
    for j=find(G.contacts(i,:))
        if j ~= i
            l = norm(G.coords(i,:)-G.coords(j,:)); % Distance between cell centers
            ri = G.R(i);
            rj = G.R(j);
            d = l - ri - rj; % Distance between cells
            GJ_areas(i,j) = 1/(1+d);
        end
    end
end

end