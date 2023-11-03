function Am = set_up_cell_area(G)
%Am = set_up_cell_area(G)
Am = zeros(G.N, 1);

for i=1:G.N
    Am(i) = 4*pi*(1e-4*G.R(i))^2;
end

end
