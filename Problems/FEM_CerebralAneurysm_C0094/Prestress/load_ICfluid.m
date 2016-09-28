function v = load_ICfluid(x, d)

load U0_Fluid U0_Fluid;

numNodes = length(x);

if d < 4
    v = U0_Fluid(1+(d-1)*numNodes:numNodes*d)';
else
    v = U0_Fluid(1+(d-1)*numNodes:end)';
end

end
