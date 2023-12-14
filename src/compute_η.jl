function η_l(l, ε_r1, ε_r2)
    if l % 2 == 0
        return 1.0
    else
        return ε_r2 / ε_r1
    end
end

function calculate_η(l_max, ε_r1, ε_r2)
    return [η_l(l, ε_r1, ε_r2) for l in 0:l_max-1]
end