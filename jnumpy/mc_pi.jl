function mc_pi(n::Integer)::Float64
    count = 0 
    for _ in 1:n
        x = rand()
        y = rand()
        if x^2 + y^2 < 1
            count += 1
        end
    end
    return 4*count/n
end