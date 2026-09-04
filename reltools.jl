module RelTools

using DoubleFloats

function get_interval(vector)
    c = Double64(3e+8)

    return -(vector[1])^2 + vector[2]^2 + vector[3]^2 + vector[4]^2
end

function is_timelike(S)
    c = Double64(3e+8)
    interval = -(c*S[1])^2 + S[2]^2 + S[3]^2 + S[4]^2

    if interval < 0
        return true
    else
        return false
    end
end

function is_spacelike(S)
    c = Double64(3e+8)
    interval = -(c*S[1])^2 + S[2]^2 + S[3]^2 + S[4]^2

    if interval > Double64(0)
        return true
    else
        return false
    end
end

function is_lightlike(S)
    c = Double64(3e+8)
    interval = -(c*S[1])^2 + S[2]^2 + S[3]^2 + S[4]^2
    #print("Interval: $(interval)\n\n")
    if isapprox(interval, Double64(0), atol=Double64(1e-5))
        return true
    else
        return false
    end
end

function check_constellation(X)
    X1,X2,X3,X4 = X[:,1], X[:,2], X[:,3], X[:,4]
    S1 = X2 - X1
    S2 = X3 - X1
    S3 = X4 - X1
    S4 = X3 - X2
    S5 = X4 - X2
    S6 = X4 - X3

    separations = (S1, S2, S3, S4, S5, S6)

    for separation in separations
        if is_lightlike(separation)
            return false
        elseif is_timelike(separation)
            return false
        end
    end

    return true
end

function check_solution(Xc, X)
    X1,X2,X3,X4 = X[:,1], X[:,2], X[:,3], X[:,4]
    S1 = X1 - Xc
    S2 = X2 - Xc
    S3 = X3 - Xc
    S4 = X4 - Xc

    separations = (S1, S2, S3, S4)

    for separation in separations
        if !is_lightlike(separation)
            return false
        end
    end

    return true
end



end