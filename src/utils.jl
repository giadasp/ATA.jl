function mean(x::Vector{Real})
    n = size(x, 1)
    if n > 0
        return sum(x) / n
    else
        return 0
    end
end

function _exp_c(x::Float64)
    ccall(:exp, Float64, (Float64,), x)
end

function _log_c(x::Float64)
    ccall(:log, Float64, (Float64,), x)
end

function _log1p_c(x::Float64)
    ccall(:log1p, Float64, (Float64,), x)
end

function _sig_c(x::Float64)
    1 / (1 + _exp_c(-x))
end

function _sig_cplus(x::Float64)
    1 / (1 + _exp_c(x))
end

function _gemmblas(A::Matrix{Float64}, B::Matrix{Float64}, ol_max::Matrix{Float64})
    sa = size(A)
    ol = copy(ol_max)
    ccall(
        ("dgemm_64_", "libopenblas64_"),
        Cvoid,
        (
            Ref{UInt8},
            Ref{UInt8},
            Ref{Int64},
            Ref{Int64},
            Ref{Int64},
            Ref{Float64},
            Ptr{Float64},
            Ref{Int64},
            Ptr{Float64},
            Ref{Int64},
            Ref{Float64},
            Ptr{Float64},
            Ref{Int64},
        ),
        'T',
        'N',
        sa[2],
        sa[2],
        sa[1],
        1.0,
        A,
        max(1, stride(A, 2)),
        B,
        max(1, stride(B, 2)),
        -one(Float64),
        ol,
        max(1, stride(ol, 2)),
    )
    return ol::Matrix{Float64}
end

function _gemvblas(
    A::Matrix{Float64},
    x::Vector{Float64},
    cons::Vector{Float64},
    sizex::Int64,
)
    ccall(
        ("dgemv_64_", "libopenblas64_"),
        Cvoid,
        (
            Ref{UInt8},
            Ref{Int64},
            Ref{Int64},
            Ref{Float64},
            Ptr{Float64},
            Ref{Int64},
            Ptr{Float64},
            Ref{Int64},
            Ref{Float64},
            Ptr{Float64},
            Ref{Int64},
        ),
        'N',
        size(A, 1),
        sizex,
        1.0,
        A,
        max(1, stride(A, 2)),
        x,
        1,
        -1.0,
        cons,
        1,
    )
    return cons
end

function _gemvblasT(A::Matrix{Float64}, x::Vector{Float64}, sizeA2::Int64)
    cons = zeros(Float64, sizeA2)
    ccall(
        ("dgemv_64_", "libopenblas64_"),
        Cvoid,
        (
            Ref{UInt8},
            Ref{Int64},
            Ref{Int64},
            Ref{Float64},
            Ptr{Float64},
            Ref{Int64},
            Ptr{Float64},
            Ref{Int64},
            Ref{Float64},
            Ptr{Float64},
            Ref{Int64},
        ),
        'T',
        size(A, 1),
        size(A, 2),
        1.0,
        A,
        max(1, stride(A, 2)),
        x,
        1,
        -1.0,
        cons,
        1,
    )
    return cons
end

function _combinations(T::Int64)
    comb = vcat([[[t, i] for i = (t+1):T] for t = 1:T]...)
end

function _myqle(x::Vector{Float64}, lx::Int64, k::Int64)
    y = ones(k) * Inf
    n = 1
    y[1] = copy(x[1])
    for i = 2:lx
        current = copy(x[i])
        k2 = 1
        while k2 <= k
            if y[k2] > current
                y[(k2+1):n] .= y[k2:(n-1)]
                y[k2] = copy(current)
                k2 = k + 1
                if n < k
                    n += 1
                end
            end
            k2 += 1
        end
    end
    return y[k]
end

function _mycopy(NH::Neighbourhood, NH_new::Neighbourhood)
    NH_new.x = copy(NH.x)
    NH_new.f = copy(NH.f)
    NH_new.obj = copy(NH.obj)
    NH_new.infeas = copy(NH.infeas)
    NH_new.ol = copy(NH.ol)
    NH_new.iu = copy(NH.iu)
    return NH_new::Neighbourhood
end

function _comp_f(NH::Neighbourhood, OptFeas::Float64)
    return (1 - OptFeas) * (sum(NH.infeas + NH.ol) + NH.iu) - OptFeas * (minimum(NH.obj))
end

function _myqle_simp(x::Vector{Float64}, ind::Vector{Float64})
    R = size(x, 1)
    alphaR = Int.(ceil.(R .* (ind)))
    alphaR[findall(alphaR .< 1)] .= 1
    alphaR[findall(alphaR .> R)] .= R
    αQle = Inf
    xsorted = sort(x)
    return sort(x)[alphaR]
end

function _cutR(
    x;
    start = "minimum",
    stop = "maximum",
    nBins = 2,
    returnBreaks = true,
    returnMidPts = false,
)
    if (start == "minimum")
        start = minimum(x)
    end
    if (stop == "maximum")
        stop = maximum(x)
    end
    bw = (stop - start) / (nBins - 1)
    midPts = zeros(nBins)
    for i = 1:nBins
        midPts[i] = start + (i - 1) * bw
    end
    breaks = collect(range(start - (bw / 2); length = nBins + 1, stop = stop + (bw / 2)))
    y = zeros(size(x, 1))
    for j = 1:size(x, 1)
        for i = 1:nBins
            if (x[j] >= breaks[i]) && (x[j] < breaks[i+1])
                y[j] = i
            end
            if i == nBins && x[j] == breaks[i+1]
                y[j] = i
            end
        end
    end
    if (returnBreaks == true || returnMidPts == true)
        if returnMidPts == false
            return (Int.(y), breaks)
        elseif returnBreaks == false
            return (Int.(y), midPts)
        else
            return (Int.(y), breaks, midPts)
        end
    else
        return Int.(y)
    end
end

function _desume_pars(
    pars::DataFrames.DataFrame;
    model = "2PL", #1PL, 2PL, 3PL
    parametrization = "at-b",
) #"at-b, at-ab, at+b, at+ab")
    n_items = size(pars, 1)
    if (model == "1PL")
        ("b" in names(pars)) ||
            ("d" in names(pars)) ||
            error("difficulty parameter b or d not defined")
        c = zeros(Float64, n_items)
        a = ones(Float64, n_items)
    elseif (model == "2PL")
        ("a" in names(pars) || "a1" in names(pars)) ||
            error("discrimination parameter a1 or a not defined")
        ("b" in names(pars)) ||
            ("d" in names(pars)) ||
            error("difficulty parameter b or d not defined")
        if "a" in names(pars)
            a = pars.a
        else
            a = pars.a1
        end
        c = zeros(Float64, n_items)
    elseif (model == "3PL")
        ("a1" in names(pars) || "a" in names(pars)) ||
            error("discrimination parameter a1 not defined")
        ("b" in names(pars)) ||
            ("d" in names(pars)) ||
            error("difficulty parameter b or d not defined")
        ("c" in names(pars)) ||
            ("g" in names(pars)) ||
            error("guessing parameter c or g not defined")
        if "a" in names(pars)
            a = pars.a
        else
            a = pars.a1
        end
        if ("c" in names(pars))
            c = pars.c
        else
            c = pars.g
        end
    elseif (model == "grm")
        ("a1" in names(pars) || "a" in names(pars)) ||
            error("discrimination parameter a1 not defined")
    end

    nb = 0
    b = zeros(Float64, size(pars, 1), 1)
    for n in names(pars)
        if startswith(string(n), "b")
            nb += 1
            b = hcat(b, pars[!, n])
        end
    end
    nd = 0
    b = zeros(Float64, size(pars, 1), 1)
    for n in names(pars)
        if startswith(string(n), "d")
            nd += 1
            d = hcat(d, pars[!, n])
        end
    end
    # TODO grm
    if nd > 0
        b = mapslices(x -> b .- x, d; dims = 2)
    end

    a2 = ones(n_items) #default: "at+b" a*theta + b
    if parametrization == "at-b" #a*theta - b
        b = -b
    elseif parametrization == "at-ab" #a*(theta - b)
        a2 = copy(a)
        b = -b
    elseif parametrization == "at+ab" #a*(theta + b)
        a2 = copy(a)
    end 
    return a, a2, b, c
end

"""
	FS_to_items(xₜ, n_items, FS_items)

# Description

Ungroup items grouped by friend sets. 

# Arguments

- **`xₜ::Vector{Float64}`** : Required. grouped items 0-1 vector.
- **`n_items::Int64`** : Required. Number if items.
- **`FS_items::Vector{Vector{Int64}}`** : Required. vector of items included in the friend sets.

# Output

- **`x_Iₜ::Vector{Float64}`** Returns a 0-1 vector at item level.

"""
function FS_to_items(xₜ::Vector{Float64}, n_items::Int64, FS_items::Vector{Vector{Int64}})
    xₜ_taken = findall(xₜ .== one(Float64))
    x_Iₜ = zeros(Float64, n_items)
    for i in xₜ_taken
        x_Iₜ[FS_items[i]] .= one(Float64)
    end
    return x_Iₜ::Vector{Float64}
end

"""
	item_char(pars, theta; model = "2PL", parametrization = "at-b", D = 1, derivatives = false) 

# Description

Compute the item characteristic function (probability of correct answer).

# Arguments

- **`pars::DataFrames.DataFrame`** : Required. DataFrame with item parameters (difficulty: b or d, discrimination: a or a1, guessing: c or g).
- **`theta::Vector{Float64}`** : Required. Vector of ability points.
- **`model`** : Optional. Default: `"2PL`". Values:  `"1PL`", `"2PL`", `"3PL`". IRT model.
- **`parametrization`** : Optional. Default: `"at-b`". Values:  `"at-b`", `"at-ab`", `"at+b`", `"at+ab`". IRT model parametrization. Ex: at-b is ``Da(\\theta-b)``.
- **`D`** : Optional. Default: 1. 
- **`derivatives`** : Optional. Default: false. If it's true compute the derivatives. Ow a `zeros(0,0)` matrix is returned. 

# Output

- **`p::Matrix{Float64}`**: Matrix of probabilites. 
- **`pder::Matrix{Float64}`**: Matrix of derivatives.
"""
function item_char(
    pars::DataFrames.DataFrame,
    theta::Vector{Float64};
    model = "2PL", #1PL, 2PL, 3PL
    parametrization = "at-b", #"at-b, at-ab, at+b, at+ab"
    D = 1,#any number (1, 1.6)
    derivatives = false,
)

    n_items = size(pars, 1)
    pder = zeros(0, 0)
    a, a2, b, c = _desume_pars(pars; model = model, parametrization = parametrization)
    nb = size(b, 2)
    N = size(theta, 1)
    p =
        _sig_c.([
            D * (a[i] * theta[n] + a2[i] * b[i, j]) for i = 1:n_items, n = 1:N, j = 1:nb
        ])
    #pder =  eachslice(((1 .- p) .* p), dims = 1) .* a
    if derivatives
        pder = mapslices(x -> (1 .- x) .* x .* a, p; dims = 1)
    end

    if model != "grm"
        p = c .+ ((1 .- c) .* p)
    else
        # TODO grm
        #generalized PCM, does not work
        # p[:, :, 1] = 1 .- p[:, :, 1]
        # pder[:, :, 1] = .- pder[:, :, 1]
        # for k = (nb - 1) :  -1 : 2
        # 	p[:, :, k] = p[:, :, k] - p[:, :, k-1]
        # 	pder[:, :, k] = pder[:, :, k] - pder[:, :, k-1]
        # end
    end
    return p, pder
end

"""
	item_info(pars, theta; model = "2PL", parametrization = "at-b", D = 1, derivatives = false) 

# Description

Compute the item Fisher information function.

# Arguments

- **`pars::DataFrames.DataFrame`** : Required. DataFrame with item parameters (difficulty: b or d, discrimination: a or a1, guessing: c or g).
- **`theta::Vector{Float64}`** : Required. Vector of ability points.
- **`model`** : Optional. Default: `"2PL`". Values:  `"1PL`", `"2PL`", `"3PL`". IRT model.
- **`parametrization`** : Optional. Default: `"at-b`". Values:  `"at-b`", `"at-ab`", `"at+b`", `"at+ab`". IRT model parametrization. Ex: at-b is ``Da(\\theta-b)``.
- **`D`** : Optional. Default: 1. 

# Output

- **`i::Matrix{Float64}`**: Matrix of the item information. 
"""
function item_info(
    pars::DataFrames.DataFrame,
    theta;
    model = "2PL", #1PL, 2PL, 3PL, grm
    parametrization = "at-b", #"at-b, at-ab, at+b, at+ab"
    D = 1,
) #true, false

    n_items = size(pars, 1)
    a, a2, b, c = _desume_pars(pars; model = model, parametrization = parametrization)
    p, pder = item_char(pars, theta; model = model, parametrization = parametrization)
    nb = size(p, 3)
    if model != "grm"
        i = (a .^ 2) .* ((1 .- p) ./ p) .* ((p .- c) ./ (1 .- c)) .^ 2
    else
        i = pder .^ 2 ./ p
        i = sum(i, dims = 3)[:, :, 1]
    end
    return i
end


function student_likelihood(
    f::Float64,
    r::Vector{Float64},
    i_index::Vector{Int64},
    a::Vector{Float64},
    b::Vector{Float64},
    theta::Float64;
    logL = false,
)
    n_items = size(a, 1)
    likel = 0
    for i in i_index
        p = b[i] + a[i] * theta
        likel += r[i] * p - _log1p_c(_exp_c(p))
    end
    return likel::Float64
end

"""
	resp_gen(f, pars, theta, design)

# Description

Generate dichotomous responses from a dataframe of item parameters and a vector of ability points (theta).

# Arguments

- **`pars::DataFrames.DataFrame`** : Required. DataFrame with item parameters (difficulty: b or d, discrimination: a or a1, guessing: c or g).
- **`theta::Vector{Float64}`** : Required. Vector of ability points.
- **`design::Matrix{Float64}`** : Optional. Default: `"at-b`". Values:  `"at-b`", `"at-ab`", `"at+b`", `"at+ab`". IRT model parametrization. Ex: at-b is ``Da(\\theta-b)``.

# Output

- **`i::Matrix{Float64}`**: Matrix of the item information. 
"""
function resp_gen(
    pars::DataFrames.DataFrame,
    theta::Vector{Float64},
    design::Matrix{Float64};
    model = "2PL", #1PL, 2PL, 3PL
    parametrization = "at-b", #"at-b, at-ab, at+b, at+ab"
    D = 1,#any number (1, 1.6)
    method = "classic uniform",
)

    a, a2, b, c = _desume_pars(pars; model = model, parametrization = parametrization)

    n_items = size(b, 1)
    N = size(theta, 1)
    n_index = Array{Array{Int64,1},1}(undef, n_items)
    i_index = Array{Array{Int64,1},1}(undef, N)
    for n = 1:N
        i_index[n] = findall(design[:, n] .== 1.0)
        if n <= n_items
            n_index[n] = findall(design[n, :] .== 1.0)
        end
    end
    p = Matrix(undef, n_items, N)
    lp = Matrix(undef, n_items, N)
    lq = Matrix(undef, n_items, N)
    for i = 1:n_items
        for n = 1:N
            pr = c[i] + ((1 - c[i]) * (1 / (1 + exp(-(b[i] + a[i] * theta[n])))))
            p[i, n] = pr
            lp[i, n] = log(pr)
            lq[i, n] = log(1 - pr)
        end
    end
    resp = Matrix(undef, n_items, N)
    if method == "cumulated students"
        @fastmath @inbounds for i = 1:n_items
            gapScore = 3
            while gapScore >= 2
                p2 = p[i, :]
                p2 = hcat((1 .- p2), p2)#2
                unif = rand(Uniform(0, 1), N)#5
                n = 1#6
                while n <= N#7
                    csum = p2[n, 1]#8
                    cat = 0#9
                    if design[i, n] == 0
                        resp[i, n] = 0
                        n = n + 1
                    else
                        while csum < unif[n]
                            cat = cat + 1
                            if (cat == 2)
                                break
                            end
                            csum = csum + p2[n, (cat+1)]
                        end
                        resp[i, n] = cat
                        n = n + 1
                    end
                end
                gapScore = abs(sum(resp[i, :]) - sum(p[i, n] for n in n_index[i]))
            end
        end
    end
    if method == "cumulated items"
        @fastmath @inbounds for n = 1:N#4
            #gapScore=3
            #while gapScore>=2
            p2 = p[:, n]
            p2 = hcat((1 .- p2), p2)#2
            unif = rand(Uniform(0, 1), n_items)#5
            StatsBase.samplei =
                StatsBase.sample(collect(1:n_items), n_items, replace = false)
            i = 1#6
            while i <= n_items#7
                csum = p2[StatsBase.samplei[i], 1]#8
                cat = 0#9
                if design[StatsBase.samplei[i], n] == 0
                    resp[StatsBase.samplei[i], n] = 0#missing
                    i = i + 1
                else
                    while csum < unif[StatsBase.samplei[i]]
                        cat = cat + 1
                        if (cat == 2)
                            break
                        end
                        csum = csum + p2[StatsBase.samplei[i], cat+1]
                    end
                    resp[StatsBase.samplei[i], n] = cat
                    i = i + 1
                end
            end
            #gapScore=abs(sum(skipmissing(resp[:,n]))-sum(p[i,n] for i in i_index[n]))
            #if gapScore>=0.5
            #	println("person ",n," gap=",gapScore)
            #end
            #end
        end
    end
    if method == "classic uniform"
        for n = 1:N#4
            #gapScore=2
            #while gapScore>=0.5
            unif = rand(Uniform(0, 1), n_items)#5
            StatsBase.samplei =
                StatsBase.sample(collect(1:n_items), n_items, replace = false)
            i = 1#6
            while i <= n_items#7
                if design[StatsBase.samplei[i], n] == 0
                    resp[StatsBase.samplei[i], n] = 0#missing
                    i = i + 1
                else
                    if unif[StatsBase.samplei[i]] < p[StatsBase.samplei[i], n]
                        resp[StatsBase.samplei[i], n] = 1
                    else
                        resp[StatsBase.samplei[i], n] = 0
                    end
                    i = i + 1
                end
            end
            #gapScore=abs(sum(skipmissing(resp[:,n]))-sum(p[i,n] for i in i_index[n]))
            #if gapScore>=0.5
            #	println("person ",n," gap=",gapScore)
            #end
            #end
        end
    end
    if method == "cumulated pattern"
        for n = 1:N#4+
            println("start person ", n)
            i_index = i_index[n]
            i_index = size(i_index, 1)
            p2 = p[i_index, n]
            p2 = hcat((1 .- p2), p2)#2
            patterns = Vector{Vector{Int64}}(undef, 1)
            nPattern = Vector{Int64}(undef, 1)
            for r = 1:1000
                respn = Vector{Int64}(undef, i_index)
                unif = rand(Uniform(0, 1), i_index)#5
                StatsBase.samplei =
                    StatsBase.sample(collect(1:i_index), i_index, replace = false)
                for i = 1:i_index
                    csum = p2[StatsBase.samplei[i], 1]#8
                    cat = 0#9
                    while csum < unif[StatsBase.samplei[i]]
                        cat = cat + 1
                        if (cat == 2)
                            break
                        end
                        csum = csum + p2[StatsBase.samplei[i], cat+1]
                    end
                    respn[i] = cat
                end
                if r == 1
                    patterns[1] = respn
                    nPattern[1] = 1
                else
                    println(size(patterns, 1))
                    corr = 0
                    for pat = 1:size(patterns, 1)
                        if patterns[pat] == respn
                            corr = pat
                        end
                    end
                    if corr == 0
                        push!(patterns, respn)
                        nPattern = hcat(nPattern, 1)
                    else
                        nPattern[corr] = nPattern[corr] + 1
                    end
                end

            end
            nPatterns = size(patterns, 1)
            probPatterns = nPattern ./ 1000
            println(probPatterns)
            resp[i_index, n] .= StatsBase.sample(patterns, pweights(probPatterns), 1)
            println("end person ", n)
        end
    end
    if method == "classic uniform pattern"
        for n = 1:N#4+
            println("start person ", n)
            i_index = i_index[n]
            i_index = size(i_index, 1)
            patterns = Vector{Vector{Int64}}(undef, 1)
            nPattern = Vector{Int64}(undef, 1)
            for r = 1:1000
                respn = Vector{Int64}(undef, i_index)
                unif = rand(Uniform(0, 1), i_index)#5
                StatsBase.samplei =
                    StatsBase.sample(collect(1:i_index), i_index, replace = false)
                for i = 1:i_index#7
                    if unif[StatsBase.samplei[i]] < p[StatsBase.samplei[i], n]
                        respn[i] = 1
                    else
                        respn[i] = 0
                    end
                end
                if r == 1
                    patterns[1] = respn
                    nPattern[1] = 1
                else
                    println(size(patterns, 1))
                    corr = 0
                    for pat = 1:size(patterns, 1)
                        if patterns[pat] == respn
                            corr = pat
                        end
                    end
                    if corr == 0
                        push!(patterns, respn)
                        nPattern = hcat(nPattern, 1)
                    else
                        nPattern[corr] = nPattern[corr] + 1
                    end
                end
            end
            nPatterns = size(patterns, 1)
            probPatterns = nPattern ./ 1000
            resp[i_index, n] .= StatsBase.sample(patterns, pweights(probPatterns), 1)
            println("end person ", n)
        end
    end
    return resp
end
