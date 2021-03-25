
"""
    resp_gen(
        pars::DataFrames.DataFrame,
        theta::Vector{Float64},
        design::Matrix{Float64};
        model = "2PL", #1PL, 2PL, 3PL
        parametrization = "at-b", #"at-b, at-ab, at+b, at+ab"
        D = 1,#any number 
        method = "classic uniform",
    )

# Description

Generate dichotomous responses from a dataframe of item parameters and a vector of ability points (theta).

# Arguments

- **`pars::DataFrames.DataFrame`** : Required. DataFrame with item parameters (difficulty: b or d, discrimination: a or a1, guessing: c or g).
- **`theta::Vector{Float64}`** : Required. Vector of ability points.
- **`design::Matrix{Float64}`** : Required. The design matrix.
- **`model`**: Optional. Default
- **`parametrization`** : Optional. Default: \"at-b\". Values:  \"at-b\", \"at-ab\", \"at+b\", \"at+ab\". IRT model parametrization. Ex: at-b is ``Da(\\theta-b)``.
- **`D`**: Optional. Default: 1.
- **`method`**: Optional. Default: \"classic uniform\".

# Output

- **`i::Matrix{Float64}`**: Matrix of binary responses. 
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

function sum_pos(x::Vector{Float64})
x = x[x.>0]
if size(x, 1) == 0
    return zero(Float64)::Float64
else
    return sum(x)::Float64
end
end

function read_settings_file(settings_file::String)
include(settings_file)
return Inputs::InputSettings
end
