using LinearAlgebra: Diagonal, diag, diagm
using Random: randn
using Statistics: std
using WAV

function genFullyDiagMs(K::Int, size::Int)

    C0 = Diagonal(ones(size))
    Cdiag = [Diagonal(randn(size)) for _ in 1:K]
    Cs = pushfirst!(Cdiag, C0)

    A = randn(size, size)

    Cx = Array{Float64,3}(undef, size, size, K + 1)
    for i in 1:K+1
        Cx[:, :, i] = A * Cs[i] * A'
    end
    return Cx, A
end


function genApproxDiagMs(K::Int, size::Int)
    C0 = Diagonal(ones(size))

    Cdiag = [Diagonal(randn(size)) for _ in 1:K]
    Cs = pushfirst!(Cdiag, C0)

    A = randn(size, size)
    Ax = [A]

    Cx = Array{Float64,3}(undef, size, size, K + 1)
 
    for i in 1:K+1
        Ak = A .+ randn(size, size) .* (10^(-2))
        Cx[:, :, i] = Ak * Cs[i] * Ak'
        push!(Ax, Ak)
    end 
    return Cx, A  
end

function genTimeCorrMs(X, lags, method = "cov", symmetrize=false)
    K = length(lags)

    T, N = size(X)

    if T < N
        @info "X may be oriented the wrong way" 
    end

    C = zeros(N, N, K)
    C0 = zeros(N, N)

    for k in 1:K
        if lags[k] >= 0
            x1 = X[1:end-lags[k], :]
            x2 = X[1+lags[k]:end, :]
        else
            x1 = X[1-lags[k]:end, :]
            x2 = X[1:end+lags[k], :]
        end

        T = size(x1, 1)

        if method == "cov"
            x1 = x1 .- repeat(sum(x1, dims=1) / T, T, 1)
            x2 = x2 .- repeat(sum(x2, dims=1) / T, T, 1)
            C[:, :, k] = x1' * x2 / (T-1)
        elseif method == "moments"
            C[:, :, k] = x1' * x2 / T
        else
            @error "Method not recognized"
        end

        if symmetrize
            C[:, :, k] = (C[:, :, k] + C[:, :, k]') / 2
        end

        if lags[k] == 0
            C0 = C[:, :, k]
        end
    end

    return C, C0
end

function generate_source(N::Int64, Ns::Int64, sig_type::String)

    if sig_type == "speech"

        s1, Fs = wavread("test/soundfiles/sentence_male_29s.wav")
        s2, Fs = wavread("test/soundfiles/sentence_female_28s.wav")
        s3, Fs = wavread("test/soundfiles/henry_theater_male_30s.wav")
        s4, Fs = wavread("test/soundfiles/numbers_female_29s.wav")
        s5, Fs = wavread("test/soundfiles/Iam_female_30s.wav")
        s6, Fs = wavread("test/soundfiles/Music_trumpet_30s.wav")
        s7, Fs = wavread("test/soundfiles/Music_piano_30s.wav")

        S = [s1[1:Ns], s2[1:Ns], s3[1:Ns], s4[1:Ns], s5[1:Ns], s6[1:Ns], s7[1:Ns]]
        S = S[:, 1:end]  # keep only N sources

    elseif sig_type == "supergauss"
        S = randn(Ns, N) .^ 5

    elseif sig_type == "sine"
        t = 0:16*pi/Ns:16*pi
        t[end] = []
        w = rand(1, N)
        S = zeros(Ns, N)
        for n in 1:N
            S[:, n] = sin(w[n] * t)
        end

    else
        error("Invalid sig_type")
    end
    return S
end

