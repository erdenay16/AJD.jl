using WAV: wavread
using Random: randn
using AJD: ffdiag, QDiag

function load_audio(filename, samples)
    y, _ = wavread(filename)
    y = y[1:samples]
    return y
end

function genSrc()

    samples = 10000
    dir = @__DIR__

    s1 = load_audio(joinpath(dir, "sentence_male_29s.wav"), samples)
    s2 = load_audio(joinpath(dir, "sentence_female_28s.wav"), samples)
    s3 = load_audio(joinpath(dir, "henry_theater_male_30s.wav"), samples)

    S = hcat(s1, s2, s3)
    return S
end

function genMix(S, lags)
    K = length(lags)
    T, N = size(S)

    if T < N
        @info "S may be oriented the wrong way"
    end

    C = zeros(N, N, K)

    for k in 1:K

        if lags[k] > 0
            x1 = S[1:end-lags[k], :]
            x2 = S[1+lags[k]:end, :]
        else
            x1 = S[(1-lags[k]):end, :]
            x2 = S[1:end+lags[k], :]
        end

        T, _ = size(x1)
        
        # moments
        C[:, :, k] = x1' * x2 / T 

        # symmetrize is true
        C[:, :, k] = (C[:, :, k] + C[:, :, k]') / 2

    end
    
    return C
end

function timeDelayedCorrMs()

    S = genSrc()
    lags = 0:5:100

    C = genMix(S, lags)

    N = size(S, 2)
    A = randn(N, N)
    
    I = size(C, 3)

    for i in 1:I
        C[:,:,i] = A * C[:,:,i] * A'
    end

    return C, C[:,:,1], A
end