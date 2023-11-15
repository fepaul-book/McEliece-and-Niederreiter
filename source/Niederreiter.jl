#=
Functions for the Niederreiter cryptosystem using Goppa codes.

Author: Felix Peter Paul
Last updated: 2023/November/15

IMPLEMENTS NIEDERREITER ONLY AS A KEM (CAN ONLY ENCRYPT RANDOM ERROR VECTOR OF WEIGHT EQUAL t 
=#

module Niederreiter

export
    generateKeysNiederreiter,
    saveKeys,
    encrypt,
    decrypt,
    haltingCheck

include("./Goppa.jl")
using .GoppaCode

include("./FiniteFieldMatrix.jl")
using .FiniteFieldMatrix

include("./Misc.jl")
using .Misc

using Random
using StatsBase
using DelimitedFiles

"""
    haltingCheck(publicKey, privateKey, samples=2000)

Checks that a given pair of keys correctly finishes, i.e. always determines whether or not a cyphertext
is valid or not, by trying to decrypt random vectors.
"""
function haltingCheck(publicKey, privateKey, samples=2000)
    t = publicKey[2]
    m = privateKey[5]
    k = size(publicKey[1], 1)
    n = size(publicKey[1], 2)

    failures = 0
    for i in 1:samples
        test = [Array{Bool, 1}([rand([0, 1]) for _ in 1:n])]
        try
            decrypt(test, publicKey[1], privateKey[1], privateKey[2], privateKey[3], privateKey[4], m, t)
        catch exc
            if isa(exc, InvalidCodewordError)
                failures += 1
                if i % 100 == 0
                    println("Already $(i) steps.")
                end
                continue
            else
                throw(exc)
            end
        end
    end

    return failures
end


"""
generateKeysNiederreiter(n, t, m)

Generates the public and private keys for a McEliece cryptosystem with length n
and with capacity to correct up to t errors.
returns (GPub, t), (GSInverse, PInverse, goppaPolynomial, codeSupport, m)
"""
function generateKeysNiederreiter(n, t, m)
    if n > 2^m
        error("Extension field can't support $(n) elements.")
    end
    if n - m*t < 1
        error("Incompatible parameters: $(n) - $(m)*$(t) < 0.")
    end

    G = GInverse = goppaPolynomial = codeSupport = k = [] # Just to make sure these variables exist in this scope.
    keysCreated = false
    while !keysCreated
        try
            # Generating keys until we have a G with full rank.
            G, goppaPolynomial, codeSupport, k = generateGoppaCodeNiederreiter(n, t, m)
            GInverse = ControlMatrixInverse!(G)
            keysCreated = true
        # To save time, it always tries to invert G right away. If G is not invertible,
        # the function throws an exception that gets catched and only indicates that another
        # goppa code is needed.
        catch exc
            # throw(exc)
            continue
        end
    end

    P = getPermutation(n)
    S, SInverse = getScrambler(k)

    GPub = (S * G) * P
    GPub = GPub .% 2 # S * GP can break out of GF(2)

    PInverse = Array{UInt32, 2}(inv(P))
    GSInverse = (GInverse * SInverse) .% 2

    return (GPub, t), (GSInverse, PInverse, goppaPolynomial, codeSupport, m)
end

"""
    saveKeys(publicKey, privateKey, publicPrefix="public", privatePrefix="private")

Saves the given keys to text files in a format friendly for Julia.
"""
function saveKeys(publicKey, privateKey, publicPrefix="public", privatePrefix="private")
    writedlm(publicPrefix * "_g.txt", publicKey[1])
    open(publicPrefix * "_t.txt", "w") do file
        write(file, string(publicKey[2]))
    end

    writedlm(privatePrefix * "_gs.txt", privateKey[1])
    writedlm(privatePrefix * "_p.txt", privateKey[2])
    writedlm(privatePrefix * "_polynomial.txt", privateKey[3])
    writedlm(privatePrefix * "_support.txt", privateKey[4])

    open(privatePrefix * "_m.txt", "w") do file
        write(file, string(privateKey[5]))
    end
end

"""
    getErrorVector(n, t)

Returns a random binary vector of n components with weigth t.
"""
function getErrorVector(n, t)
    positions = sample(1:n, t, replace=false)
    error = falses(n)

    for pos in positions
        error[pos] = true
    end

    return error
end

"""
    encrypt(plaintext, HPub, t)

Takes an arbitrarily long plaintext and encrypts with a public key from a Niederreiter
cryptosystem. Returns size(plaintext)/k codewords with t errors.
"""
function encrypt(HPub, t)
    n = size(HPub)[2]
    message = getErrorVector(n, t)
    cipher = (HPub * message) .%2
    return cipher
end

"""
    decryptWord!(cyphertext, H, HInverse, P, goppaPolynomial, codeSupport, m, t)

Decodes and decrypts a 1xn array of bits. It expects that P is the
invers of its respective matrix.
"""
function decryptWord!(cyphertext, H, HInverse, P, goppaPolynomial, codeSupport, m, t)
    codeword = findError(cyphertext, goppaPolynomial, codeSupport, m, t)
    plaintext = (P * permutedPlaintext) .% 2
    if !isPlaintextValid(plaintext, cyphertext, H, 0)
        error("The message is invalid.")
    end
    return plaintext
end

"""
    decrypt(msg, G, GInverse, P, goppaPolynomial, codeSupport, m, t)

Decrypts an array of codewords (msg) using a private key, and
returns the plaintext as a string.
"""
function decrypt(msg, H, HInverse, P, goppaPolynomial, codeSupport, m, t)
    plainbits = [decryptWord!(word, H, HInverse, P, goppaPolynomial, codeSupport, m, t) for word in msg]
    plainbits = [bit for bit in Iterators.flatten(plainbits)]
    return bits2string(plainbits)
end

end # Module Niederreiter