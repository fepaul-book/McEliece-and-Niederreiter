using DelimitedFiles

using .Niederreiter

println("start")
println("processing")
a = generateKeysNiederreiter(16, 3, 4)
println(a)
println("finished")