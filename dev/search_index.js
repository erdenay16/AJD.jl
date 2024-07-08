var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = AJD","category":"page"},{"location":"#AJD","page":"Home","title":"AJD","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for AJD.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [AJD]","category":"page"},{"location":"#AJD.QDiag","page":"Home","title":"AJD.QDiag","text":"QDiag(C_0, C, weights, approach, tolerance, maximum_iteration, \n    random_number_generator = Xoshiro()) -> AbstractArray{<:Real}\n\nThis is an implementation of the algorithm introduced in:      Vollgraf, Roland & Obermayer, Klaus. (2006).      Quadratic optimization for simultaneous matrix diagonalization.      Signal Processing, IEEE Transactions on. 54. 3270 - 3278. 10.1109/TSP.2006.877673. Namings of the parameters are consistent with the paper for easier understanding. \n\nArguments\n\nC_0::AbstractArray{<:Real}: This is the sphering matrix introduced as C^(0).\nC::AbstractArray{<:Real}: This is the set of matrices that is introduced as C\nweights::AbstractArray{<:Real}:: This is the weights vector mathbflpha. It will be   later normalized such that sumk=1K lpha_k = 1.\napproach::String: This is the flag for approaches \"NK3\" and \"N5\" introduced in the paper.\ntolerance::Real: This is the tolerance for the error.\nmaximum_iteration::Integer: Maximum number of iterations.\nrandomnumbergenerator::Xoshiro: This is random number generator for matrix W. Default   generator generates random numbers based on default_rng() but a seed can be introduced    by the user. This argument is added for testing purposes.\n\n\n\n\n\n","category":"function"}]
}