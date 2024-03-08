module DensityMatrixPropagation

using LinearAlgebra, StaticArrays, SparseArrays
using Quantica
using DifferentialEquations
using FFTW

export densitymatrix, PeierlsHamiltonian, dimension, pdistance, ham, current, propagate
export MonocromaticSinField, SinSquaredPulseField, LorentzianPulseField, GaussianPulseField, StepFunctionField


include("peierls.jl")
include("densitymatrix.jl")
include("applyham.jl")
include("laserfields.jl")
include("timepropagation.jl")

end
