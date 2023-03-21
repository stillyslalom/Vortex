# Postprocessing algorithms for PIV data

## Conversion of Prana data to Julia data structures


## Spurious vector detection


abstract type VectorMedianMeasure end
const VMM = VectorMedianMeasure

struct MedianMagnitude <: VectorMedianMeasure end

struct MedianComponents <: VectorMedianMeasure end

abstract type SpuriousVectorDetector end

struct MedianFilter{MM <: VMM} <: SpuriousVectorDetector
    δ::Int
    ϵ::Float64
    mm::MM
end

function (f::MedianFilter)(u, v)

end

struct NormalizedMedianFilter{MM <: VMM} <: SpuriousVectorDetector
    δ::Int
    ϵ::Float64
    ϵ₀::Float64
    mm::MM
end

struct NeighborDifference <: SpuriousVectorDetector
    δ::Int
    ϵ::Float64
end

struct GlobalHistogram <: SpuriousVectorDetector
    N::Int # number of clusters
end

struct DynamicMean <: SpuriousVectorDetector
    δ::Int
    C₁::Float64
    C₂::Float64
end
