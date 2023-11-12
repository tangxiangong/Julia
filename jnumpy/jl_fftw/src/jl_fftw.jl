module jl_fftw

using TyPython
using TyPython.CPython
using MKL
import FFTW

const FFT1D_PlanType = Union{Base.return_types(FFTW.plan_fft, (StridedVector{ComplexF64}, Int))...}
const plan1Ds = Dict{Tuple{DataType, NTuple{N, Int} where N}, FFT1D_PlanType}()

"""
1-D FFT
"""
@export_py function fft(x::StridedVector)::Vector{ComplexF64}
    x = collect(ComplexF64, x)
    plan = get!(plan1Ds, (typeof(x), size(x))) do
        FFTW.plan_fft(copy(x), 1; flags=FFTW.MEASURE)
    end
    return plan * x
end

@export_py function mc_pi(n::Integer)::Float64
    count = 0 
    for _ in 1:n
        x = rand()
        y = rand()
        if x^2 + y^2 < 1
            count += 1
        end
    end
    return 4*count/n
end

function init()
    @export_pymodule _fast_fft begin
        fft = Pyfunc(fft)
        mc_pi = Pyfunc(mc_pi)
    end
end

precompile(init, ())

end # module jl_fftw
