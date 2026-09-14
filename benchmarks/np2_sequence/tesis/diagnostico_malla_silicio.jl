# =============================================================================
#  Diagnostico: convergencia con la malla del SCF 3P del silicio
#
#  Frente a Froese Fischer, el silicio tiene las diferencias mayores de la secuencia en <r^-3>
#  (0.34 %), en T (4.4e-3 Ha) y en el cociente virial (1.5e-5), mientras que su energia y sus F^k
#  coinciden a ~1e-6 Ha y ~1e-5. Este diagnostico prueba si la causa es la malla: repite el SCF de
#  silicon_rohf.jl con N_elems = 100, 200 y 300. El script se copia en memoria cambiando solo
#  N_elems; no se guarda ningun .jld2 y silicon_rohf.jl no se toca.
#
#  Resultado (2026-09-13): convergido desde N = 200, y <r^-3> se aleja 2e-5 de Froese Fischer. No
#  es la malla (docs/claude/REVISION-RESULTADOS-2026-09-13.md, seccion 8).
#
#  Uso: julia --project=. benchmarks/np2_sequence/tesis/diagnostico_malla_silicio.jl   (~30 s)
# =============================================================================

isdefined(Main, :CABECERA_CI) || include(joinpath(@__DIR__, "comun.jl"))
isdefined(Main, :DIISState) || include(joinpath(NP2, "rohf_diis.jl"))

function main_malla_silicio(Ns = (100, 200, 300))
    src = read(joinpath(NP2, "silicon_rohf.jl"), String)
    occursin("N_elems = 100", src) || error("silicon_rohf.jl ya no tiene 'N_elems = 100'.")
    i = INFO["Si"]
    ff = FF["Si"]
    @printf("\n%-8s %16s %14s %14s %10s %12s %12s\n", "malla", "E (Ha)", "T (Ha)", "virial", "<r^-3>", "F0", "F2")
    for N in Ns
        code = replace(src, "N_elems = 100" => "N_elems = $N",
                       "function solve_silicon_rohf(" => "function solve_silicon_rohf_N$(N)(")
        code = replace(code, r"Pkg\.activate\([^\n]*\n" => "")
        code = replace(code, r"if abspath\(PROGRAM_FILE\) == @__FILE__.*"s => "")
        include_string(Main, code, "silicon_rohf_N$(N).jl")
        r = Base.invokelatest(getfield(Main, Symbol("solve_silicon_rohf_N$(N)")), R_MAX;
                              estado = "3P", save = false, verbose = false, use_diis = false)
        ws = cached_init_scf_workspace(R_MAX, N, Val(i.K), i.Z; γ = i.gamma, calc_R_matrices = true)
        orbs = r.orbitals
        T = sum(o.occ * dot(o.coeffs, (ws.T .+ (o.l * (o.l + 1) / 2) .* ws.R_inv2) * o.coeffs)
                for o in orbs if o.occ > 0)
        v = orbs[end]
        @printf("N = %-4d %16.8f %14.6f %14.9f %10.6f %12.8f %12.8f\n", N, r.E_total, T, -(r.E_total - T) / T,
                dot(v.coeffs, ws.R_inv3 * v.coeffs), compute_Rk(ws, v, v, v, v, 0), compute_Rk(ws, v, v, v, v, 2))
    end
    @printf("%-8s %16s %14s %14s %10s %12s %12s\n", "F.F.", ff.E, ff.T, ff.virial, ff.r3, ff.F0, ff.F2)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_malla_silicio()
end
