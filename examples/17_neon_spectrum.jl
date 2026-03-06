using JLD2
using LinearAlgebra
using Printf

function compute_ns_2p_spectrum(data_file::String)
    # 1. Load the data
    data = load(data_file)
    
    R_mat = data["R_mat"]
    evals_fs = data["evals_fs"]
    evecs_fs = data["evecs_fs"]
    evals_fp = data["evals_fp"]
    evecs_fp = data["evecs_fp"]
    active_s = data["active_s"]
    active_p = data["active_p"]
    n = data["num_splines"]

    # 2. Isolate the occupied 2p state
    c_2p = zeros(Float64, n)
    c_2p[active_p] .= evecs_fp[:, 1]
    E_2p = evals_fp[1]

    println("\n=== Synthetic Emission Spectrum: ns -> 2p ===")
    @printf("%-8s | %-12s | %-15s | %-15s\n", "Upper", "Energy (Ha)", "Wavelength (nm)", "Radial Dipole Sq")
    println("-"^58)

    # 3. Loop through virtual s-states (3s, 4s, 5s...)
    for i in 3:10
        E_ns = evals_fs[i]
        
        # Stop if we hit the continuum (positive energies)
        # if E_ns >= 0.0
        #     break 
        # end
        
        c_ns = zeros(Float64, n)
        c_ns[active_s] .= evecs_fs[:, i]
        
        dE_Ha = E_ns - E_2p
        dE_eV = dE_Ha * 27.211386
        wavelength_nm = 1239.84193 / dE_eV
        
        # Calculate the Radial Dipole Matrix Element
        radial_dipole = dot(c_ns, R_mat * c_2p)
        dipole_sq = radial_dipole^2
        
        state_label = "$(i)s"
        @printf("%-8s | %12.6f | %15.2f | %15.6e\n", state_label, E_ns, wavelength_nm, dipole_sq)
    end
end

compute_ns_2p_spectrum("sodium_core_results.jld2")
