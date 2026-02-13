# For hydrogenic atoms
function solve_orbital!(orb::Orbital, atom::Atom, basis)
    H, S = assemble_hamiltonian(basis, atom.Z, orb.l)
    _solve_and_update!(orb, H, S)
end

# For SCF (takes an interaction matrix V_ext)
function solve_orbital!(orb::Orbital, atom::Atom, basis, V_ext::AbstractMatrix)
    H_core, S = assemble_hamiltonian(basis, atom.Z, orb.l)
    H_total = H_core + V_ext # we sum J (in the future K)
    _solve_and_update!(orb, H_total, S)
end

# Auxiliar function
function _solve_and_update!(orb, H, S)
    state_idx = orb.n - orb.l
    evals, evecs = solve_eigen(H, S, nev=state_idx)
    orb.energy = evals[state_idx]
    orb.coeffs = evecs[:, state_idx]
end


