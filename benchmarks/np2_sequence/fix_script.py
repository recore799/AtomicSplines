import re

with open('germanium_rohf_vpol.jl', 'r') as f:
    content = f.read()

# Fix the broken sed replacements
content = re.sub(r'# estado = strip\(readline\(\)\)\n\s*# if estado \["av", "3P"\]\n\s*# println\. Usando \'av\' por defecto\."\)\n\s*# estado = "av"\n\s*# end', '', content)
content = re.sub(r'filename = "germanium_rohf_results_\$\(estado\)_R\$\(R_max\)\.jld2"', 'filename = "germanium_rohf_results_$(estado)_R$(R_max)_ad$(@sprintf("%.3f", alpha_d)).jld2"', content)
content = content.replace('# end', 'end')
content = content.replace('if estado ["av", "3P"]', '')
content = content.replace('println. Usando \'av\' por defecto.")', '')
content = content.replace('estado = "av"', '')
content = content.replace('solve_germanium_rohf(30.0)', '')
content = content.replace('ws = cached_init_scf_workspace(R_max, N_elems, Val(8), Z; γ=3.0, calc_R_matrices=true)', 'ws = cached_init_scf_workspace(R_max, N_elems, Val(8), Z; γ=3.0, calc_R_matrices=true, alpha_d=alpha_d, r_c=r_c)')

with open('germanium_rohf_vpol.jl', 'w') as f:
    f.write(content)
