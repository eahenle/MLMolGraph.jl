COFs from [Mercado](https://archive.materialscloud.org/record/2018.0003/v3) were reduced to their primitive cells
using pymatgen.

The primitive cells were loaded into Xtals and bonds were inferred; crystals with errors during
loading or bond inference were removed from the set (filter_and_bond_cif2jld.jl;
Mercado_COFs/bond_issues) and the remaining structures were saved as JLDs 
(filter_and_bond_cif2jld.jl; Mercado_COFs/bonded_xtal_jlds).


