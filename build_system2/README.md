- align_two_proteins.py
  
  > This python code will translate and rotate a protein (i.e. A.pdb) to align with another protein (i.e. opm.pdb). The aligned protein (e.g. A_opm.pdb) will be saved.

        e.g. python align_two_proteins.py -t opm.pdb -f A.pdb -s A_opm.pdb 





- build_system.py
  
  > Add a membrane (dppe), water, and salt (NaCl). The length of the box must be a multiple of 15 A.

        e.g. python build_system.py -i input_build.json

```
{
    "fname_protein" : "A_opm.pdb",
    "protein_charge": -30,
    "fname_system" : "A_opm_system.pdb",
    "box" : [150, 150, 150],
    "membrane" : {
        "bool" : true,
        "fname_unit" : "dppe_box15.pdb"
    },
    "solvent" : {
        "bool" : true,
        "fname_unit" : "water_box15.xyz",
        "NaCl(Molarity)" : 0.1
    }
}
```


