# Add a Membrane

With the optimized protein (`PaSDH_opt_H.pqr`) and the equilibrated Memembrane snapshot (bomdrr.sav.box75), we can add a membrane:



`python build_membrane.py`



# Add a Water

With the equilibrated water system (902 water molecules in 30 A^3), we add water molecules to the system consisting of PaSDH + membrane (DPPE).

`python build_water.py`



__
