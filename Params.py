# Params.py
# control all global variables

## materials definition. (atomic num, atomic weight, density (g/cm3), radiation length (m))
materials = { "fe" : (26, 55.845, 7.874, .01757), 
              "si" : (14, 28.0855, 2.329, .0937),
              "air": (7.34, 14.719, 1.205e-3, 3.04e2),  
              "c"  : (6, 12.0107, 2.0, .2135)  }  

Q = 1  ## in units of e
m = 105.658  ## in MeV 
solRad = 7.3  ## in m
solLength = 21.6   ## in m
BMag = 4.0 ## in T
MSCtype = 'kuhn'
