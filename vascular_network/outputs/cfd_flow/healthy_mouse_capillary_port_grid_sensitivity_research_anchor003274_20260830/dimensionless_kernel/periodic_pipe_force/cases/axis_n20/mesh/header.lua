 version = 1
 label = 'fluidTree'
 comment = 'axis_n20; wall_libb continuous q; periodic axial direction'
 boundingbox = {
    origin = {
           0.000000000000000E+00,
           0.000000000000000E+00,
           0.000000000000000E+00 
    },
    length =   40.388853197829045E-06 
}
 nElems = 37920
 minLevel = 8
 maxLevel = 8
 nProperties = 2
 effBoundingbox = {
    origin = {
          18.616737020874326E-06,
          18.616737020874326E-06,
          10.728289130673340E-06 
    },
    effLength = {
           3.155379156080392E-06,
           3.155379156080392E-06,
          18.932274936482361E-06 
    } 
}
 property = {
    {
        label = 'has boundaries',
        bitpos = 3,
        nElems = 9600 
    },
    {
        label = 'has qVal',
        bitpos = 8,
        nElems = 9120 
    } 
}
