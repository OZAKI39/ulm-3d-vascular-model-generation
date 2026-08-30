 version = 1
 label = 'fluidTree'
 comment = 'axis_n16; wall_libb continuous q; periodic axial direction'
 boundingbox = {
    origin = {
           0.000000000000000E+00,
           0.000000000000000E+00,
           0.000000000000000E+00 
    },
    length =   50.486066497286306E-06 
}
 nElems = 19968
 minLevel = 8
 maxLevel = 8
 nProperties = 2
 effBoundingbox = {
    origin = {
          23.665343670602957E-06,
          23.665343670602957E-06,
          15.776895780401969E-06 
    },
    effLength = {
           3.155379156080396E-06,
           3.155379156080396E-06,
          18.932274936482361E-06 
    } 
}
 property = {
    {
        label = 'has boundaries',
        bitpos = 3,
        nElems = 6056 
    },
    {
        label = 'has qVal',
        bitpos = 8,
        nElems = 5760 
    } 
}
