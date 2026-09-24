 version = 1
 label = 'fluidTree'
 comment = 'fixed physical Hagen-Poiseuille pipe, N=16'
 boundingbox = {
    origin = {
         -25.243033248643153E-06,
         -25.243033248643153E-06,
         -25.243033248643153E-06 
    },
    length =   50.486066497286306E-06 
}
 nElems = 19640
 minLevel = 8
 maxLevel = 8
 nProperties = 2
 effBoundingbox = {
    origin = {
          -1.577689578040197E-06,
          -1.577689578040197E-06,
          -9.466137468241182E-06 
    },
    effLength = {
           3.155379156080394E-06,
           3.155379156080394E-06,
          18.932274936482365E-06 
    } 
}
 property = {
    {
        label = 'has boundaries',
        bitpos = 3,
        nElems = 6024 
    },
    {
        label = 'has qVal',
        bitpos = 8,
        nElems = 5728 
    } 
}
