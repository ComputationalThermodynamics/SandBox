# General data (Igneous database, Holland et al., 2018)
struct global_variable{N,_T, _I, _C}
     oxides      ::  NTuple{N,NTuple{5, _C}}  
     apo         ::  SVector{N,_T};
     R           ::  _I
end

function init_sys_infos()

    oxides = ("SiO2 ", "Al2O3", "CaO  ", "MgO  ", "FeO  ", "K2O  ", "Na2O ", "TiO2 ", "O    ", "Cr2O3", "H2O  ")
    apo    = SVector{11}([3., 5., 2., 2., 2., 3., 3., 3., 1., 5., 3.]);
    R      = 0.0083144;                            # gas constant

    gv =  global_variable(      str2char.(oxides),
                                apo,
                                R       )

    return gv
end

