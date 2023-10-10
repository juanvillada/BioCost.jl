# -------------------
# Nucleic Acids Cost
# -------------------
# A Adenine  : C5 H5 N5 O0
# C Cytosine : C4 H5 N3 O1
# G Guanine  : C5 H5 N5 O1
# T Thymine  : C5 H6 N2 O2
# U Uracil   : C4 H4 N2 O2

dict_DNA_C = Dict(
    DNA_A => 5,
    DNA_C => 4,
    DNA_G => 5,
    DNA_T => 5
)

dict_DNA_H = Dict(
    DNA_A => 5,
    DNA_C => 5,
    DNA_G => 5,
    DNA_T => 6
)

dict_DNA_N = Dict(
    DNA_A => 5,
    DNA_C => 3,
    DNA_G => 5,
    DNA_T => 2
)

dict_DNA_O = Dict(
    DNA_A => 0,
    DNA_C => 1,
    DNA_G => 1,
    DNA_T => 2
)