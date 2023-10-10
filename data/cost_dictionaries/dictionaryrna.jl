# -------------------
# Nucleic Acids Cost
# -------------------
# A Adenine  : C5 H5 N5 O0
# C Cytosine : C4 H5 N3 O1
# G Guanine  : C5 H5 N5 O1
# T Thymine  : C5 H6 N2 O2
# U Uracil   : C4 H4 N2 O2

dict_RNA_C = Dict(
    RNA_A => 5,
    RNA_C => 4,
    RNA_G => 5,
    RNA_U => 4
)

dict_RNA_H = Dict(
    RNA_A => 5,
    RNA_C => 5,
    RNA_G => 5,
    RNA_U => 4
)

dict_RNA_N = Dict(
    RNA_A => 5,
    RNA_C => 3,
    RNA_G => 5,
    RNA_U => 2
)

dict_RNA_O = Dict(
    RNA_A => 0,
    RNA_C => 1,
    RNA_G => 1,
    RNA_U => 2
)