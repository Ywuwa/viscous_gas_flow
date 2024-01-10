To build the code you need to create makefile by yourself or run it in an IDE



DESCRIPTION FILE:

hcond_vgas_2.dat 



PROGRAM FILES:

otl_glad.cpp - main file

input.cpp - read hcond_vgas_2.dat and init structures

setka.cpp - init not rectangular domain

hcond_vgas_2.cpp - compute functions

tex_init.cpp

tex_tab.cpp

tex_end.cpp - creates tex table with the results



ADDITIONAL:

cartoon_warmmap_not_rectangular_area.py - visualize gas dynamic (the visualized function is set up inside hcond_vgas_2.cpp in print section)
