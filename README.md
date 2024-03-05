To build the code you need to create makefile by yourself or run it in an IDE
----
## DESCRIPTION FILE:

**hcond_vgas_2.dat** - if u gonna try to run the program, the file has to be in the same directory 

## PROGRAM FILES:

**otl_glad.cpp** - main file <br>
**input.cpp** - read hcond_vgas_2.dat and init structures <br>
**setka.cpp** - init not rectangular domain <br>
**hcond_vgas_2.cpp** - compute functions <br>


**tex_init.cpp** <br>
**tex_tab.cpp** <br>
**tex_end.cpp** - creates tex table with the results <br>
----

When program is finished there are many .txt files appeared in the directory, you can vizualize the result using:

## ADDITIONAL:

cartoon_warmmap_not_rectangular_area.py - visualize gas dynamic (the print into files is setting up inside hcond_vgas_2.cpp in print section)
