#include "tex_init.h"
#include <stdio.h>

void tex_init(char *filename)
{
  FILE *fi1 = fopen(filename,"a");
  fprintf(fi1,"\\documentstyle{article}\n");
  fprintf(fi1,"\\begin{document}\n");
  fclose(fi1);
}

