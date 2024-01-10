#include "tex_end.h"
#include <stdio.h>

void tex_end(char *filename)
{
  FILE *fi1 = fopen(filename,"a");
  fprintf(fi1,"\\end{document}\n");
  fclose(fi1);

}

