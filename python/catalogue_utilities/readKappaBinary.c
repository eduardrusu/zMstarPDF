#include <stdio.h>
#include <stdlib.h>

int main (int argc, char** argv) {
  
  FILE *fpt;
  float* kappa;
  
  const size_t number_of_elements = 4096 * 4096;
  const size_t size_of_elements   = sizeof(float);
  const char kappa_file_name[]   = "kappa_maps_z_1.3857/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.kappa";
  
  kappa = calloc(number_of_elements, size_of_elements);
  
  if (!kappa)
    {
      fprintf(stderr, "Cannot allocate memory for kappa array\n");
      exit(-1);
    }

  //fprintf(stdout, "%s\n", kappa_file_name);
	 
  fpt = fopen (kappa_file_name, "r");
  
  if (!fpt)
    {
      fprintf(stderr, "Cannot open kappa file\n");
      exit(-1);
    }
	 
  fread(kappa, size_of_elements, number_of_elements, fpt);
  fclose (fpt);


  fpt = fopen ("kappa_values.dat", "w");
  int i;
  for (i=0; i<number_of_elements; i++)
    {
      fprintf(fpt, "%d\t%lf\n", i, kappa[i]);
    }
  fclose(fpt);

  fprintf(stdout, "test kappa value:\n %f  %f  %f\n", kappa[0], kappa[1000+2000*4096], kappa[4095+4095*4096]);

  free(kappa);
  return 0;
  
}
