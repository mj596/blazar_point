#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

int
save(double col1[], double col2[], const char prefix[], double suffix, int n)
     /* saves vectors col1[1..n] and col2[1..n] in a file "prefix_suffix" */
{
  char *filename,*rad,*filename2;
  FILE *output,*output2;
  int i;
  
  filename=(char *)malloc((size_t)(23*sizeof(char)));
  filename2=(char *)malloc((size_t)(23*sizeof(char)));
  rad=(char *)malloc((size_t)(23*sizeof(char)));
  strcpy(filename,INPUT_FILE_ID);
  strcat(filename,"/");
  strcat(filename,prefix);
  strcat(filename,"_");
  sprintf(rad,"%.6f",suffix);
  strcat(filename,rad);
  strcpy(filename2,filename);
  strcat(filename2,"_normal");
 
  output = fopen(filename,"w");
  output2 = fopen(filename2,"w");

  for(i=1;i<=n;i++)
    {
      fprintf(output,"%e %e\n",(col1[i]<=0.0) ? 0.0 : log10(col1[i]), (col2[i]<=0.0) ? 0.0 :log10(col2[i]));
      fprintf(output2,"%e %e\n",(col1[i]<=0.0) ? 0.0 : col1[i], (col2[i]<=0.0) ? 0.0 :col2[i]);
    }

  fclose(output);
  fclose(output2);
  //free(rad);
  free(filename);
  //free(filename2);
  return 0;
}

int
save_normal(double col1[], double col2[], const char prefix[], int n)
     /* saves vectors col1[1..n] and col2[1..n] in a file "prefix" */
{
  char *filename;
  FILE *output;
  int i;
  
  filename=(char *)malloc((size_t)(23*sizeof(char)));
  strcpy(filename,INPUT_FILE_ID);
  strcat(filename,"/");
  strcat(filename,prefix);
  output = fopen(filename,"w");
  for(i=0;i<n;i++)
    fprintf(output,"%e %e\n",(col1[i]<=0.0) ? 0.0 : col1[i], (col2[i]<=0.0) ? 0.0 :col2[i]);
  fclose(output);
  free(filename);
  return 0;
}

int
save_0(double col1[], double col2[], const char prefix[], int n)
     /* saves vectors col1[1..n] and col2[1..n] in a file "prefix" */
{
  char *filename;
  FILE *output;
  int i;
  
  filename=(char *)malloc((size_t)(23*sizeof(char)));
  strcpy(filename,INPUT_FILE_ID);
  strcat(filename,"/");
  strcat(filename,prefix);
  output = fopen(filename,"w");
  for(i=0;i<n;i++)
    fprintf(output,"%e %e\n",(col1[i]<=0.0) ? 0.0 : log10(col1[i]), (col2[i]<=0.0) ? 0.0 :log10(col2[i]));
  fclose(output);
  free(filename);
  return 0;
}



int save_normal_radius(double col1[], double col2[], const char prefix[], double suffix, int n)
     /* saves vectors col1[1..n] and col2[1..n] in a file "prefix_suffix" */
{
  char *filename,*rad;
  FILE *output;
  int i;
  
  filename=(char *)malloc((size_t)(23*sizeof(char)));
  rad=(char *)malloc((size_t)(23*sizeof(char)));
  strcpy(filename,INPUT_FILE_ID);
  strcat(filename,"/");
  strcat(filename,prefix);
  strcat(filename,"_");
  sprintf(rad,"%.6f",suffix);
  strcat(filename,rad);
  output = fopen(filename,"w");
  for(i=0;i<n;i++)
    fprintf(output,"%e %e\n",(col1[i]<=0.0) ? 0.0 : col1[i], (col2[i]<=0.0) ? 0.0 :col2[i]);
  fclose(output);
  free(rad);
  free(filename);
  return 0;
}
