/*
  Program
    cal-xzsidesy20-wall

  Description
    Interactive program to approximate the xz deformation distance (wall at y=20) so as to match the density of a bulk system (no walls).
    Independant on the bond length, only scaling is applied to reach at a better prediction. Three simulation runs are mandatory to match
    the density of a bulk system (no walls): First run is an equilibrum run, with yhi = 20; second being a run with yhi = 20 and deforming x and z.
    Third run is with similar settings as the latter, and makes the approximation better.
    Scaling is independant of most material characterestics, so I guess its safe to work with this code at any chain length.

  Variables
    

  Keywords
    density, xzdeform, LAMMPS, polymers,

  Related
    slabxyt

  Author
    gg nov 29 2016
*/
#include <stdio.h>
#include <math.h> //gcc test.c -o test -lm
#include <stdlib.h>
int main() {
double sas, sa, sa2, m, y=20.0, va, v, paim, pr, pdr, xa1, xa2, xa3, pdr1, pdr2; //y set at 20
int choice, hit=0;
printf("[Read Description]\n");
printf("Enter nof particles (atoms):\n");
scanf("%lf", &m); //normally 16*384=6144
printf("Enter the density obtained from slabxyt (after fit [5:14]) for equilibration run, wall @ yhi=20:\n");
scanf("%lf", &pr); //scanning double requires %lf; equilibrum run had all sides at 20,20,20
printf("Enter the volume obtained from thermo data for bulk, no walls:\n");
scanf("%lf", &v); //scanning double requires %lf
//printf("Enter fixed position of y (yhi):\n");
//scanf("%lf", &y); //normally 20
/*calculation*/
double i; //search from midpoint
paim=m/v; //bulk density = aim
va=m/pr; //apparent volume
sas=va/y; //product of xz is sas
sa=sqrt(sas); //get x or z  : use -lm for compiling in gcc
/*searching for match*/
//double f=(m/(sa*sa*y)); //test-works
//printf("%.4lf\n", f); //test
Y:
printf("\nStatus: 1 or 2?\n 1 - After initial equilibration run with yhi=20\n 2 - Further scaling with xz deform, yhi=20 in comparison to equilib run\n");
scanf("%d", &choice);
if (choice==1)
{
  for(i=(y/2);i<=sa;i+=0.00001) 
  { //weird for loop
    if (((m/(i*i*y))-paim)<=0.00001) 
    {     //precision control
      sa=i;//update x or z with the find
      hit++;
      //printf("\nHit!\n");
      break; //first match
    }
  }    
  if(hit==1) {
  printf("\nFirst prediction: use x side as %.4lf\n", sa);
  printf("PS. Rem: x = z, deform equally. Fix y=%lf. \n", y);
  printf("Use option 2 after doing slabxyt for this 'xz deform =%.4lf, yhi = %lf' system to optimize\n", sa, y);
  return 0;
  }
  if(hit==0) 
  {
     printf("\nNo hits!\n");
     return 0;
  }
}

if (choice==2) 
{
  int ch;
  X:
  printf("\nChoose 1 or 2: \n 1 - Manual entry of hi(x,z) deform position\n 2 - Use initial prediction by program\n 3 - Compare two deformation (defxz, yhi=20) runs and scale\n");
  scanf("%d", &ch);
  if (ch==1) 
  {  //maual entry of first prediction that you used - useful in cases where we went badass and decided a diff value
    printf("Enter the previous approximation used for the 'xz deform yhi=20' run:\n");
    scanf("%lf", &sa);
    printf("Enter the density obtained from slabxyt (after fit [5:14]) for xz deformation run, wall @ yhi=20, x=z=%.4lf:\n", sa);
    scanf("%lf", &pdr); //rho deformation run
    sa2=((20*(pdr-paim))+(sa*(paim-pr)))/(pdr-pr); //20 is x or z
    printf("\nMore precise prediction: use x side as %.4lf\n", sa2);
    printf("PS. Rem: x = z, deform equally. Fix y=%lf. \n", y);
    printf("So, xz deform =%.4lf, yhi = %lf' system to optimize\n", sa2, y);
    return 0;
  }
  if (ch==2) 
  { //if you followed the same prediction suggested by computer for xz deform run
    for(i=(y/2);i<=sa;i+=0.00001) 
    { //repeat to find sa again
      if (((m/(i*i*y))-paim)<=0.00001) 
      { //precision control
        sa=i;
        break; //first match break out of loop
      }
    }
    //read the value from above 
   printf("Enter the density obtained from slabxyt (after fit [5:14]) for xz deformation run, wall @ yhi=20, x=z=%.4lf:\n", sa);
   scanf("%lf", &pdr); //rho deformation run
   sa2=((20*(pdr-paim))+(sa*(paim-pr)))/(pdr-pr); //20 is x or z
   printf("\nMore precise prediction: use x side as %.4lf\n", sa2);
   printf("PS. Rem: x = z, deform equally. Fix y=%lf. \n", y);
   printf("So, xz deform =%.4lf, yhi = %lf' system to optimize\n", sa2, y);
   return 0;
  }
  if (ch==3) 
  { //xa1,xa2,xa3 used only here
    printf("Enter the one of the previous approximation used for the 'xz deform yhi=20' run:\n");
    scanf("%lf", &xa1);
    printf("Enter the density obtained from slabxyt (after fit [5:14]) for xz deformation run, wall @ yhi=20, x=z=%.4lf:\n", xa1);
    scanf("%lf", &pdr1); //rho deformation run
    printf("Enter the second approximation used for the 'xz deform yhi=20' run:\n");
    scanf("%lf", &xa2);
    printf("Enter the density obtained from slabxyt (after fit [5:14]) for xz deformation run, wall @ yhi=20, x=z=%.4lf:\n", xa2);
    scanf("%lf", &pdr2); //rho deformation run
    double dxa12= xa1-xa2;
    double dpdr1aim=pdr1-paim;
    double dpdr12=pdr1-pdr2;
    double dpdr21=pdr2-pdr1;
    xa3=((dxa12*dpdr1aim)-(xa1*dpdr12))/dpdr21;   
    printf("\nMore precise prediction: use x side as %.4lf\n", xa3);
    printf("PS. Rem: x = z, deform equally. Fix y=%lf. \n", y);
    printf("So, xz deform =%.4lf, yhi = %lf' system to optimize\n", xa3, y);
    return 0;
  } //continue after ch=2 match 
   
  else //for ch
  {
    printf("Invalid option\n");
    goto X;
    return 0;
  }
} //end of choice 2 loop

else 
  { //for choice
    printf("Invalid option\n");
    goto Y;
    return 0;
  }
return 0;
}