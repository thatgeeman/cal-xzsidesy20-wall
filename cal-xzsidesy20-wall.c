/*
  Program
    cal-xzsidesy20-wall

  Description


  Keywords


  Related


  Author
    gg nov 29 2016
*/
#include <stdio.h>
#include <math.h> //gcc test.c -o test -lm
#include <stdlib.h>
int main() {
double sas, sa, sa2, m, y=20.0, va, v, paim, pr, pdr; //y set at 20
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

printf("Status: 1 or 2?\n 1 - After initial equilibration run with yhi=20\n 2 - Further scaling with xz deform, yhi=20\n");
scanf("%d", &choice);
if (choice==1)
{
  for(i=(y/2);i<=sa;i+=0.00001) { //weird for loop
    if (((m/(i*i*y))-paim)<=0.00001) { //precision control
      sa=i;//update x or z with the find
      hit++;
      //printf("\nHit!\n");
      break; //first match
      }
      //sa+=0.01; //precision factor - kills core in while!!
      else
      {
        printf(".\n"); //waiting for match
      }
    }
    if(hit==1) {
      printf("\nFirst prediction: use x side as %.4lf\n", sa);
      printf("PS. Rem: x = z, deform equally. Fix y=%lf. \n", y);
      printf("Use option 2 after doing slabxyt for this 'xz deform =%.4lf, yhi = %lf' system to optimize\n", sa, y);
      return 0;
    }
    if(hit==0) {
      printf("\nNo hits!\n");
      return 0;
    }
}
if (choice==2) {
  int ch;
  printf("Choose 1 or 2: \n 1 - Manual entry of hi(x,z) deform position\n 2 - Use initial prediction by program\n");
  scanf("%d", &ch);
  if (ch==1) {  //maual entry of first prediction that you used - useful in cases where we went badass and decided a diff value
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
  if (ch==2) { //if you followed the same prediction suggested by computer for xz deform run
  for(i=(y/2);i<=sa;i+=0.00001) { //repeat to find sa again, cos the computer forgets :P
    if (((m/(i*i*y))-paim)<=0.00001) { //precision control
      sa=i;
      break; //first match
      }
    }
    printf("Enter the density obtained from slabxyt (after fit [5:14]) for xz deformation run, wall @ yhi=20, x=z=%.4lf:\n", sa);
    scanf("%lf", &pdr); //rho deformation run
    sa2=((20*(pdr-paim))+(sa*(paim-pr)))/(pdr-pr); //20 is x or z
    printf("\nMore precise prediction: use x side as %.4lf\n", sa2);
    printf("PS. Rem: x = z, deform equally. Fix y=%lf. \n", y);
    printf("So, xz deform =%.4lf, yhi = %lf' system to optimize\n", sa2, y);
    return 0;
  }
  else {
    printf("Invalid option");
    return 0;
  }
} //end of choice 2 loop

else {
  printf("Invalid option\n");
  return 0;
}
return 0;
}
