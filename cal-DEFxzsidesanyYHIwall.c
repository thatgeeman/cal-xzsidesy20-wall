/*
  Program
    cal-DEFxzsidesanyYHIwall

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
    gg apr 09 2017 [any yhi]
    gg mai 22 2017 [replication]
*/
#include <stdio.h>
#include <math.h> //gcc test.c -o test -lm
#include <stdlib.h>
int main() {
double sas, sa, sa2, m, y, va, v, paim,ptaim, pr, pdr;
double xa1, xa2, xa3, pdr1, pdr2, repl;
double x,z;
int choice, hit=0;
printf("[read description] ~not harmful\n");
printf("enter nof particles (atoms) in bulk:\n");
scanf("%lf", &m); //normally nch*ncl
printf("enter yhi used for equilib run with wall\n");
scanf("%lf",&y);
printf("enter x and z used for equilib run with wall\n");
scanf("%lf %lf",&x, &z);
printf("Enter the density obtained from slabxyt (after fit [5:14]) for equilibration run, wall @ yhi:%lf\n", y);
scanf("%lf", &pr); //scanning double requires %lf; equilibrum run had all sides at 20,20,20
printf("Enter the volume obtained from thermo.data for bulk:\n");
scanf("%lf", &v); //scanning double requires %lf
printf("Replication case? 1=No 0=Yes:\n");
scanf("%d", &repl);
if (repl==1)
{
/*calculation*/
double i; //search from midpoint
paim=m/v; //bulk density = aim
va=m/pr; //apparent volume
sas=va/y; //product of apparent xz is sas
sa=sqrt(sas); //get x or z  : use -lm for compiling in gcc
/*searching for match*/
//double f=(m/(sa*sa*y)); //test-works
//printf("%.4lf\n", f); //test
Y:
printf("\n .. status: 1 or 2?\n 1 - after initial equilibration run with yhi %lf\n 2 - further scaling with xz deform, yhi %lf in comparison to equilib run\n", y, y);
scanf("%d", &choice);
if (choice==1)
{
  for(i=(y/2);i<=sa;i+=0.0001) 
  { //weird for loop
    if (((m/(i*i*y))-paim)<=0.0001) 
    {     //precision control
sa=i;//update x or z with the find
hit++;
//printf("\nHit!\n");
break; //first match
    }
  }    
  if(hit==1) {
  printf("\nFirst prediction: use x and z side as %.4lf\n", sa);
  printf("deform xz equally to %lf, fix y=%lf. \n", sa, y);
  printf("rerun program with option 2 after doing slabxyt for this 'xz deform =%.4lf, yhi = %lf' system to optimize\n", sa, y);
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
  printf("\nChoose 1 or 2: \n 1 - manual entry of hi (x,z) deform position\n 2 - use initial prediction by program\n 3 - Compare two deformation (defxz, yhi %lf) runs and scale (independant of equilib run)\n", y);
  scanf("%d", &ch);
  if (ch==1) 
  {  //maual entry of first prediction that you used - useful in cases where we went badass and decided a diff value
    printf("Enter the previous approximation used for the 'xz deform yhi=20' run:\n");
    scanf("%lf", &sa);
    printf("Enter the density obtained from slabxyt (after fit [5:14]) for xz deformation run, wall @ yhi %lf, x=z=%.4lf:\n",y,sa);
    scanf("%lf", &pdr); //rho deformation run
    sa2=((x*(pdr-paim))+(sa*(paim-pr)))/(pdr-pr); //20 is x or z
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
   sa2=((x*(pdr-paim))+(sa*(paim-pr)))/(pdr-pr); //20 is x or z
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
}

// if replication present requires scaling
if (repl==0)
{
int nrepl, nch2, chl2, m2;

printf("Enter planned no of Replication along one dim:\n");
scanf("%d", &nrepl);
printf("Enter nchains in minimized case:\n");
scanf("%d", &nch2);
printf("Enter ch length in minimized case\n");
scanf("%d", &chl2);
m2=nch2*chl2;
/*calculation*/
double i; //search from midpoint
ptaim=m/v; //bulk density = aim
paim=ptaim/nrepl;
va=m2/pr; //apparent volume pr is dens from wall case equilib
sas=va/y; //product of apparent xz is sas
sa=sqrt(sas); //get x or z  : use -lm for compiling in gcc
/*searching for match*/
//double f=(m/(sa*sa*y)); //test-works
//printf("%.4lf\n", f); //test
Yd:
printf("\n .. status: 1 or 2?\n 1 - after initial equilibration run with yhi %lf\n 2 - further scaling with xz deform, yhi %lf in comparison to equilib run\n", y, y);
scanf("%d", &choice);
if (choice==1)
{
  for(i=(y/2);i<=sa;i+=0.0001) 
    { //weird for loop
if (((m2/(i*i*y))-paim)<=0.0001) 
   {     //precision control
    sa=i;//update x or z with the find
    hit++;
    //printf("\nHit!\n");
    break; //first match
   }
     }    
   if(hit==1) {
printf("\nFirst prediction: use x and z side as %.4lf\n", sa);
printf("deform xz equally to %lf, fix y=%lf. \n", sa, y);
printf("rerun program with option 2 after doing slabxyt for this 'xz deform =%.4lf, yhi = %lf' system to optimize\n", sa, y);
return 0;
  }
  if(hit==0) 
  {   printf("\nNo hits!\n");
return 0;
   }
  }

  if (choice==2) 
  {
    int ch;
Xd:
  printf("\nChoose 1 or 2: \n 1 - manual entry of hi (x,z) deform position\n 2 - use initial prediction by program\n 3 - Compare two deformation (defxz, yhi %lf) runs and scale (independant of equilib run)\n"
  , y);
    scanf("%d", &ch);
if (ch==1) 
  {  //maual entry of first prediction that you used - useful in cases where we went badass and decided a diff value
printf("Enter the previous approximation used for the 'xz deform yhi=20' run:\n");
    scanf("%lf", &sa);
  printf("Enter the density obtained from slabxyt (after fit [5:14]) for xz deformation run, wall @ yhi %lf, x=z=%.4lf:\n",y,sa);
scanf("%lf", &pdr); //rho deformation run
    sa2=((x*(pdr-paim))+(sa*(paim-pr)))/(pdr-pr); //20 is x or z
  printf("\nMore precise prediction: use x side as %.4lf\n", sa2);
printf("PS. Rem: x = z, deform equally. Fix y=%lf. \n", y);
    printf("So, xz deform =%.4lf, yhi = %lf' system to optimize\n", sa2, y);
  return 0;
    }
if (ch==2) 
  { //if you followed the same prediction suggested by computer for xz deform run
     for(i=(y/2);i<=sa;i+=0.00001) 
   { //repeat to find sa again
  if (((m2/(i*i*y))-paim)<=0.00001) 
 { //precision control
    sa=i;
     break; //first match break out of loop
}
    }
  //read the value from above 
    printf("Enter the density obtained from slabxyt (after fit [5:14]) for xz deformation run, wall @ yhi=20, x=z=%.4lf:\n", sa);
  scanf("%lf", &pdr); //rho deformation run
     sa2=((x*(pdr-paim))+(sa*(paim-pr)))/(pdr-pr); //20 is x or z
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
  goto Xd;
return 0;
  }
  } //end of choice 2 loop

  else 
    { //for choice
  printf("Invalid option\n");
goto Yd;
    return 0;
}

}
//end of repl option 0
return 0;
}
