#include<stdio.h>
#include<math.h>
#include<stdlib.h>

/*function to calculate f(x,y) in y'=f(x,y)*/

double deriv(float x,float y)
{
    float fd;
    fd=y+2*x-x*x;
    return fd;
}

//function to calculate g(x)//
double  funyx(double x)
{
	return  x*x + exp(x);
}

/*function to calculate y'f(x,y) by euler method*/

double  Euler(double x0, double y0, double L, double N)
{
	double ye, y_der;
	double h=L/N;
	y_der=deriv(x0,y0);
    ye=y0+h*y_der;
    return ye;
}

/*function to calculate y'=f(x,y) by Mid point method*/

double  Midpt(double x0, double y0, double y1, double L, double N)
{
	double h=L/N;
	double y2,y1d;
	double x1=x0+h;
    y1d=deriv(x1,y1);
    y2=y0+2*h*y1d;
	return y2;
}

/*function to calculate y'=f(x,y) by Runge kutta second order method*/

double  RK2(double x0, double y0, double L, double N)
{
	double h=L/N;
	double k1, k2, yrk2;
	k1=h*deriv(x0,y0);
	k2=h*deriv(x0+h,y0+k1);
	yrk2=y0+(k1+k2)/2;
	return yrk2;
}

/*function to calculate y'f(x,y) by Runge kutta fourth order method*/

double  RK4(double x0, double y0, double L, double N)
{
	double h=L/N;
	double p1, p2, p3, p4, yrk4;
	double p=x0, q=y0;
	p1=h*deriv(p,q);
	p2=h*deriv(p+h/2,q+p1/2);
	p3=h*deriv(p+h/2,q+p2/2);
	p4=h*deriv(p+h,q+p3);
	yrk4=q+(p1+2*p2+2*p3+p4)/6;
	return yrk4;
}


int main()
{
	FILE *out;
	out = fopen("Output_problem1.txt","w");
	double x0, y0, i, h, L, steps[]={2,5,10,15,20,25};                              
	double errE[6], errMP[6], errRk2[6], errRk4[6];
	int   k,l=1, N;
	
	  x0 = 0 , y0 = 1, N = 5, L = 5 ;      
	
	 for(k=1; k<=6; k++)
	 {
     	N=steps[k-1];
     
	    h=(L - x0)/N ;	  
 
     double yval[N];  yval[0]=funyx(x0);	
     double yeul[N],  erroE=0; yeul[0] = y0; int j=1;
     for (i=h; i<=L; i=i+h)
	 {
	            
	            yeul[j] = Euler(i-h,yeul[j-1],L,N);
             	yval[j] = funyx(i);
             	erroE        = erroE + pow((yval[j]-yeul[j]),2);
		        j++;				
        }
        errE[k] = sqrt(erroE)/N;
        
    double yvalMP[N];   yvalMP[0] = y0;  yvalMP[1] = yeul[1];  j=1; erroE = pow((yval[j]-yeul[j]),2);
    
        for (i = h; i<L; i=i+h)
		{
	
	            yvalMP[j+1] = Midpt(i-h,yvalMP[j-1],yvalMP[j],L,N); 
	            erroE      = erroE + pow((yval[j+1]-yvalMP[j+1]),2);
		        j++;
        } 
		errMP[k] = sqrt(erroE)/N;  
	   
	double yvalRk2[N],   e_rk2=0;    yvalRk2[0] = y0;    j = 1;
    
        for (i = h; i<=L; i=i+h)
		{
	            
				yvalRk2[j] = RK2(i-h,yvalRk2[j-1],L,N);
				e_rk2      = e_rk2 + pow((yval[j]-yvalRk2[j]),2);
		        j++;
        }
        errRk2[k] = sqrt(e_rk2)/N;

	
	double yvalRk4[N],    errk4=0;    yvalRk4[0] = y0;    j = 1;
     
         for (i = h; i<=L; i=i+h)
		 {
	           
			    yvalRk4[j] = RK4(i-h,yvalRk4[j-1],L,N);
			    errk4      = errk4 + pow((yval[j]-yvalRk4[j]),2);
	            j++;
    	}
    	errRk4[k] = sqrt(errk4)/N;
	
	if (N == steps[1] || N == steps[2] || N == steps[4])
	{
 
    	fprintf(out,"\n  x\t	Analytical(x)\t	yn(euler))\t yn(Midpoint)\t	 yn(Rk2)\t	yn(Rk4)\n");
        j=0;  
	for (i=0; i<=L; i=i+h)
	{
      	 fprintf(out,"\n %.4f\t\t %.5f \t   %.5f\t    %.5f\t       %.5f\t         %.5f", i , yval[j], yeul[j], yvalMP[j], yvalRk2[j], yvalRk4[j]);
    	j++;
    }
    fprintf(out,"\n\n\n"); 
    }
	 
	    
}
fprintf(out,"\n\n \t (N) \t\t  Euler\t\t  MidPoint \t\tRK2 \t\t  RK4 \n");
for (k=1; k<=6; k++)
{
	fprintf(out,"\n\t %.0f\t\t  %f \t    %f\t        %f\t %f", steps[k-1], errE[k], errMP[k], errRk2[k],errRk4[k] );
	
}
	fclose(out);
	printf("output is in Output_problem1 txt file");
	return 0;
	
}  
    	
	
	


	