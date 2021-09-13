#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define pi 3.141592

/*function to calculate f(x,y,z) in y'= u */

double fder(double x,double y,double z)
{
    double fdx;
    fdx = z;
    return fdx;
}

/*function to calculate f(x,y,z) in y''= u' */

double gder(double p,double q,double r)
{
	double vdx;
	vdx = 6*cos(p)+2-q;
	return vdx;
}

/*function to calculate analytical  y(x) */

float lx(double f)
{
	double funy;
	funy = 2*cos(f)-pi*sin(f)+3*f*sin(f)+2;
	return funy;
}

/*function to calculate analytical in y'(x) */

double mx(double g)
{
	double funm;
	funm = sin(g)+cos(g)*(3*g-pi);
	return funm;
}

/*function which call above four function for calculation ynRK4, yn'RK4   */

double Rk4system(double x0, double u0, double v0, double h, double *un1, double *vn1)
{	
double un=u0, vn=v0, xn=x0;
double k1, k2, k3, k4, l1, l2, l3, l4;
	k1=h*fder(xn, un, vn);
	l1=h*gder(xn, un, vn);
	k2=h*fder(xn+h/2, un+k1/2, vn+l1/2);
	l2=h*gder(xn+h/2, un+k1/2, vn+l1/2);
	k3=h*fder(xn+h/2, un+k2/2, vn+l2/2);
	l3=h*gder(xn+h/2, un+k2/2, vn+l2/2);
	k4=h*fder(xn+h, un+k3, vn+l3);
	l4=h*gder(xn+h, un+k3, vn+l3);
	*un1=un+(k1+2*k2+2*k3+k4)/6;
	*vn1=vn+(l1+2*l2+2*l3+l4)/6;
}

int main()
{
	int j, k, N, M=6;
	double x0, y0, z0, i, h, L, steps[]= {2,5,10,15,20,25};                              
	double y1, z1, yErr[M], zErr[M], erysum, erzsum;
	FILE*out = fopen("Output_problem3.txt","w");
	for( M = 0; M <= 5; M++ )
	{
		x0 = 0, y0 = 4, z0 = -pi, N = 5, L = pi ;
		N =  steps[M]; 
	    h = (L - x0) / N; 
	
	
	
	        
	
	 
		 
	double  yval[N], zval[N], yvalrk4[N], zvalrk4[N];                   
	
  	yvalrk4[0] = yval[0] = y0; 
	zvalrk4[0] = zval[0] = z0;
    j=1;
    erysum = erzsum = 0;
    
    for (i=x0+h; i<=L+h/2; i=i+h)
     {
       
	       Rk4system(x0, yvalrk4[j-1], zvalrk4[j-1], h, &y1 , &z1 );
            yvalrk4 [j] = y1;
            zvalrk4 [j] = z1;
           yval[j] = lx(x0+h);
              erysum  = erysum + pow(yval[j]-yvalrk4[j],2);
           zval[j] = mx(x0+h);
		      erzsum  = erzsum + pow(zval[j]-zvalrk4[j],2);
			  x0 = x0+h;        
           j++;          
          } 
                      
        if (N == 5 || N == 10 || N == 20)
		{
    	fprintf(out,"\n\n\t   (x)\t      Analytical,y(x)\t\t yn (RK4)\t    Analytical,z(x)\t\t zn (RK4)\n");		
	
    	j=0;   
	    
		  for (i=0; i <= L+h/2; i=i+h)
		  {
      	 
        	fprintf(out,"\n\t  %.4f \t     %f\t\t    %f\t\t   %f\t\t   %f\t ", i , yval[j], yvalrk4[j], zval[j], zvalrk4[j]);
        	j++;
   	      } fprintf(out,"\n\n\n"); 
        } 

            yErr[k] =  sqrt(erysum)/N;
            zErr[k] =  sqrt(erzsum)/N;
            k++;
              
  }
        
    
fprintf(out,"\n\n \t (N) \t\t   eyL2\t\t\t       ezL2 \n");
for (k=1; k<=6; k++)
	{	   
		fprintf(out,"\n\t %.0f\t\t    %f \t        %f", steps[k-1], yErr[k], zErr[k]);	
	}
   	printf("\n\n The ouput is in Output_problem3.txt file");   

    fclose(out);
	return 0;
	
} 


        	
