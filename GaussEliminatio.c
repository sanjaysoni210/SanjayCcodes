#include<stdio.h>
#include <stdlib.h>
#include<math.h>

int main()
{
	int i,j, k,p, q,x,y;
	int r=0,c=1;

	char ch1; //ch1 is used to read the file
	FILE *fp1;  //  initializing a file pointer
	fp1= fopen("file.txt","r"); 
	
	while((ch1=fgetc(fp1))!=EOF)  
	{
		if (ch1=='\t'&& r==0) 
		c++;  
		else if (ch1=='\n')
		r++; 
	}
	printf("For the augmented matrix A|b Rows=%d\tColumns=%d\n",r,c);
	printf("\nThe augmented matrix A|b is\n\n");
	fclose(fp1);
	float A[r][c]; 
	fp1=fopen("file.txt","r"); //opening file in read mode where the matrix data was present
	while(!feof(fp1)) //while the end of file is not reached it will be scanned
	{
		for(i=1;i<=r;i++)
		for(j=1;j<=c;j++)
		fscanf(fp1,"%f",&A[i][j]); 
	}
	for(i=1;i<=r;i++)
	{
		for(j=1;j<=c;j++)
		{
		printf("%f\t",A[i][j]);  //matrix A|b is shown in the screen
	}
	printf("\n");
	}
	fclose(fp1);
	float temp[20],ratio,solution[20],sum;
	for( j = 1; j <= c-1; j++){                                                   
        for( i = 1; i <= r; i++){
        	k=i+1; y=i;
     while (k <= r){
    if (fabs(A[k][i]) > fabs(A[y][i])){                          
        y=k;}
         k++;}
		x=i; q=y;
	 for (p = 1; p <= c; p++){                                 
    	 temp[p]  =  A[q][p];
        	A[q][p]  =  A[x][p];
    	    A[x][p]  =  temp[p];	}
	   if ( i>j){
	   ratio = A[i][j] / A[j][j];
    	 for( k = 1; k <= c; k++){
	     A[i][k] = A[i][k] - ratio * A[j][k];
		           }
		     	}
		}
	}
	solution[r] = A[r][c] / (A[r][c-1]);
    
    printf("\n Matrix Reduced to Upper Triangular Form: \n\n");    

	for( i = 1; i <= r; i++){
    	for(j = 1; j <= c; j++){
    	    printf("%f\t", A[i][j]);
	    }
	    printf("\n");
	}
	if (A[r][c-1]==0 && A[r][c]!=0){

	    printf("\n Solution does not exist\n");                             

	} else {

	for(i = r-1; i >= 1; i--){                                                     
		sum = 0;
		for(j = i+1; j<=c-1; j++){
			sum = sum + A[i][j] * solution[j];
		}
		solution[i] = ( A[i][c]-sum) / A[i][i];
	}
	printf("\nSolution: \n");                          
    for (i = 1; i <= r; i++)
	printf("\n The value of x(%d)= %f\t",i,solution[i]);
    }
	return 0;
}

