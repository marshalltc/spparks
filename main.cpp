#include "spparks/spparks_files/library.h"
#include "lammps/library.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <cstdlib>
using namespace std;

int main(){
 
/*-----------run spparks and open lammps------------*/
	//Initialize an instance of spparks and return pointer to it
	char **argv;
	void *ptr;
	spparks_open_no_mpi(0,argv,&ptr);

	//Initialize an instance of lammps and return pointer to it
	char **argv2;
	void *ptr2;
	lammps_open_no_mpi(0,argv2,&ptr2);
	
	//Process input script from file 
	char filename[]="spparks/in.spparks";
	spparks_file(ptr,filename);

	//deposit layers using modified deposition command where the last number of the command is the particle id
	for(int i = 1; i < 101; i++){
		char depo[]  = "deposition 25 0 0 -1 2 9 26 2";
		char run[]   = "run 1";

		if(i%30 == 0) depo[sizeof(depo)-2] ='2'; 
		else depo[sizeof(depo)-2] = '1';

		spparks_command(ptr,depo);
		spparks_command(ptr,run);
	}

	//Close spparks
	spparks_close(ptr);
/*------------------------------------------------------------------------------------------*/
	
/*------------make a data file to load into lammps from spparks atom positions-------------*/	
	FILE * readfile = fopen("data/dumpster100.txt","r");
	FILE * writefile = fopen("data/data.lammps","w");
	
	//write the header 
	char myheader[] = "\n                   \n2 atom types\n\n0 50 xlo xhi\n0 50 ylo yhi\n0 50 zlo zhi\n\nMasses\n\n1 2\n2 3\n\nAtoms\n\n";
	fputs(myheader,writefile);
	
	//to make the id number we count atoms read from file that are a certain atom type
	int count = 1;

	//copy lines and count atoms
	while(!feof(readfile)){
		char newline[50];
		char oneline[50];
		char *x;
		char *y;
		char *z;
		char *i2;
		
		fgets(oneline,50,readfile);
		
		strtok(oneline, " ");
		strtok(NULL, " ");
		x = strtok(NULL, " ");
		y = strtok(NULL," ");
		z = strtok(NULL, " ");
		i2 = strtok(NULL, " ");
		
		int g;
		if(i2 != NULL)g	= atoi(i2);
		if(g == 1 || g == 2){
			char id[7];
			sprintf(id,"%i",count);
			strcpy(newline,id);
			strcat(newline," ");	
			strcat(newline,i2);
			strcat(newline," ");
			strcat(newline,x);
			strcat(newline," ");
			strcat(newline,y);
			strcat(newline," ");
			strcat(newline,z);
			strcat(newline," ");
			strcat(newline,"\n");

			fputs(newline,writefile);
			count++;   
		}
	}
	
	//write the final atom count to the header for lammps data file
	fseek(writefile,0,SEEK_SET);
	char natoms[50];
	sprintf(natoms,"\n%i atoms \n", count-1);
	fputs(natoms,writefile);
	
	fclose(readfile);
	fclose(writefile);
/*------------------------------------------------------------------------*/



/*-----------execute lammps script------------*/
	//Process input script from file 
	char filename2[]="lammps/in.lammps";
	lammps_file(ptr2,filename2);

	//Close lammps
	lammps_close(ptr2);
/*---------------------------------------------*/


/*--------filter xrd data----------------------*/

FILE * readfile2 = fopen("data/xrd_data","r");
FILE * writefile2 = fopen("data/xrd_data.csv","w");

while(!feof(readfile2)){
		char newline[50];
		char oneline[50];
		char *a;
		char *b;
		char *c;
		char *d;
		
		fgets(oneline,50,readfile2);
		
		a = strtok(oneline, " ");
		b = strtok(NULL, " ");
		c = strtok(NULL, " ");
		d = strtok(NULL," ");
		
		 
		int g = 0;
		if(a != NULL)g = atoi(a); 
		
		if(g < 20 && g > 0){
			
			strcpy(newline,a);
			strcat(newline,",");	
			strcat(newline,b);
			strcat(newline,",");
			strcat(newline,c);
			strcat(newline,",");
			strcat(newline,d);
			

			fputs(newline,writefile2);
		}
	}
	
	
	fclose(readfile2);
	fclose(writefile2);

/*----------------------------------------------*/







}
