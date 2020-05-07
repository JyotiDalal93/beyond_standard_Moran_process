#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<fstream>
#include<iomanip>
#include<math.h>

using namespace std;

float U_lb, U_ub, pow_exp;
 
float power_law_distributed()                   // A function that gives random numbers generated from a power-law function 
{

float term1, term2, term3, term4, y, model_para_mod, x;      // pow_exp: a parameter of the model

y = (rand()/(1.0*RAND_MAX));

model_para_mod = (pow_exp + 1);

term1 = pow(U_lb, (model_para_mod));
term2 = pow(U_ub, (model_para_mod));
term3 = (1/(model_para_mod));

term4 = (term2 - term1)*y + term1;
x = pow(term4, term3);

return x;
}

 
int main()
{

ofstream out;
out.open("fixation_pop_size_1000_s_0.5_1.txt");
//out.open("fixation_pop_size_500_s_0.1_2.txt");
int ind_death[1000];                // ind_death[], types[], w[], f[]: maximum size = 'pop_size'
char types[50000]; 
                     
int tot_A[110000];                  // tot_A[]: maximum size = 't_max'  

float w[10000], f[10000];                

int i, j, k, Init_A, t_max, t, n, no_exp, fix_A, Pop_size, U, count, ind_repro, ref;
float s, ran, l_bound, u_bound, sum, denom, sum_fit, model_para; 
char ch;

int no_fixation;
srand (time(NULL));

//******** Can also be taken as an input from users

t_max = 55000;                   // time (or number of generations) given for population to evolve 
Pop_size = 1000;                  // population size
Init_A = 1;                      // number of allele/type A one wants to start with, i.e., present at t = 0

no_exp = 10000;                  // number of times the same population was allowed to evolve, keeping initial conditions same

s = 0.5;                         // selection coefficient

out<<setw(15)<<" model_parameter "<<setw(15)<<" fix. prob. "<<" \n ";

//***************************************************************************

//for(model_para = 4.75; model_para <= 4.0; model_para+=0.25)
for(model_para = 6.0; model_para <= 6.0; model_para+=0.5)
{
fix_A = 0;
f[0]=0;
tot_A[0] = Init_A; 
no_fixation = 0;

U_lb = 2.0;
U_ub = (Pop_size);
pow_exp = -1.0*(model_para);            // power_exponent of the progeny distribution

for(n=1; n<(no_exp + 1); n++)
{
for(i=1; i<=(Pop_size); i++)
{

//********** 'A' allele present in 'Init_A' copies in a population of size 'Pop_size' at t = 0 ********

if(i<=Init_A)   
{
types[i] = 'A';     
w[i] = 1 + s;                    // w[i] are the absolute fitnesses of the individuals
}
else
{
types[i] = 'a'; 
w[i] = 1;
}
}

//***************************************************************************

for(t=1; t<(t_max+1); t++)       
tot_A[t]=0;
                      
for(t=1; t<(t_max+1); t++)
{
denom = (Pop_size + s*tot_A[t-1]);
sum_fit=0;

for(i=1; i<=(Pop_size); i++)
{
sum_fit+=w[i];                     // ADDINg up all the w[i] 
f[i]=sum_fit;

if(types[i]=='A')
tot_A[t]+=1;                       // tot_A[t]: total number of allele/type 'A' present after 't' number of generations
}

ran = (rand()/(1.0*RAND_MAX));

for(j=1; j<=(Pop_size); j++)
{
l_bound = f[j-1]/(denom); 
u_bound = f[j]/(denom); 

if(((l_bound < ran)&&(ran < u_bound))||(ran==l_bound))
{
ind_repro = j;                    // THE individual chosen for reproduction
break;
}
}

//********************************************************************************************************************

U = int(power_law_distributed());     // U = 2, 3, 4 .....(Pop_size) 

/* Following part: sampling/choosing 'U - 1' distinct individuals RANDOMLY from 'N - 1' other individuals (except for the individual chosen for the reproduction process) 
// the 'U - 1' distinct individuals chosen will not make it to the next generation, i.e., will DIE in this step

//U = 2;*/
ind_death[1] = ind_repro; 

i = 2;
while(i<=U)
{
ch='F';

while(ch=='F')
{
count = 0;

ind_death[i] = int(((rand()/(1.0*RAND_MAX)))*(Pop_size + 1));  
                   
if((ind_death[i]>0)&&(ind_death[i]<(Pop_size + 1)))
{
for(k=1; k<(i); k++)
{

if(ind_death[i]!=ind_death[k])
{
count++;
if(count==(i-1))
{
ch = 'T';
i++;
break;
}
}

}
}
}
}

// Sampling: DONE

/* Updating the number of allele 'A' present*/

for(i=2; i<=U; i++)
{
ref = ind_death[i];
if((types[ind_repro]=='A')&&(types[ref]=='a'))
tot_A[t]+=1;

else if((types[ind_repro]=='a')&&(types[ref]=='A'))
tot_A[t]-=1;

w[ref] = w[ind_repro];
types[ref] = types[ind_repro];
}


/* Updating the number of allele 'A' present: done */


if((tot_A[t]==0)||(tot_A[t]==(Pop_size)))   // two alleles has reached fixation, given enough time
{
no_fixation++;                              // no_fixation = no_exp SHOULD HOLD, to make sure that one of two alleles has reached fixation
break;
}

} // t


if((tot_A[t]==(Pop_size)))   // allele 'A' has reached fixation
fix_A+=1;

} // no_exp

out<<setw(15)<<model_para<<setw(15)<<(1.0*fix_A)/(1.0*no_exp)<<" \n ";
cout<<setw(15)<<model_para<<setw(15)<<(1.0*fix_A)/(1.0*no_exp)<<" \n ";
cout<<" no_fixation "<<no_fixation<<" \n ";
}
out<<"\n\n";
out<<" #allele A (at t = 0) "<<Init_A<<"\n";
out<<" population size: "<<Pop_size<<"\n";
out<<" selection coefficient: " <<s<<"\n";
out<<" #averagings done "<<no_exp<<"\n";
out.close();
}// main function


