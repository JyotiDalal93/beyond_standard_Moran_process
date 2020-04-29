// A population evolving as per the STANDARD MORAN process, with two allelic types (A and a) present at t = 0

#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<fstream>
#include<iomanip>
#include<math.h>

using namespace std;

//********************************************************************************************************************

// A function that gives random numbers generated from a power-law function
float U_lb, U_ub, nodel_para;
 
float power_law_distributed()   // g is the power-exponent here, a parameter of the model
{

float term1, term2, term3, term4, y, nodel_para_mod, x;

y = (rand()/(1.0*RAND_MAX));

nodel_para_mod = (nodel_para + 1);

term1 = pow(U_lb, (nodel_para_mod));
term2 = pow(U_ub, (nodel_para_mod));
term3 = (1/(nodel_para_mod));

term4 = (term2 - term1)*y + term1;
x = pow(term4, term3);

return x;
}
 
int main()
{

ofstream out;
//***************************************************************************
out.open("fixation_compre_500.txt");
int ind_death[1000], tot_A[110000];      // ind_death[], types[], w[], f[]: maximum size = 'pop_size'
char types[50000];                       // tot_A[]: maximum size = 't_max'  THINK OF A BETTER WAY 
float w[10000], f[10000];                // THINK OF A BETTER WAY 

int i, j, k, Init_A, t_max, t, n, n_exp, fix_A, Pop_size, U, count, ind_repro, ref, m;
float s, ran, l_bound, u_bound, sum, denom, sum_fit; 
char ch;

int df;
srand (time(NULL));

//***************************************************************************

//** Can also be taken as an input from users

t_max = 50000;       // time (or number of generations) given for population to evolve 
Pop_size = 500;     // population size
Init_A = 1;          // number of allele/type A one wants to start with, i.e., present at t = 0

n_exp = 10000;       // number of times the same population was allowed to evolve, keeping initial conditions same

s = 0.1;             // selection coefficient

out<<" m "<<" fix. prob. "<<" \n ";

//***************************************************************************
for(m = 2.0; m <= 5.0; m+=1.0)
{
fix_A = 0;
f[0]=0;
tot_A[0] = Init_A; 
df = 0;

U_lb = 2.0;
U_ub = (Pop_size);
nodel_para = -1.0*(m);

//***************************************************************************


for(n=1; n<(n_exp + 1); n++)
{
for(i=1; i<=(Pop_size); i++)
{

//***************************************************************************

//********** 'A' allele present in 'Init_A' copies in a population of size 'Pop_size' at t = 0 ********

if(i<=Init_A)
{
types[i] = 'A';    
w[i] = 1 + s;      // w[i] are the absolute fitnesses of the individuals
}
else
{
types[i] = 'a'; 
w[i] = 1;
}
}

//***************************************************************************

for(t=1; t<(t_max+1); t++)       // THINK OF A BETTER WAY
tot_A[t]=0;                      // tot_A[t]: total number of allele/type 'A' present after 't' number of generations

for(t=1; t<(t_max+1); t++)
{
denom = (Pop_size + s*tot_A[t-1]);
sum_fit=0;

for(i=1; i<=(Pop_size); i++)
{
sum_fit+=w[i];                  // ADDINg up to the w[i] 
f[i]=sum_fit;

if(types[i]=='A')
tot_A[t]+=1;
}

ran = (rand()/(1.0*RAND_MAX));

for(j=1; j<=(Pop_size); j++)
{
l_bound = f[j-1]/(denom); 
u_bound = f[j]/(denom); 

if(((l_bound < ran)&&(ran < u_bound))||(ran==l_bound))
{
ind_death[1] = j;
break;
}
}
 
//cout<<"ind_death[1]"<<ind_death[1]<<"\n";

//********************************************************************************************************************

// sampling/choosing 'U - 1' distinct individuals RANDOMLY who will not make it to the next generation 

 
U = int(power_law_distributed());     // U = 2, 3, 4 .....(Pop_size) 

//out<<" U: "<<U<<"    "<<"\n";

//U = 2;
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

//cout<<"ind_death[2]"<<ind_death[2]<<"\n";

/*cout<<"\n WITHOUT ##################################################################\n";

for(i=1; i<=(Pop_size); i++)
cout<<"  w[i] "<<w[i]<<" types[i] "<<types[i]<<"\n";

cout<<"\n WITHOUT ##################################################################\n";


cout<<"\n t: "<<t;
cout<<" b4  "<<" ind_death[i] "<<"   "<<" w[i] "<<"    "<<" types[i] "<<"\n";

for(i=1; i<=U; i++)
cout<<"   "<<ind_death[i]<<"   "<<w[ind_death[i]]<<"    "<<types[ind_death[i]]<<"\n";cout<<"ind_death[1]"<<ind_death[1]<<"\n";*/


ind_repro = ind_death[1];

for(i=2; i<=U; i++)
{
ref = ind_death[i];
if((types[ind_repro]=='A')&&(types[ref]=='a'))
tot_A[t]+=1;

else if((types[ind_repro]=='a')&&(types[ref]=='A'))
tot_A[t]-=1;

//else if((types[ind_repro]=='A')&&(types[ref]=='A'))
//tot_A[t]+=0;

//else
//tot_A[t]+=0;


w[ref] = w[ind_repro];
types[ref] = types[ind_repro];

//cout<<" ref = "<<ref<<" w[ref] "<<w[ref]<<" types[ref] "<<types[ref]<<"\n";

}

//cout<<" after  "<<" ind_death[i] "<<"   "<<" w[i] "<<"    "<<" types[i] "<<"\n";

//for(i=1; i<=(Pop_size); i++)
//cout<<"  w[i] "<<w[i]<<" types[i] "<<types[i]<<"\n";

//cout<<"\n";

//********************************************************************************************************************

//cout<<t<<" :  "<<"   "<<tot_A[t]<<" \n ";

if((tot_A[t]==0)||(tot_A[t]==(Pop_size)))   // one of the two alleles has reached fixation
{
df++;
break;
}
} // t
if((tot_A[t]==(Pop_size)))   // allele 'A' has reached fixation
fix_A+=1;

//cout<<"fix_A"<<fix_A<<"  ";
} // n_exp
//cout<<s<<"  "<<Pop_size<<"  "<<nodel_para<<" Fixation fixability "<<(1.0*fix_A)/(1.0*n_exp)<<"    "<<" \n ";

out<<m<<"  "<<(1.0*fix_A)/(1.0*n_exp)<<" \n ";
cout<<m<<"  "<<(1.0*fix_A)/(1.0*n_exp)<<" \n ";

cout<<" df "<<df<<" \n ";
}
out<<" N: "<<Pop_size<<" s " <<s<<"\n";
out.close();
}// main function

// Solution: TRY TO PRINT AFTER EACH STEP because one step is FUCKED UP!

// try it TOMORROW!

