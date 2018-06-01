#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstring>
#include <map>
#include <fstream>
#include "BOBHASH32.h"
using namespace std;
ifstream fin("data",ios::in|ios::binary);
char a[105];
string Read()
{

    fin.read(a,4);
    a[4]='\0';
    string tmp=a;
    return tmp;
}
int M,m;
const int G=8;
struct node {int C,odd,even; string FP; int count;} HK[35][500005];
int TT[35][500005];
BOBHash32 *bobhash;
int X,Y,Last_T,i;

const int ODD=1;
int WIN;
int co_p,co_f,top_k,pw,MEM;
double co;
void Clear()
{
    for (int i=0; i<G; i++)
      for (int j=0; j<M; j++) HK[i][j].C=HK[i][j].odd=HK[i][j].even=HK[i][j].count=0;
    Last_T=0; X=Y=0;
}
void Insert(string x,int T)
{
    int P_idx = Last_T / WIN % 2;
    if (T/WIN == Last_T/WIN)
    {
        while (TT[X][Y] <= T % WIN && Y!=M) {if (P_idx==ODD && HK[X][Y].even) HK[X][Y].C++,HK[X][Y].even=0; if (P_idx!=ODD && HK[X][Y].odd) HK[X][Y].C++,HK[X][Y].odd=0; X++; if (X==G) X=0,Y++;}
        if (Y==M) Y=0;
    }
    if (T/WIN != Last_T/WIN)
    {
        for (int i = Last_T/WIN; i<T/WIN; i++)
        {
            while (Y!=M)
            {
                if (P_idx==ODD && HK[X][Y].even) HK[X][Y].C++,HK[X][Y].even=0;
                if (P_idx!=ODD && HK[X][Y].odd) HK[X][Y].C++,HK[X][Y].odd=0;
                X++; if (X==G) X=0,Y++;
            }
            Y=0;
            P_idx ^= 1;
        }
        while (TT[X][Y] <= T % WIN && Y!=M) {if (P_idx==ODD && HK[X][Y].even) HK[X][Y].C++,HK[X][Y].even=0; if (P_idx!=ODD && HK[X][Y].odd) HK[X][Y].C++,HK[X][Y].odd=0; X++; if (X==G) X=0,Y++;}
        if (Y==M) Y=0;
    }
    int p=bobhash->run(x.c_str(),x.size()) % M;
    int now=0,NOW=0; bool FLAG=false,FLAG2=false;
    for (int k=0; k<G; k++)
    if (HK[k][p].FP==x) {if (P_idx==ODD) HK[k][p].odd=1; else HK[k][p].even=1; HK[k][p].count++; FLAG=true;} else
    {
        if ((HK[k][p].C+HK[k][p].odd+HK[k][p].even)*co_p+HK[k][p].count*co_f<(HK[now][p].C+HK[now][p].even+HK[now][p].odd)*co_p+HK[now][p].count*co_f) now=k;
    }
    if (!FLAG)
    {
        HK[now][p].C--; HK[now][p].count--;
        if (HK[now][p].C<0) HK[now][p].C=0;
        if (HK[now][p].count<0) HK[now][p].count=0;
        if (HK[now][p].C*co_f + HK[now][p].count*co_p<=0)
        {
            HK[now][p].FP=x;
            int MINC=10000,MINcount=10000;
            for (int k=0; k<G; k++)
                if (k!=now) MINC=min(MINC,HK[k][p].C),MINcount=min(MINcount,HK[k][p].count);
            HK[now][p].C=max(0,MINC-1);
            HK[now][p].count=max(1,MINcount-1);
            HK[now][p].odd=HK[now][p].even=0;
            if (P_idx==ODD) HK[now][p].odd=1; else HK[now][p].even=1;
            FLAG2=true;
        }
    }
    Last_T=T;
}

void do_it()
{
    for (int k=0; k<G; k++)
      for (int p=0; p<M; p++)
    {
        HK[k][p].C+=HK[k][p].even+HK[k][p].odd;
        HK[k][p].even=HK[k][p].odd=0;
    }
}
map <string,int> mp,End,ans,MP,Value,Value2;
int j;

struct NODE{string x;int y;} t[5000005],tt[5000005];
int CNT,CNTT;
int CMP(NODE i, NODE j) {return i.y>j.y;}
char ID[10000005][12];
int Time[10000005];
const int pd=1000;
int main()
{
   char F[1005];
    bobhash = new BOBHash32(1000);
    m=10000000; // the number of items
	cin>>MEM>>top_k>>co_p>>co_f;
    for (i=1; i<=m; i++)
    {
        string SS=Read();
        for (j=0; j<4; j++) ID[i][j]=SS[j]; ID[i][SS.size()]='\0';
        scanf("%s",ID[i]);
        Time[i]=i;
    }
    for (i=m; i>=1; i--) Time[i]-=Time[1];
    WIN = Time[m] / pd;

    int cnt=0;
    M=MEM*1024*8/G/64;
    for (i=0; i<M; i++)
      for (j=0; j<G; j++)
      {
         TT[j][i]=int((WIN)/(M*G+0.0)*cnt);
         cnt++;
      }
    Clear();
    for (i=1; i<=m; i++)
    {
      string s=ID[i]; int T=Time[i];
      Insert(s,T);
    }
    do_it();
    CNTT=0;
    for (i=0; i<M; i++)
	  for (j=0; j<G; j++)
      {
        tt[++CNTT].x=HK[j][i].FP;
        tt[CNTT].y=HK[j][i].C*co_p+HK[j][i].count*co_f;
      }
    sort(tt+1,tt+CNTT+1,CMP);

    double PRE=0; double ARE=0.0;
    for (i=1; i<=top_k; i++)
    {
      if (Value[tt[i].x]) PRE+=1.0;
      ARE+=fabs(Value2[tt[i].x]-tt[i].y)/(Value2[tt[i].x]+0.0);
    }
    ARE/=top_k; PRE/=top_k;
    printf("PRE: %.10f\nARE: %.10f\n",PRE,ARE);
    return 0;
}
