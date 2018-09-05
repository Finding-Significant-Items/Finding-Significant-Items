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
int M,m,MEM,pd,co_p,co_f,top_k,WIN;
const int G=8;

struct node {int C,odd,even; string FP; int count;} HK[35][500005];
int TT[35][500005];
BOBHash32 *bobhash; // hash function
int X,Y,Last_T,i,j;

const int ODD=1;

void Clear() // Initialization
{
    for (int i=0; i<G; i++)
      for (int j=0; j<M; j++) HK[i][j].C=HK[i][j].odd=HK[i][j].even=HK[i][j].count=0;
    Last_T=0; X=Y=0;
}

void Insert(string x,int T)
{
    int P_idx = Last_T / WIN % 2; // judge whether this period is an even-numbered period or an odd-numbered period
    if (T/WIN == Last_T/WIN) // if the item x does not cause the increment of period
    {
        // the pointer p moves clockwise
        while (TT[X][Y] <= T % WIN && Y!=M) {if (P_idx==ODD && HK[X][Y].even) HK[X][Y].C++,HK[X][Y].even=0; if (P_idx!=ODD && HK[X][Y].odd) HK[X][Y].C++,HK[X][Y].odd=0; X++; if (X==G) X=0,Y++;}
        if (Y==M) Y=0;
    }
    if (T/WIN != Last_T/WIN) // if the item x casues the increment of period
    {
        // similarly, the pointer p moves clockwise
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
    int p=bobhash->run(x.c_str(),x.size()) % M; // calculate the index of the mapped bucket
    int now=0,NOW=0; bool FLAG=false,FLAG2=false;
    for (int k=0; k<G; k++)
    if (HK[k][p].FP==x) {if (P_idx==ODD) HK[k][p].odd=1; else HK[k][p].even=1; HK[k][p].count++; FLAG=true;} // if x is found in that bucket
    else // calculate the smallest cell
    {
        if ((HK[k][p].C+HK[k][p].odd+HK[k][p].even)*co_p+HK[k][p].count*co_f<(HK[now][p].C+HK[now][p].even+HK[now][p].odd)*co_p+HK[now][p].count*co_f) now=k;
    }
    if (!FLAG) // if x is not found in that bucket
    {
        HK[now][p].C--; HK[now][p].count--; // Singificance Decrementing operation
        if (HK[now][p].C<0) HK[now][p].C=0;
        if (HK[now][p].count<0) HK[now][p].count=0;
        if (HK[now][p].C*co_f + HK[now][p].count*co_p<=0) // judge whether that item could be replaced with x
        {
            HK[now][p].FP=x;
            int MINC=10000,MINcount=10000;
            for (int k=0; k<G; k++)  // Long-tail Replacement technique
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
        HK[k][p].C+=HK[k][p].even+HK[k][p].odd; // the estimated persistency should be incremented if one of the flags is not equal to 0
        HK[k][p].even=HK[k][p].odd=0;
    }
}

map <string,int> mp,End,ans,MP,Value,Value2;

int CNT,CNTT;
char ID[10000005][12];
int Time[10000005];

struct NODE{string x;int y;} t[5000005],tt[5000005];
int CMP(NODE i, NODE j) {return i.y>j.y;}
void get_real()
{
    CNT=0; Value.clear(); Value2.clear();
    for (map<string,int>::iterator sit = MP.begin(); sit!=MP.end(); sit++)  // get the significances of all items
    {
        t[++CNT].x=sit->first; // means the ID
        t[CNT].y=End[sit->first]*co_p+sit->second*co_f; // means the significance
    }

    sort(t+1,t+CNT+1,CMP); // sort all values
    for (i=1; i<=top_k; i++) Value[t[i].x]=t[i].y;  // the real signifiacnt items are recorded in the Value
    for (i=1; i<=CNT; i++) Value2[t[i].x]=t[i].y;  // the significances of all items are recorded in the Value2
}
int main()
{
    bobhash = new BOBHash32(1000);
    // m denotes the number of flows.
    cin>>m;
    // co_f corresponds to the alpha in the paper,  co_p corresponds to the beta in the paper.
    cin>>co_f>>co_p;
    // top_k means the value of k in the paper.
    cin>>top_k;
    // MEM means the memomy size (KB)
    cin>>MEM;
    // pd means the number of periods
    cin>>pd;
    // there are m flows in total, every flow contains an ID and its arriving time.
    freopen("stack-new.txt","r",stdin);
    for (i=1; i<=m; i++)
    {
        scanf("%s",ID[i]);
        scanf("%d",&Time[i]);
    }
    // initialize the time.
    for (i=m; i>=1; i--) Time[i]-=Time[1];
    // WIN denotes the time span of one period
    WIN = Time[m] / pd;
    // Last denotes the arriving time of the previous flow
    int Last = 0;

    // prepare for the real top-k significant items
    for (i=1; i<=m; i++)
    {
        string s=ID[i]; int T=Time[i];
        if (T / WIN != Last / WIN)  // meet a new period
        {
            for (map<string,int> :: iterator sit = mp.begin(); sit!=mp.end(); sit++) End[sit->first]++;  //End[i] denotes the persistency of i
            mp.clear();
        }
        Last = T;
        mp[s]++; MP[s]++; // MP[i] denotes the real frequency of i, and mp[i] denotes whether i has appeared in this period
    }
    for (map<string,int> :: iterator sit = mp.begin(); sit!=mp.end(); sit++) End[sit->first]++;  // execute the last period
    // the process of preparing ends

    // the algorithm starts
    get_real(); // get the real top-k significant items
    int cnt=0;
    M=MEM*1024*8/G/64;  //the number of buckets, where G means the number of cells in each bucket
    for (i=0; i<M; i++)
        for (j=0; j<G; j++)
        {
            TT[j][i]=int((WIN)/(M*G+0.0)*cnt);  // TT[j][i] means the scanned time of j^{th} cell in the i^{th} bucket
            cnt++;
        }
    Clear(); // initialization

    for (i=1; i<=m; i++)
    {
        string s=ID[i];
        int T=Time[i];
        Insert(s,T); // an item S arrives at time T
    }
    do_it(); // execute the last period

    // get all the significances recorded in the lossy table
    CNTT=0;
    for (i=0; i<M; i++)
        for (j=0; j<G; j++)
        {
            tt[++CNTT].x=HK[j][i].FP;  // ID
            tt[CNTT].y=HK[j][i].C*co_p+HK[j][i].count*co_f; // significance
        }
    sort(tt+1,tt+CNTT+1,CMP); // get the estimated top-k significant items

    double PRE=0; double ARE=0.0; // calculate the precision and AAE of LTC
    for (i=1; i<=top_k; i++)
    {
        if (Value[tt[i].x]) PRE+=1.0; // means it is indeed a significant item
        ARE+=fabs(Value2[tt[i].x]-tt[i].y)/(Value2[tt[i].x]+0.0); // accumulate the value of the ARE
    }
    ARE/=top_k; PRE/=top_k;
    printf("precision: %.5f\nARE: %.5f\n",PRE,ARE);
    return 0;
}
