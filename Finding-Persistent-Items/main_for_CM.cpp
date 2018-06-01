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
#include "params.h"
#include "Sketchpheap.h"
using namespace std;
map <string ,int> B,C,End,Value,Value2,mp;
struct node {string x;int y;} t[10000005];
ifstream fin("data",ios::in|ios::binary);
char a[105];
string Read()
{
    fin.read(a,4);
    a[4]='\0';
    string tmp=a;
    return tmp;
}
string s[10000005];
int WIN,Time[10000005];
int pd,CNT;
int CMP(node i,node j) {return i.y>j.y;}
int main()
{
    int MEM,K;
	cin>>MEM>>K>>pd;
    MEM/=2;
    int m=10000000;  // the number of flows in stream
    // prepare for sketch plus a min-heap
    Sketchpheap *SH;
    // Insertion
    for (int i=1; i<=m; i++)
    {
        s[i]=Read();
        Time[i]=i;
    }
    cout<<"preparing true flow"<<endl;
    for (int i=m; i>=1; i--) Time[i]-=Time[1];
    WIN = Time[m] / pd;
    int Last=0;
    for (int i=1; i<=m; i++)
    {
        string S=s[i]; int T=Time[i];
        if (T / WIN != Last / WIN)  // meet a new priod 
        {
            for (map<string,int> :: iterator sit = mp.begin(); sit!=mp.end(); sit++) End[sit->first]++; 
            mp.clear();
        }
        Last = T;
        mp[S]++;
    }
    for (map<string,int> :: iterator sit = mp.begin(); sit!=mp.end(); sit++) End[sit->first]++;  // execute the last period
    CNT=0; Value.clear(); Value2.clear();
	// the persistencies of all items are stored in 'End'
    for (map<string,int>::iterator sit = End.begin(); sit!=End.end(); sit++)  
    {
        t[++CNT].x=sit->first;
        t[CNT].y=sit->second;
    }
    cout<<CNT<<endl;
    sort(t+1,t+CNT+1,CMP);
    for (int i=1; i<=K; i++) Value[t[i].x]=t[i].y;  // get top-k persistent items
    for (int i=1; i<=CNT; i++) Value2[t[i].x]=t[i].y; 


    int sh_M;

    {
        for (sh_M=1; 432*K+sh_M*4*16<=MEM*1024*8; sh_M++);
        SH=new Sketchpheap(sh_M,K,MEM*1024*8/(int(log(pd)/log(2))));
        for (int i=1; i<=m; i++)
        {
            SH->Insert(s[i],Time[i]/WIN+1);
        }
        int SH_sum=0,SH_AAE=0; double SH_ARE=0;
        string SH_string; int SH_num;
        for (int i=0; i<K; i++)
        {
            SH_string=(SH->Query(i)).first; SH_num=(SH->Query(i)).second;
            SH_AAE+=abs(Value2[SH_string]-SH_num); SH_ARE+=abs(Value2[SH_string]-SH_num)/(Value2[SH_string]+0.0);
            if (Value[SH_string]) SH_sum++;
        }

        printf("%dKB   top %d:\nAccepted: %d/%d  %.10f\nARE: %.10f\nAAE: %.10f\n\n",MEM*2,K,SH_sum,K,(SH_sum/(K+0.0)),SH_ARE/K,SH_AAE/(K+0.0));
    }
    return 0;
}
