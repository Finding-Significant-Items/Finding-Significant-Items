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
map <string ,int> B,C;
struct node {string x;int y;} p[10000005];
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
int cmp(node i,node j) {return i.y>j.y;}
int main(int argv, char **argc)
{

    int MEM,K;
	cin>>MEM>>K;
    int m=10000000; 
    // prepare for CM sketch plus a min-heap
    Sketchpheap *SH;
    // Insertion
    for (int i=1; i<=m; i++)
    {
		char aa[105];
        s[i]=Read();
        B[s[i]]++;
    }

	// prepare for real top-k frequent items
	int cnt=0;
    for (map <string,int>::iterator sit=B.begin(); sit!=B.end(); sit++)
    {
        p[++cnt].x=sit->first;
        p[cnt].y=sit->second;
    }
    sort(p+1,p+cnt+1,cmp);
    for (int i=1; i<=K; i++) C[p[i].x]=p[i].y;


    int sh_M;
    {
        for (sh_M=1; 432*K+sh_M*4*16<=MEM*1024*8; sh_M++);
        SH=new Sketchpheap(sh_M,K);
        for (int i=1; i<=m; i++)
        {
            SH->Insert(s[i]);
        }
        int SH_sum=0,SH_AAE=0; double SH_ARE=0;
        string SH_string; int SH_num;
        for (int i=0; i<K; i++)
        {
            SH_string=(SH->Query(i)).first; SH_num=(SH->Query(i)).second;
            SH_AAE+=abs(B[SH_string]-SH_num); SH_ARE+=abs(B[SH_string]-SH_num)/(B[SH_string]+0.0);
            if (C[SH_string]) SH_sum++;
        }

        printf("%dKB   top %d:\nAccepted: %d/%d  %.10f\nARE: %.10f\nAAE: %.10f\n\n",MEM,K,SH_sum,K,(SH_sum/(K+0.0)),SH_ARE/K,SH_AAE/(K+0.0));
    }
    return 0;
}
