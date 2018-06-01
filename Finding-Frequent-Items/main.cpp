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
#include "ssummary.h"
#include "spacesaving.h"
#include "LossyCounting.h"
#include "CSS.h"
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
int cmp(node i,node j) {return i.y>j.y;}
int main()
{
    int MEM,K;
    cin>>MEM>>K;

    cout<<"MEM="<<MEM<<"KB"<<endl<<"Find top"<<K<<endl<<endl;
    cout<<"prepare for all algorithms"<<endl;
    int m=10000000;  // the number of flows

    // spacesaving
    int ss_M;
    for (ss_M=1; 432*ss_M<=MEM*1024*8; ss_M++);
    spacesaving *ss; ss=new spacesaving(ss_M,K);

    // LossyCounting
    int LC_M;
    for (LC_M=1; 227*LC_M<=MEM*1024*8; LC_M++);
    LossyCounting *LC; LC=new LossyCounting(K);

    // CSS
    int css_M;
    for (css_M=1; 179*css_M+4*css_M*log(css_M)/log(2)<=MEM*1024*8; css_M++);
    CSS *css; css=new CSS(css_M,K); css->clear();
   
   // Insertion
    for (int i=1; i<=m; i++)
	{
		string s=Read();
		B[s]++;
		ss->Insert(s);
		LC->Insert(s,i/LC_M); if (i%LC_M==0) LC->clear(i/LC_M);
		css->Insert(s);
	}
	ss->work();
	LC->work();
	css->work();

	// prepare for true flow
	int cnt=0;
    for (map <string,int>::iterator sit=B.begin(); sit!=B.end(); sit++)
    {
        p[++cnt].x=sit->first;
        p[cnt].y=sit->second;
    }
    sort(p+1,p+cnt+1,cmp);
    for (int i=1; i<=K; i++) C[p[i].x]=p[i].y;

    // calculate PRE, ARE, AAE

    int LC_sum=0,LC_AAE=0; double LC_ARE=0;
    string LC_string; int LC_num;
    for (int i=0; i<K; i++)
    {
        LC_string=(LC->Query(i)).first; LC_num=(LC->Query(i)).second;
        LC_AAE+=abs(B[LC_string]-LC_num); LC_ARE+=abs(B[LC_string]-LC_num)/(B[LC_string]+0.0);
        if (C[LC_string]) LC_sum++;
    }

    int ss_sum=0,ss_AAE=0; double ss_ARE=0;
    string ss_string; int ss_num;
    for (int i=0; i<K; i++)
    {
        ss_string=(ss->Query(i)).first; ss_num=(ss->Query(i)).second;
        ss_AAE+=abs(B[ss_string]-ss_num); ss_ARE+=abs(B[ss_string]-ss_num)/(B[ss_string]+0.0);
        if (C[ss_string]) ss_sum++;
    }

    int css_sum=0,css_AAE=0; double css_ARE=0;
    string css_string; int css_num;
    for (int i=0; i<K; i++)
    {
        css_string=(css->Query(i)).first; css_num=(css->Query(i)).second;
        css_AAE+=abs(B[css_string]-css_num); css_ARE+=abs(B[css_string]-css_num)/(B[css_string]+0.0);
        if (C[css_string]) css_sum++;
    }
    printf("LossyCounting:\nAccepted: %d/%d  %.10f\nARE: %.10f\nAAE: %.10f\n\n",LC_sum,K,(LC_sum/(K+0.0)),LC_ARE/K,LC_AAE/(K+0.0));
    printf("spacesaving:\nAccepted: %d/%d  %.10f\nARE: %.10f\nAAE: %.10f\n\n",ss_sum,K,(ss_sum/(K+0.0)),ss_ARE/K,ss_AAE/(K+0.0));
    printf("CSS:\nAccepted: %d/%d  %.10f\nARE: %.10f\nAAE: %.10f\n\n",css_sum,K,(css_sum/(K+0.0)),css_ARE/K,css_AAE/(K+0.0));
    return 0;
}
