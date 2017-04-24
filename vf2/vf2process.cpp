#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <bits/stdc++.h>
using namespace std;
vector<vector<int> > tree;
vector<bool> vis;
vector<vector<int> > adjlist;
int verc,avgoutdeg,i,j,k,l,m,ei=0;

int main(int argc,char** argv)
{
 float ran;
 int maxlettersize;
 srand(time(0)); 
 int a,b;
 map<int,int> vertexid;
  vector<int> temp;
 int verc=0;
 int lc=0;
 while(1){
  if(scanf("%d\t%d ",&a,&b)<2)
    break;
  lc++;
  if(vertexid.find(a)==vertexid.end())
    adjlist.push_back(temp),vertexid[a]=verc++;
  if(vertexid.find(b)==vertexid.end())
    adjlist.push_back(temp),vertexid[b]=verc++;
  adjlist[vertexid[a]].push_back(vertexid[b]);
 }
 //printf("%d %d\n",lc,stoi(argv[1]));
 printf("%d\n",verc );
 int ei=0;
 for(i=0;i<=verc;i++){
  ei+=adjlist[i].size();
 }
 printf("%d\n",ei);
 for(i=0;i<verc;i++)
  for(j=0;j<adjlist[i].size();j++)
    printf("%d %d\n",i,adjlist[i][j]);

// ran=rand();
  // l=(int)((((int)ran)%verc));
  // printf("%d\n",l);

 return 0;
}