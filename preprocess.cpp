#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <bits/stdc++.h>
using namespace std;
vector<vector<int> > tree;
vector<bool> vis;
vector<set<int> > adjlist;
set<int>::iterator it;
int verc,avgoutdeg,i,j,k,l,m,ei=0;
void dfs(int i){
  int j;
  std::queue<int> q;
  q.push(i);
  vis[i]=true;
  while(!q.empty()){
    i=q.front();
    q.pop();
    for (it=adjlist[i].begin(); it!=adjlist[i].end(); ++it){
      if(!vis[*it]){
        vis[*it]=true;
        q.push(*it);
        tree[i][j]=1;
      }
    }
  }
}
int main(int argc,char** argv)
{
 float ran;
 int maxlettersize;
 srand(time(0)); 
 int a,b;
 map<int,int> vertexid;
  set<int> temp;
  vector<int> temp2;
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
  adjlist[vertexid[a]].insert(vertexid[b]);
  adjlist[vertexid[b]].insert(vertexid[a]);
 }
 //printf("%d %d\n",lc,stoi(argv[1]));
 printf("%d\n",verc );
 vis.resize(verc);
 //scanf("%d",&maxlettersize);
 //return 0;
 maxlettersize=stoi(argv[1]);
 maxlettersize=std::max(maxlettersize,1); 
 maxlettersize=std::min(maxlettersize,26);
 int ei=0;
 for(i=0;i<=verc;i++){
  temp2.resize(adjlist[i].size());
  tree.push_back(temp2);
  printf("%d ",ei);
  ei+=adjlist[i].size();
 }
 printf("\n");
 for(i=0;i<verc;i++){
  ran=(int)rand();  
  printf("%c ",65+(int)ran%maxlettersize);
 }
 printf("\n");
 for(i=0;i<verc;i++)
  for (it=adjlist[i].begin(); it!=adjlist[i].end(); ++it)
    printf("%d ",*it);
 printf("\n");
 dfs(0);
 int cc=0;
for(i=0;i<verc;i++)
  for(j=0;j<adjlist[i].size();j++){
    printf("%d ",tree[i][j]);
    if(tree[i][j])
      cc++;
  }
 printf("\n");
// if(cc!=verc-1)
 // printf("\n\nFalse\n");
// ran=rand();
  // l=(int)((((int)ran)%verc));
  // printf("%d\n",l);

 return 0;
}