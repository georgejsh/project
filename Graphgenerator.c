#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <bits/stdc++.h>
int *ver,*edge,*count,*tree;
bool * vis;
int verc,avgoutdeg,i,j,k,l,m,ei=0;
void dfs(int i){
  int j;
  std::queue<int> q;
  q.push(i);
  vis[i]=true;
  while(!q.empty()){
    i=q.front();
    q.pop();
    for(j=ver[i];j<ver[i+1];j++){
      if(!vis[edge[j]]){
        vis[edge[j]]=true;
        q.push(edge[j]);
        tree[j]=1;
      }
    }
  }
}
int main()
{
 float ran,var,maxvar,curravg=0,currvar=0;
 int maxlettersize;
 srand(time(0)); 
 //printf("Enter the Vertex Count: ");
 scanf("%d",&verc);
 ver=(int *)malloc(sizeof(int)*(verc+1));
 //tree=(int *)malloc(sizeof(int)*(verc+1));
 vis=(bool *)malloc(sizeof(bool)*(verc+1));
 count=(int *)malloc(sizeof(int)*(verc+1));
 memset(count,0,sizeof(int)*(verc+1));
  memset(vis,false,sizeof(int)*(verc+1));
// printf("Enter the Average OutDegree: ");
 scanf("%d",&avgoutdeg);
 
 if(avgoutdeg>=verc)
  avgoutdeg=verc-1;
 //printf("Variance: ");
 scanf("%f",&var);
 //var*=verc;
 scanf("%d",&maxlettersize);
 maxlettersize=std::max(maxlettersize,1); 
 maxlettersize=std::min(maxlettersize,26);
 edge=(int *)malloc(sizeof(int)*(verc*(avgoutdeg+var)));
 tree=(int *)malloc(sizeof(int)*(verc*(avgoutdeg+var)));
 memset(tree,0,sizeof(int)*(verc*(avgoutdeg+var)));
 maxvar=var;//floor(sqrt(var));
 printf("%d\n",verc);
 for(i=0;i<verc;i++)
 {
  ran=(int)rand();
  //printf("%f",ran);
  if(maxvar)
   j=avgoutdeg-maxvar+(int)ran%((int)(2*maxvar));
  else
   j=avgoutdeg;
  j=std::min(j,verc);
  j=std::max(j,0);
  ver[i]=curravg;
  curravg+=j;
  //currvar+=(j-avgoutdeg)*(j-avgoutdeg);
  //printf("%d %d\n",i,j);
  for(k=0;k<j;)
  {
   ran=rand();
   l=(int)((((int)ran)%verc));//%verc;
   if(count[l]==0 && l!=i)
    edge[ei++]=l,k++,count[l]=1;//,printf("%d ",l);;
   
  }
  for(k=ver[i];k<curravg;k++)
   count[edge[k]]=0;
  //printf("\n");
 }
 ver[i]=curravg;
// printf("%d\n",ei);
 for(i=0;i<=verc;i++)
  printf("%d ",ver[i]);
 printf("\n");
 for(i=0;i<verc;i++){
  ran=(int)rand();  
  printf("%c ",65+(int)ran%maxlettersize);
 }
 printf("\n");
 
 for(i=0;i<ei;i++)
  printf("%d ",edge[i]);
 printf("\n");
 dfs(0);
 for(i=0;i<ei;i++)
  printf("%d ",tree[i]);
 printf("\n");
// ran=rand();
  // l=(int)((((int)ran)%verc));
  // printf("%d\n",l);

 return 0;
}