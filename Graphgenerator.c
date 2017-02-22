#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <bits/stdc++.h>
int *ver,*edge,*count,*tree;
bool * vis;
int verc,avgoutdeg,i,j,k,l,m,ei=0;
void dfs(int i){
  int j;
  vis[i]=true;
  for(j=ver[i];j<ver[i+1];j++){
    if(!vis[edge[j]]){
      dfs(edge[j]);
      tree[j]=1;
    }
  }
}
int main()
{
 float ran,var,maxvar,curravg=0,currvar=0;
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
 var*=verc;
 edge=(int *)malloc(sizeof(int)*(verc*(avgoutdeg+var)));
 tree=(int *)malloc(sizeof(int)*(verc*(avgoutdeg+var)));
 memset(tree,0,sizeof(int)*(verc*(avgoutdeg+var)));
 maxvar=floor(sqrt(var));
 printf("%d\n",verc);
 for(i=0;i<verc;i++)
 {
  ran=(int)rand();
  //printf("%f",ran);
  if(maxvar)
   j=avgoutdeg-maxvar+(int)ran%((int)(maxvar*2+1));
  else
   j=avgoutdeg;
  ver[i]=curravg;
  curravg+=j;
  currvar+=(j-avgoutdeg)*(j-avgoutdeg);
 // printf("%d\n",j);
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