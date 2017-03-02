#include <stdio.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <stdlib.h>
#include "device_launch_parameters.h"
#include <thrust/scan.h>
const int BASE1 = 10000 + 7;
const int BASE2 = 100000 + 3;
const int MOD1 = 1000000 + 3;
const int MOD2 = 1000000 + 37;
__global__ void findhash(int *d_qvert,int *d_qverc,int *d_qvid,int *d_qelist,bool *d_over,bool *d_qtree,int *d_hash1,int *d_hash2)
{
	int i;
	int ver=threadIdx.x*blockDim.y+threadIdx.y+blockDim.x*blockDim.y*(blockIdx.x*gridDim.y+blockIdx.y);
	if(ver>=d_qverc[0])
		return ;
	if(d_qvid[ver]!=0)
		return;
	int l=d_qvert[ver+1];
	int hash1=1,hash2=1;
	int flag=0;
	for(i=d_qvert[ver];i<l;i++)
	{
		int m=d_qelist[i];
		bool treeedge=d_qtree[i];
		if(treeedge){
			int tt=d_qvid[m];
			if(tt==0)
				return;
			flag=1;
			hash1=(hash1*1L*BASE1)*tt % MOD1;
			hash2=(hash2*1L*BASE2)*tt % MOD2;
		}
	}
	if(flag==0)
		return;
	if(flag==1){
		*d_over=false;
		d_hash1[hash1]=1;
		d_hash2[hash2]=1;
	}
	
}
__global__ void setdeg1(int *d_qvert,int *d_qverc,int *d_qvid,bool *d_qtree)
{
	int i;
	int ver=threadIdx.x*blockDim.y+threadIdx.y+blockDim.x*blockDim.y*(blockIdx.x*gridDim.y+blockIdx.y);
	if(ver>=d_qverc[0])
		return ;
	if(d_qvid[ver]!=0)
		return;
	
	int l=d_qvert[ver+1];
	bool treeedge;
	
	for(i=d_qvert[ver];i<l;i++)
	{
		treeedge=d_qtree[i];
		
		if(treeedge)
			return;
	}
//printf("%d %d\n",ver,i);
	d_qvid[ver]=1;	
}
/*__global__ void alignhash(bool *d_hash1,bool *d_hash2)
{
	int ver=threadIdx.x*blockDim.y+threadIdx.y+blockDim.x*blockDim.y*(blockIdx.x*gridDim.y+blockIdx.y);
	if(ver>=1000038)
		return ;
	if(d_hash1[ver] || d_hash2[ver]){
		d_hash1=true;
	}
}*/
__global__ void puttoid(int *d_qvert,int *d_qverc,int *d_qvid,int *d_qelist,bool *d_qtree,int *d_loc,int * d_qidtov)
{
	int i;
	int ver=threadIdx.x*blockDim.y+threadIdx.y+blockDim.x*blockDim.y*(blockIdx.x*gridDim.y+blockIdx.y);
	if(ver>=d_qverc[0])
		return ;
	if(d_qvid[ver]!=0)
		return;
	int l=d_qvert[ver+1];
	int hash1=1,hash2=1;
	int flag=0;
	for(i=d_qvert[ver];i<l;i++)
	{
		int m=d_qelist[i];
		bool treeedge=d_qtree[i];
		if(treeedge){
			int tt=d_qvid[m];
			if(tt==0)
				return;
			flag=1;
			hash1=(hash1*1L*BASE1)*tt % MOD1;
			hash2=(hash2*1L*BASE2)*tt % MOD2;
		}
	}
	//printf("%d %d %d \n",ver,flag,d_loc[hash1]);
	if(flag==0)
		return;
	int id=d_loc[hash1];
	d_qvid[ver]=id;
	d_qidtov[id]=ver;
}
__device__ bool chechall(int ver,bool *check,int i,int dfrom,int dto,int *d_delist,int *d_qelist,int *d_qvid,int qfrom,int qto,int ** d_dcvslist){
	//int ql=qfrom-qto;
	int ql=qto-qfrom;
	int j,k,l;
	//d_dcvslist[2][ql]=true;
	if(i==ql){
		
		k=d_qelist[i+qfrom-1];
		k=d_qvid[k];
		if(k>=d_qvid[ver])
			return true;
		for(j=dfrom;j<dto;j++){
			l=d_delist[j];
			if(check[j])
				continue;
			if(!d_dcvslist[k][l])
				continue;
			return true;
		}
	}
	else{
		int res=false;
		k=d_qelist[i+qfrom-1];
		k=d_qvid[k];
		if(k>=d_qvid[ver])
			return chechall(ver,check,i+1,dfrom,dto,d_delist,d_qelist,d_qvid,qfrom,qto,d_dcvslist);
		
		for(j=dfrom;j<dto;j++){
			l=d_delist[j];
			if(check[j])
				continue;
			if(!d_dcvslist[k][l])
				continue;
			check[j]=true;
			res|=chechall(ver,check,i+1,dfrom,dto,d_delist,d_qelist,d_qvid,qfrom,qto,d_dcvslist);
			if(res==true)
				return true;
			check[j]=false;
		}
	}
	return false;
}
__global__ void findcvs(bool *temp,int ver,int *d_dvert,int *d_dverc,int *d_delist,int *d_qvert,int *d_qelist,int *d_qvid,int ** d_dcvslist )
{
	//int i;
	int dver=threadIdx.x*blockDim.y+threadIdx.y+blockDim.x*blockDim.y*(blockIdx.x*gridDim.y+blockIdx.y);
	if(dver>=d_dverc[0])
		return ;
	int ql=d_qvert[ver+1]-d_qvert[ver];
	int dl=d_dvert[dver+1]-d_dvert[dver];
	if(ql>dl)
		return;
	//bool *checked=(bool*)malloc(sizeof(bool)*d_dverc[0]);
	//bool *checked=new bool[d_dverc[0]];
	//memset(checked,false,sizeof(bool)*d_dverc[0]);
	//chechall(bool *check,int i,int dfrom,int dto,int *d_delist,int *d_qelist,int *d_qvid,int qfrom,int qto,bool ** d_dcvslist)
	if(chechall(ver,temp,1,d_dvert[dver],d_dvert[dver+1],d_delist,d_qelist,d_qvid,d_qvert[ver],d_qvert[ver+1],d_dcvslist))
		d_dcvslist[d_qvid[ver]][dver]=true;
	//free(checked);
}
__global__ void puttolist(int *d_dverc,int *d_loc,int * d_dcvslist )
{
	int dver=threadIdx.x*blockDim.y+threadIdx.y+blockDim.x*blockDim.y*(blockIdx.x*gridDim.y+blockIdx.y);
	if(dver>=d_dverc[0])
		return ;
	if(d_loc[dver]!=d_loc[dver+1])
		d_dcvslist[d_loc[dver]]=dver;
}
__global__ void checkperm(int *found,int * qdmap,int * d_qverc,int * d_qelist,int * d_qvert,int * d_dvert,int *d_delist,bool *d_qtree){
	int i;
	//found[0]=false;
	int ver=threadIdx.x*blockDim.y+threadIdx.y+blockDim.x*blockDim.y*(blockIdx.x*gridDim.y+blockIdx.y);
	if(ver>=d_qverc[0])
		return ;
	int n,p,j,k,flag=0;
	//for(ver=0;ver<d_qverc[0];ver++){
		int l=d_qvert[ver+1];
		int dver=qdmap[ver];
		
		n=d_dvert[dver+1];
		for(i=d_qvert[ver];i<l;i++)
		{
			flag=0;
			j=d_qelist[i];
			//if(!d_qtree[i])
			//	continue;
			p=d_dvert[dver];
			k=qdmap[j];
			for(;p<n;p++){
				if(k==d_delist[p]){
					flag=1;
					break;
				}
			}
			if(!flag){
				//*found=false;
				if(d_qtree[i]){				
					found[0]=found[1]=-1;		
					return;
				}
				else
					found[1]=0,found[0]++;
			}
		}
	//}
}

int * qdmap;
int *d_qverc,*d_dverc;
int *d_qvid,*d_qidtov,*h_qidtov,*h_qvid;
int *d_qvert,*d_qelist,*d_dvert,*d_delist;//,*d_dvelist,*d_qvelist;
bool *d_qtree,*d_over;	
int *d_qdmap;
bool h_over;
bool *h_qtree;
int *d_size_cvs,*h_size_cvs,ansi=0,treeansi=0;
long long int * h_anslist,*d_anslist;
long long int * h_treeanslist,*d_treeanslist;
int *h_treeremlist,*d_treeremlist;
__global__ void processoperation(int type,int a,int b,int nans,long long int *anslist,int *qverc,int *dverc,int *qvert,int *qelist,int *dvert,int *delist,int **cvsverlist,int *size_cvs,int *qvid){
	int ansi=threadIdx.x*blockDim.y+threadIdx.y+blockDim.x*blockDim.y*(blockIdx.x*gridDim.y+blockIdx.y);
	if(ansi>=nans)
		return ;
	long long int indexperm=anslist[ansi];
	if(indexperm==-1)
		return;
	//int *d_qdmap=new int[d_qverc[0]];//&d_mapans[d_qverc[0]*threadId];//new int[d_qverc[0]];
        int mapvera,mapverb=-1;
	int *aedges=NULL,till,i;
	//printf("%d\n ",indexperm);
	//anslist[ansi]=-1;
	for(i=0;i<qverc[0];i++){
                int j=qvid[i];
//                printf("j%d %d %d %dj ",i,j,size_cvs[j],mapvera);
		if(type==0 && i==a)
			mapvera=cvsverlist[j][indexperm%size_cvs[j]],aedges=&delist[dvert[mapvera]],till=dvert[mapvera+1];
                else if(type==0 && i==b)
			mapverb=cvsverlist[j][indexperm%size_cvs[j]];
		else if(type==1 && cvsverlist[j][indexperm%size_cvs[j]]==a)
			mapvera=i,aedges=&qelist[qvert[mapvera]],till=qvert[mapvera+1];
		else if(type==1 && cvsverlist[j][indexperm%size_cvs[j]]==b)
			mapverb=i;
                indexperm/=size_cvs[j];
        }
	bool flag=false;
//		anslist[ansi]=-1;
	if(aedges==NULL || mapverb==-1 || indexperm>0)
		return;
//	printf("j%d %dj",aedges[0],mapverb);
	
	for(i=0;i<till;i++){
		if(aedges[i]==mapverb){
			if(type==1){
				anslist[ansi]=-1;
				break;
			}
			flag=true;
			break;
		}
	}
	if(!flag && type==0)
		anslist[ansi]=-1;
}

__global__ void processqdnontree(int type,int a,int b,int ntans,long long int *tanslist,int *tremlist,int *qverc,int *dverc,int *qvert,int *qelist,int *dvert,int *delist,int **cvsverlist,int *size_cvs,int *qvid){
	int ansi=threadIdx.x*blockDim.y+threadIdx.y+blockDim.x*blockDim.y*(blockIdx.x*gridDim.y+blockIdx.y);
	if(ansi>=ntans)
		return ;
	long long int indexperm=tanslist[ansi];
	if(tremlist[ansi]==0)
		return;
        int mapvera,mapverb=-1;
	int *aedges=NULL,till,i;
	for(i=0;i<qverc[0];i++){
                int j=qvid[i];
		if(i==a)
			mapvera=cvsverlist[j][indexperm%size_cvs[j]],aedges=&delist[dvert[mapvera]],till=dvert[mapvera+1];
		else if(i==b)
			mapverb=cvsverlist[j][indexperm%size_cvs[j]];
                indexperm/=size_cvs[j];
        }
	
	bool flag=false;
//		anslist[ansi]=-1;
	if(aedges==NULL || mapverb==-1 || indexperm>0)
		return;
//	printf("j%d %dj",aedges[0],mapverb);
	for(i=0;i<till;i++){
		if(aedges[i]==mapverb){
			flag=true;
			break;
		}
	}
	if(!flag)
		 atomicDec((unsigned int *)&tremlist[ansi],tremlist[ansi]);
}

__global__ void doquery(int *nquery,int * type,int * vera,int * verb,int *ntans,long long int * tanslist,int *tremlist,int *nans,long long int * anslist,int ** cvsmatrix,int **cvslist,int *qverc,int *dverc,int *qvert,int *qelist,int *dvert,int *delist,int *size_cvs,int *qvid,bool *qtree){
	
	int qi=threadIdx.x*blockDim.y+threadIdx.y+blockDim.x*blockDim.y*(blockIdx.x*gridDim.y+blockIdx.y);
	if(qi>=nquery[0])
		return ;
	//nans[0]=1;
	dim3 blocks((sqrtf(nans[0])/16 )+ 1,(sqrtf(nans[0])/16)+1);
	dim3 threads(16,16);
	if(type[qi]==0){
		processoperation<<<blocks,threads>>>(type[qi],vera[qi],verb[qi],nans[0],anslist,qverc,dverc,qvert,qelist,dvert,delist,cvslist,size_cvs,qvid);
	}
	else if(type[qi]==1){
		processoperation<<<blocks,threads>>>(type[qi],vera[qi],verb[qi],nans[0],anslist,qverc,dverc,qvert,qelist,dvert,delist,cvslist,size_cvs,qvid);
	}
	else if(type[qi]==3){
		int l=qvert[vera[qi]+1],i;
		bool flag=false,istree=false;
		for(i=qvert[vera[qi]];i<l;i++){
			if(qelist[i]==verb[qi]){
				flag=true;
				if(qtree[i])
					istree=true;
				break;
			}
		}
		if(!flag)
			return;
		if(!istree){
			processqdnontree<<<blocks,threads>>>(type[qi],vera[qi],verb[qi],ntans[0],tanslist,tremlist,qverc,dverc,qvert,qelist,dvert,delist,cvslist,size_cvs,qvid);
		}
		else{
			//processqdtree<<<blocks,threads>>>(type[qi],vera[qi],verb[qi],ntans[0],tanslist,tremlist,qverc,dverc,qvert,qelist,dvert,delist,cvslist,size_cvs,qvid);
		}
	}
	else{
	}

}

void callforallperm(bool * check,int ** cvslist,int i,int max,int dmax,long long int mapid){
	int j,k,l;
	l=h_qvid[i-1];
	//printf("i%d %di",i,l);
	if(i==max){
		for(j=0;j<dmax;j++)
			if(cvslist[l][j] && !check[j]){
				qdmap[i-1]=j;
				//mapid+=j*h_size_cvs[l];
				dim3 blocks((sqrt(max)/16 )+ 1,(sqrt(max)/16)+1);
				dim3 threads(16,16);
				int h_found[2]={0,1};
				int *d_found;
				cudaMalloc(&d_found,sizeof(int)*2);
				cudaMemcpy(d_found, h_found, sizeof(int)*2, cudaMemcpyHostToDevice) ;
				cudaMemcpy(d_qdmap, qdmap, sizeof(int)*(max+1), cudaMemcpyHostToDevice);
				checkperm<<<blocks,threads>>> (d_found,d_qdmap,d_qverc,d_qelist,d_qvert,d_dvert,d_delist,d_qtree);
				//checkperm(bool *found,int * qdmap,int * d_qverc,int * d_qelist,int * d_qvert,int * d_dvert,int *d_delist)
				cudaError_t err = cudaGetLastError();
					if(err!=cudaSuccess)
					{
						printf("Error: %s\n", cudaGetErrorString(err));
						printf("Not Ok");
					}
				cudaMemcpy(h_found, d_found, sizeof(int)*2, cudaMemcpyDeviceToHost) ;
				if(h_found[1]==1){
					//for(k=0;k<max;k++)
					//	printf("%d ",qdmap[k]);
					//printf("\n");
					//	printf("OK\n");
					h_anslist[ansi++]=mapid*h_size_cvs[l]+j;
				}
				else if(h_found[0]!=-1){
					h_treeanslist[treeansi++]=mapid*h_size_cvs[l]+j;
					h_treeremlist[treeansi++]=h_found[0]/2;
				}
				//printf("\n");
			}
	}
	else{
		for(j=0;j<dmax;j++){
				//printf("%d %d %d\n",j,check[j],cvslist[l][j]);
				if(cvslist[l][j] && !check[j]){
				//ansi++;
				check[j]=true;
				qdmap[i-1]=j;
				//mapid+=l*h_size_cvs[l];
				callforallperm(check,cvslist,i+1,max,dmax,mapid*h_size_cvs[l] +j );
				check[j]=false;
			}
		}
	}
}

int main(int argc, char **argv)
{
	int deviceId = 4;
	cudaSetDevice(deviceId);
	int h_qverc,h_dverc;
	
	
	int *h_qvert,*h_qelist,*h_dvert,*h_delist;//,*h_dvelist,*h_qvelist;
	
	int *d_hash1,*d_hash2;
	
	int i,j;
	int **h_cvslist,**d_cvslist,**h_tem;
	scanf("%d",&h_qverc);
	h_qvert=(int *)malloc(sizeof(int)*(h_qverc+1));
	h_qvid=(int *)malloc(sizeof(int)*(h_qverc+1));
	h_qidtov=(int *)malloc(sizeof(int)*(h_qverc+1));
	h_tem=(int **)malloc(sizeof(int*)*(h_qverc+1));
	h_cvslist=(int **)malloc(sizeof(int*)*(h_qverc+1));
	for(i=0;i<=h_qverc;i++){
		scanf("%d",&h_qvert[i]);
	}
	h_qelist=(int *)malloc(sizeof(int)*h_qvert[h_qverc]);
	for(i=0;i<h_qvert[h_qverc];i++)
		scanf("%d",&h_qelist[i]);
	h_qtree=(bool *)malloc(sizeof(bool)*h_qvert[h_qverc]);
	for(i=0;i<h_qvert[h_qverc];i++){
		scanf("%d",&j);
		if(j==1)
			h_qtree[i]=true;
		else
			h_qtree[i]=false;
	}
	
	scanf("%d",&h_dverc);
	h_dvert=(int *)malloc(sizeof(int)*(h_dverc+1));
	for(i=0;i<=h_dverc;i++){
		scanf("%d",&h_dvert[i]);
	}
	for(i=0;i<=h_qverc;i++)
		h_cvslist[i]=(int *)malloc(sizeof(int)*(h_dverc+1));
	
	h_delist=(int *)malloc(sizeof(int)*h_dvert[h_dverc]);
	for(i=0;i<h_dvert[h_dverc];i++)
		scanf("%d",&h_delist[i]);
	cudaMalloc(&d_qverc,sizeof(int));
	cudaMalloc(&d_over,sizeof(bool));
	cudaMalloc(&d_qvert,sizeof(int)*(h_qverc+1));
	cudaMalloc(&d_qidtov,sizeof(int)*(h_qverc+1));
	//cudaMalloc(&d_loc,sizeof(int)*(h_qverc+1));
	cudaMalloc(&d_qelist,sizeof(int)*h_qvert[h_qverc]);
	cudaMalloc(&d_qtree,sizeof(bool)*h_qvert[h_qverc]);
	cudaMalloc(&d_hash1,sizeof(int)*1000038);
	cudaMalloc(&d_hash2,sizeof(int)*1000038);
	
	cudaMalloc(&d_qvid,sizeof(int)*(h_qverc+1));
	cudaMemcpy(d_qverc,&h_qverc,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d_qvert,h_qvert,sizeof(int)*(h_qverc+1),cudaMemcpyHostToDevice);
	cudaMemcpy(d_qelist,h_qelist,sizeof(int)*h_qvert[h_qverc],cudaMemcpyHostToDevice);
	cudaMemcpy(d_qtree,h_qtree,sizeof(bool)*h_qvert[h_qverc],cudaMemcpyHostToDevice);
	cudaMemset(d_hash1,0,sizeof(int)*1000038);
	cudaMemset(d_hash2,0,sizeof(int)*1000038);
	//cudaMemset(d_loc,0,sizeof(int)*(h_qverc+1));
	cudaMemset(d_qidtov,-1,sizeof(int)*(h_qverc+1));
	cudaMemset(d_qvid,0,sizeof(int)*(h_qverc+1));
	int *h_hash1=(int *)malloc(sizeof(int)*1000038);
	int *h_hash2=(int *)malloc(sizeof(int)*1000038);
	dim3 blocks((sqrt(h_qverc)/16 )+ 1,(sqrt(h_qverc)/16)+1);
	dim3 threads(16,16);
	
	//int *d_qvert,int *d_dverc,int *d_qvid,int *d_qelist,bool *d_over,bool *d_hash1,bool *d_hash2)
	h_over=true;
	//h_qvid[1]=1;
	//h_qvid[3]=1;
	//cudaMemcpy(d_qvid,h_qvid,sizeof(int)*(h_qverc+1),cudaMemcpyHostToDevice);
//printf("qt%d %dqt\n",h_qtree[0],h_qtree[1]);
	setdeg1<<<blocks,threads>>>(d_qvert,d_qverc,d_qvid,d_qtree);

	h_over=false;
	int maxval=2;
	while(!h_over)
	{
		h_over=true;
		cudaMemcpy(d_over, &h_over, sizeof(bool), cudaMemcpyHostToDevice) ;
		cudaMemset(d_hash1,0,sizeof(int)*1000038);
		findhash <<<blocks,threads>>> (d_qvert,d_qverc,d_qvid,d_qelist,d_over,d_qtree,d_hash1,d_hash2);
		//(int *d_qvert,int *d_dverc,int *d_qvid,int *d_qelist,bool *d_over,bool *d_hash1,bool *d_qtree,bool *d_hash2)
		cudaError_t err = cudaGetLastError();
		if(err!=cudaSuccess)
		{
			printf("Error: %s\n", cudaGetErrorString(err));
			printf("Not Ok");
		}
		cudaMemcpy(h_hash1,d_hash1,sizeof(int)*1000038,cudaMemcpyDeviceToHost);
		h_hash1[0]+=maxval;
		thrust::exclusive_scan(h_hash1,h_hash1+1000038,h_hash1);
		maxval=h_hash1[1000037];
		cudaMemcpy(d_hash1,h_hash1,sizeof(int)*1000038,cudaMemcpyHostToDevice);
		puttoid<<<blocks,threads>>>(d_qvert,d_qverc,d_qvid,d_qelist,d_qtree,d_hash1,d_qidtov);
///		cudaMemcpy(h_hash2,d_hash2,sizeof(bool)*1000038,cudaMemcpyDeviceToHost);
		cudaMemcpy(&h_over, d_over, sizeof(bool), cudaMemcpyDeviceToHost) ;
		
		//printf("over flag:%d ",h_over);
		/*for(i=0;i<h_qverc;i++){
			//if()
				printf("%d ",h_qvid[i]);
//			if(h_hash2[i])
//				printf("h2 %d ",i);
//			if(h_hash1[i] || h_hash2[i])
//				printf("\n");
		}
		printf("\n");*/

	}	
	cudaMemcpy(h_qvid,d_qvid,sizeof(int)*h_qverc,cudaMemcpyDeviceToHost);
	cudaMemcpy(h_qidtov,d_qidtov,sizeof(int)*(h_qverc+1),cudaMemcpyDeviceToHost);
	for(i=0;i<=h_qverc;i++){
		printf("%d ",h_qidtov[i]);
	}
	printf("\n");
	for(i=0;i<=h_qverc;i++){
		printf("%d ",h_qvid[i]);
	}
	printf("\n");
	cudaFree(d_qtree);
	cudaFree(d_hash1);
	cudaFree(d_hash2);
	free(h_hash1);
	free(h_hash2);
	free(h_qtree);
	
	cudaMalloc(&d_cvslist,sizeof(int*)*(h_qverc+1));
	for(i=0;i<=h_qverc;i++){
		cudaMalloc(&h_tem[i],sizeof(int)*(h_dverc+1));
		cudaMemset(h_tem[i],0,sizeof(int)*(h_dverc+1));
	}
	cudaMemset(h_tem[1],1,sizeof(int)*(h_dverc+1));
	cudaMemcpy(d_cvslist,h_tem,sizeof(int*)*(h_qverc+1),cudaMemcpyHostToDevice);
	cudaMalloc(&d_dvert,sizeof(int)*(h_dverc+1));
	cudaMalloc(&d_dverc,sizeof(int));
	cudaMalloc(&d_delist,sizeof(int)*h_dvert[h_dverc]);	
	cudaMemcpy(d_dverc,&h_dverc,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d_dvert,h_dvert,sizeof(int)*(h_dverc+1),cudaMemcpyHostToDevice);
	cudaMemcpy(d_delist,h_delist,sizeof(int)*h_dvert[h_dverc],cudaMemcpyHostToDevice);
	dim3 dblocks((sqrt(h_dverc)/16 )+ 1,(sqrt(h_dverc)/16)+1);
	dim3 dthreads(16,16);
	int **d_cvsverlist,**d_temverlist;
	memset(h_cvslist[1],1,sizeof(int)*(h_dverc+1));
	h_size_cvs=(int *)malloc(sizeof(int)*(h_qverc+1));
	memset(h_size_cvs,0,sizeof(int)*(h_qverc+1));
	cudaMalloc(&d_size_cvs,sizeof(int)*(h_qverc+1));
	cudaMemset(d_size_cvs,0,sizeof(int)*(h_qverc+1));
	cudaMalloc(&d_cvsverlist,sizeof(int*)*(h_qverc+1));
	d_temverlist=(int **)malloc(sizeof(int*)*(h_qverc+1));
	for(i=0;i<=h_qverc;i++){
		cudaMalloc(&d_temverlist[i],sizeof(int)*(h_dverc+1));
		cudaMemset(d_temverlist[i],0,sizeof(int)*(h_dverc+1));
	}
	cudaMemcpy(d_cvsverlist,d_temverlist,sizeof(int*)*(h_qverc+1),cudaMemcpyHostToDevice);
	long long int totalthreads=1;
	int *h_temploc;
	h_temploc=(int *)malloc(sizeof(int)*(h_dverc+1));
	for(i=0;i<h_dverc;i++)
		h_temploc[i]=i;
	cudaMemcpy(d_temverlist[1],h_temploc,sizeof(int)*(h_dverc+1),cudaMemcpyHostToDevice);
	h_size_cvs[1]=h_dverc;
	bool *d_tempcheck;	
	cudaMalloc(&d_tempcheck,sizeof(bool)*(h_dvert[h_dverc]+1));	
	printf("Starting cvs find\n");
	for(i=0;i<=h_qverc;i++)
	{
		if(h_qidtov[i]!=-1)
		{
		cudaMemset(d_tempcheck,false,sizeof(bool)*h_dvert[h_dverc]);	

		//findcvs(int ver,int *d_dvert,int *d_dverc,int *d_delist,int *d_qvert,int *d_qelist,int *d_qvid,bool ** d_dcvslist )
		findcvs<<<dblocks,dthreads>>>(d_tempcheck,h_qidtov[i],d_dvert,d_dverc,d_delist,d_qvert,d_qelist,d_qvid,d_cvslist);
		cudaError_t err = cudaGetLastError();
	
		printf("id %d \n",i);
		cudaMemcpy(h_cvslist[i],h_tem[i],sizeof(int)*(h_dverc+1),cudaMemcpyDeviceToHost);
		//printf("%d ",h_qidtov[i]);
		thrust::exclusive_scan(h_cvslist[i],h_cvslist[i]+h_dverc+1,h_temploc);
		h_size_cvs[i]=h_temploc[h_dverc];
		cudaMemcpy(h_tem[i],h_temploc,sizeof(int)*(h_dverc+1),cudaMemcpyHostToDevice);
		puttolist<<<dblocks,dthreads>>>(d_dverc,h_tem[i],d_temverlist[i]);
		//for(j=0;j<=h_dverc;j++)
		//	if(h_cvslist[i][j])
		//		printf("%d ",j);
		//printf("\n");
	//	cudaMemcpy(h_cvslist[i],d_temverlist[i],sizeof(int)*(h_dverc+1),cudaMemcpyDeviceToHost);
		//cudaMemcpy(h_temploc,h_tem[i],sizeof(int)*(h_size_cvs[i]),cudaMemcpyDeviceToHost);
		//printf("On list");
		//for(j=0;j<h_size_cvs[i];j++)
		//	printf("%d ",h_temploc[j]);
		//printf("\n");	
		printf("size %d\n",h_size_cvs[i]);
		}	
		
	}
	cudaMemcpy(d_size_cvs,h_size_cvs,sizeof(int)*(h_qverc+1),cudaMemcpyHostToDevice);
	//cudaMemcpy(h_delist,d_delist,sizeof(int)*(h_dvert[h_dverc]),cudaMemcpyDeviceToHost);
	//	for(j=0;j<h_dvert[h_dverc];j++)
	//		printf("%d ",h_delist[j]);	
	cudaFree(d_tempcheck);
	free(h_temploc);
	bool * check=(bool *)malloc(sizeof(bool)*(h_dverc+1));
	memset(check,false,sizeof(bool)*(h_dverc+1));
	qdmap=(int *)malloc(sizeof(int)*(h_qverc+1));
	cudaMalloc(&d_qdmap,sizeof(int)*(h_qverc+1));
	h_anslist=(long long int *)malloc(sizeof(long long int)*1000001);
	cudaMalloc(&d_anslist,sizeof(long long int)*(1000001));
	h_treeanslist=(long long int*)malloc(sizeof(long long int)*1000001);
	cudaMalloc(&d_treeanslist,sizeof(long long int)*(1000001));
	h_treeremlist=(int *)malloc(sizeof(int)*1000001);
	cudaMalloc(&d_treeremlist,sizeof(int)*(1000001));
	ansi=0;
	callforallperm(check,h_cvslist,1,h_qverc,h_dverc,0);
	
	printf("Final:%d\n",ansi);
	//answers found
	int * d_ansi,*d_treeansi;
	cudaMalloc(&d_ansi,sizeof(int));
	cudaMemcpy(d_ansi,&ansi,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d_anslist,h_anslist,sizeof(long long int)*(ansi),cudaMemcpyHostToDevice);
	cudaMalloc(&d_treeansi,sizeof(int));
	cudaMemcpy(d_treeansi,&treeansi,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d_treeanslist,h_treeanslist,sizeof(long long int)*(treeansi),cudaMemcpyHostToDevice);
	cudaMemcpy(d_treeremlist,h_treeremlist,sizeof(int)*(treeansi),cudaMemcpyHostToDevice);
	

/*	int nqueries,*d_nqueries;
	int *h_vera,*h_verb;
	int *d_vera,*d_verb;
	int *h_type,*d_type;
	scanf("%d",&nqueries);	
	h_vera=(int*) malloc(sizeof(int)*nqueries);
	h_verb=(int*) malloc(sizeof(int)*nqueries);
	h_type=(int *) malloc(sizeof(int)*nqueries);
	
		
	cudaMalloc(&d_vera,sizeof(int)*nqueries);
	cudaMalloc(&d_verb,sizeof(int)*nqueries);
	cudaMalloc(&d_type,sizeof(int)*nqueries);

	for(i=0;i<nqueries;i++){
		scanf("%d%d%d",&h_type[i],&h_vera[i],&h_verb[i]);
	}
	cudaMemcpy(d_vera,h_vera,sizeof(int)*nqueries,cudaMemcpyHostToDevice);
	cudaMemcpy(d_verb,h_verb,sizeof(int)*nqueries,cudaMemcpyHostToDevice);
	cudaMemcpy(d_type,h_type,sizeof(int)*nqueries,cudaMemcpyHostToDevice);
	
	cudaMalloc(&d_nqueries,sizeof(int));
	cudaMemcpy(d_nqueries,&nqueries,sizeof(int),cudaMemcpyHostToDevice);
	
	dim3 qblocks((sqrt(nqueries)/16 )+ 1,(sqrt(nqueries)/16)+1);
	dim3 qthreads(16,16);

	doquery<<<qblocks,qthreads>>>(d_nqueries,d_type,d_vera,d_verb,d_treeansi,d_treeanslist,d_treeremlist,d_ansi,d_anslist,d_cvslist,d_cvsverlist,d_qverc,d_dverc,d_qvert,d_qelist,d_dvert,d_delist,d_size_cvs,d_qvid,d_qtree);

	
	cudaMemcpy(h_anslist,d_anslist,sizeof(long long int)*(ansi),cudaMemcpyDeviceToHost);
	for(i=0;i<ansi;i++)
		if(h_anslist[i]==-1)
			printf(" %d ",i);
*/	cudaFree(d_over);
	cudaFree(d_qverc);
	cudaFree(d_qvert);
	cudaFree(d_qelist);
	cudaFree(d_qvid);
	cudaFree(d_qidtov);
	cudaFree(d_dvert);
	cudaFree(d_delist);
	cudaFree(d_dverc);
	cudaFree(d_cvslist);
	cudaFree(d_cvsverlist);
	cudaFree(d_size_cvs);
	cudaFree(d_anslist);
	/*free(h_qvid);
	free(h_qvert);
	//free(h_qelist);
	free(h_qidtov);
	free(h_cvslist);
	free(h_dvert);
	free(h_delist);*/
}	
