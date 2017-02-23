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
			if(check[l])
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
			if(check[l])
				continue;
			if(!d_dcvslist[k][l])
				continue;
			check[l]=true;
			res|=chechall(ver,check,i+1,dfrom,dto,d_delist,d_qelist,d_qvid,qfrom,qto,d_dcvslist);
			if(res==true)
				return true;
			check[l]=false;
		}
	}
	return false;
}
__global__ void findcvs(int ver,int *d_dvert,int *d_dverc,int *d_delist,int *d_qvert,int *d_qelist,int *d_qvid,int ** d_dcvslist )
{
	//int i;
	int dver=threadIdx.x*blockDim.y+threadIdx.y+blockDim.x*blockDim.y*(blockIdx.x*gridDim.y+blockIdx.y);
	if(dver>=d_dverc[0])
		return ;
	int ql=d_qvert[ver+1]-d_qvert[ver];
	int dl=d_dvert[dver+1]-d_dvert[dver];
	if(ql>dl)
		return;
	bool *checked=(bool*)malloc(sizeof(bool)*d_dverc[0]);
	//bool *checked=new bool[d_dverc[0]];
	memset(checked,false,sizeof(bool)*d_dverc[0]);
	//chechall(bool *check,int i,int dfrom,int dto,int *d_delist,int *d_qelist,int *d_qvid,int qfrom,int qto,bool ** d_dcvslist)
	if(chechall(ver,checked,1,d_dvert[dver],d_dvert[dver+1],d_delist,d_qelist,d_qvid,d_qvert[ver],d_qvert[ver+1],d_dcvslist))
		d_dcvslist[d_qvid[ver]][dver]=true;
	free(checked);
}
__global__ void puttolist(int *d_dverc,int *d_loc,int * d_dcvslist )
{
	int dver=threadIdx.x*blockDim.y+threadIdx.y+blockDim.x*blockDim.y*(blockIdx.x*gridDim.y+blockIdx.y);
	if(dver>=d_dverc[0])
		return ;
	if(d_loc[dver]!=d_loc[dver+1])
		d_dcvslist[d_loc[dver]]=dver;
}
__global__ void checkperm(bool *found,int * qdmap,int * d_qverc,int * d_qelist,int * d_qvert,int * d_dvert,int *d_delist){
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
			p=d_dvert[dver];
			k=qdmap[j];
			for(;p<n;p++){
				if(k==d_delist[p]){
					flag=1;
					break;
				}
			}
			if(!flag){
				*found=false;
				return;
			}
		}
	//}
}

__global__ void findall(int *d_mapans,int *d_cans,int *d_qvid,int * d_qverc,int * d_qelist,int * d_qvert,int * d_dvert,int *d_delist,int ** d_cvsverlist,int * d_size_cvs)
{
	bool found[1]={true};
	long long int blockId = blockIdx.x 
			 + blockIdx.y * gridDim.x 
			 + gridDim.x * gridDim.y * blockIdx.z; 
	long long int threadId = blockId * (blockDim.x * blockDim.y * blockDim.z)
			  + (threadIdx.z * (blockDim.x * blockDim.y))
			  + (threadIdx.y * blockDim.x)
			  + threadIdx.x;
	int i=0;
	long long int indexperm=threadId;
	for(i=0;i<d_qverc[0];i++){
		int j=d_qvid[i];
		indexperm/=d_size_cvs[j];
	}
	if(indexperm)
		return;
	indexperm=threadId;	
	int *d_qdmap=&d_mapans[d_qverc[0]*threadId];
	for(i=0;i<d_qverc[0];i++){
		int j=d_qvid[i];
		d_qdmap[i]=d_cvsverlist[j][indexperm%d_size_cvs[j]];
		indexperm/=d_size_cvs[j];
	}
	
	//dim3 blocks((max/16 )+ 1,(max/16)+1);
	//dim3 threads(16,16);
	//found[0]=true;
	//checkperm<<<blocks,threads>>> (found,d_qdmap,d_qverc,d_qelist,d_qvert,d_dvert,d_delist);
	int n,p,j,k,flag=0,ver;
	for(ver=0;ver<d_qverc[0];ver++){
		int l=d_qvert[ver+1];
		int dver=d_qdmap[ver];
		
		n=d_dvert[dver+1];
		for(i=d_qvert[ver];i<l;i++)
		{
			flag=0;
			j=d_qelist[i];
			p=d_dvert[dver];
			k=d_qdmap[j];
			for(;p<n;p++){
				if(k==d_delist[p]){
					flag=1;
					break;
				}
			}
			if(!flag){
				*found=false;
				return;
			}
		}
	}
	if(found[0]){
		//printf("%d ", threadId);
	}
}
int * qdmap;
int *d_qverc,*d_dverc;
int *d_qvid,*d_qidtov,*h_qidtov,*h_qvid;
int *d_qvert,*d_qelist,*d_dvert,*d_delist;//,*d_dvelist,*d_qvelist;
bool *d_qtree,*d_over;	
int *d_qdmap;
bool h_over;
/*void callforallperm(bool * check,int ** cvslist,int i,int max,int dmax){
	int j,k,l;
	l=h_qvid[i-1];
	//printf("i%d %di",i,l);
	if(i==max){
		for(j=0;j<dmax;j++)
			if(cvslist[l][j] && !check[j]){
				qdmap[i-1]=j;
				dim3 blocks((max/16 )+ 1,(max/16)+1);
				dim3 threads(16,16);
				h_over=true;
				//for(k=0;k<max;k++)
				//		printf("%d ",qdmap[k]);
				
				cudaMemcpy(d_over, &h_over, sizeof(bool), cudaMemcpyHostToDevice) ;
				cudaMemcpy(d_qdmap, qdmap, sizeof(int)*(max+1), cudaMemcpyHostToDevice);
				checkperm<<<blocks,threads>>> (d_over,d_qdmap,d_qverc,d_qelist,d_qvert,d_dvert,d_delist);
				//checkperm(bool *found,int * qdmap,int * d_qverc,int * d_qelist,int * d_qvert,int * d_dvert,int *d_delist)
				cudaError_t err = cudaGetLastError();
					if(err!=cudaSuccess)
					{
						printf("Error: %s\n", cudaGetErrorString(err));
						printf("Not Ok");
					}
				cudaMemcpy(&h_over, d_over, sizeof(bool), cudaMemcpyDeviceToHost) ;
				if(h_over){
					for(k=0;k<max;k++)
						printf("%d ",qdmap[k]);
					//printf("\n");
						printf("OK\n");
				}
				//printf("\n");
			}
	}
	else{
		for(j=0;j<dmax;j++){
				//printf("%d %d %d\n",j,check[j],cvslist[l][j]);
				if(cvslist[l][j] && !check[j]){
				check[j]=true;
				qdmap[i-1]=j;
				callforallperm(check,cvslist,i+1,max,dmax);
				check[j]=false;
			}
		}
	}
}*/

int main(int argc, char **argv)
{
	int deviceId = 4;
	cudaSetDevice(deviceId);
	int h_qverc,h_dverc;
	
	
	int *h_qvert,*h_qelist,*h_dvert,*h_delist;//,*h_dvelist,*h_qvelist;
	
	bool *h_qtree;
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
	dim3 blocks((h_qverc/16 )+ 1,(h_qverc/16)+1);
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
	dim3 dblocks((h_dverc/16 )+ 1,(h_dverc/16)+1);
	dim3 dthreads(16,16);
	int **d_cvsverlist,**d_temverlist;
	int *d_size_cvs,*h_size_cvs;
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
	for(i=0;i<h_dverc;i++)
		h_cvslist[1][i]=i;
	cudaMemcpy(d_temverlist[1],h_cvslist[1],sizeof(int)*(h_dverc+1),cudaMemcpyHostToDevice);
	h_size_cvs[1]=h_dverc;
	for(i=0;i<=h_qverc;i++)
	{
		if(h_qidtov[i]!=-1)
		{

		//findcvs(int ver,int *d_dvert,int *d_dverc,int *d_delist,int *d_qvert,int *d_qelist,int *d_qvid,bool ** d_dcvslist )
		findcvs<<<dblocks,dthreads>>>(h_qidtov[i],d_dvert,d_dverc,d_delist,d_qvert,d_qelist,d_qvid,d_cvslist);
		cudaError_t err = cudaGetLastError();
	
		cudaMemcpy(h_cvslist[i],h_tem[i],sizeof(int)*(h_dverc+1),cudaMemcpyDeviceToHost);
		for(j=0;j<=h_dverc;j++)
			if(h_cvslist[i][j])
				printf("%d ",j);
		printf("\n");
		//printf("%d ",h_qidtov[i]);
		thrust::exclusive_scan(h_cvslist[i],h_cvslist[i]+h_dverc+1,h_cvslist[i]);
		h_size_cvs[i]=h_cvslist[i][h_dverc];
		cudaMemcpy(h_tem[i],h_cvslist[i],sizeof(int)*(h_dverc+1),cudaMemcpyHostToDevice);
		puttolist<<<dblocks,dthreads>>>(d_dverc,h_tem[i],d_temverlist[i]);
	//	cudaMemcpy(h_cvslist[i],d_temverlist[i],sizeof(int)*(h_dverc+1),cudaMemcpyDeviceToHost);
	//	for(j=0;j<=h_dverc;j++)
	//		printf("%d ",h_cvslist[i][j]);	
		}	
		
	}
	//	cudaMemcpy(h_delist,d_delist,sizeof(int)*(h_dvert[h_dverc]),cudaMemcpyDeviceToHost);
	//	for(j=0;j<h_dvert[h_dverc];j++)
	//		printf("%d ",h_delist[j]);	
			
	for(i=0;i<h_qverc;i++)
		if(h_size_cvs[h_qvid[i]])
		totalthreads*=h_size_cvs[h_qvid[i]];
	printf("Start %lld\n",totalthreads);
	cudaMemcpy(d_size_cvs,h_size_cvs,sizeof(int)*(h_qverc+1),cudaMemcpyHostToDevice);
	//totalthreads=1000;
	dim3 dpblocks(((int)(totalthreads/8 )+ 1),((int)(totalthreads/8)+1),((int)(totalthreads/8)+1));
	dim3 dpthreads(8,8,8);
	 int *d_mapans,*d_countans,*h_countans;
	cudaMalloc(&d_mapans,sizeof( int)*totalthreads*(h_qverc+1));
	cudaMalloc(&d_countans,sizeof( int)*(totalthreads+1));
	cudaMemset(d_countans,0,sizeof( int)*(totalthreads+1));
	h_countans=(int *)malloc(sizeof(int)*(totalthreads+1));
	//h_countans=0;
	//cudaMemcpy(d_countans, &h_countans, sizeof( int), cudaMemcpyHostToDevice) ;
	//cudaMemcpy(&h_countans, d_qverc, sizeof(int), cudaMemcpyDeviceToHost) ;
	//printf("%d\n",h_countans);
	findall<<<dpblocks,dpthreads>>> (d_mapans,d_countans,d_qvid,d_qverc,d_qelist,d_qvert,d_dvert,d_delist,d_cvsverlist,d_size_cvs);
	cudaMemcpy(h_countans, d_countans, sizeof(int)*(totalthreads+1), cudaMemcpyDeviceToHost) ;
	thrust::exclusive_scan(h_countans,h_countans+totalthreads+1,h_countans);
	printf("%d\n",h_countans[totalthreads] );
	//printf("%d\n",h_countans);
	/*j=0;
	for(i=0;i<totalthreads;i++)
		if(h_countans[i])
			j++;
	printf("%d\n",j);
	*///bool * check=(bool *)malloc(sizeof(bool)*(h_dverc+1));
	//memset(check,false,sizeof(bool)*(h_dverc+1));
	//qdmap=(int *)malloc(sizeof(int)*(h_qverc+1));
	//cudaMalloc(&d_qdmap,sizeof(int)*(h_qverc+1));
	//callforallperm(check,h_cvslist,1,h_qverc,h_dverc);

	cudaFree(d_over);
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
	/*free(h_qvid);
	free(h_qvert);
	//free(h_qelist);
	free(h_qidtov);
	free(h_cvslist);
	free(h_dvert);
	free(h_delist);*/
}
	
