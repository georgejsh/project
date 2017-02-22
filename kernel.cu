#include <stdio.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <stdlib.h>
#include "device_launch_parameters.h"
#include <thrust/scan.h>
/*
__global__ void dfs_parallel(int *d_frontier1,int *d_vertex,int *d_loc,int *d_edge,int *d_frontier2)
{
	int i = blockIdx.x * gridDim.y + blockIdx.y,j=i*blockDim.x*blockDim.y+threadIdx.x*blockDim.y+threadIdx.y;
	if(d_frontier1[j]==-1)
		return;
	int ver=d_frontier1[j],k=d_loc[j],l=d_vertex[ver+1];
    for(i=d_vertex[ver];i<l;i++)
	{
		d_frontier2[k++]=d_edge[i];
	}
}

__global__ void cull(int *d_frontier1,int *d_dist,int *d_frontier2,int *d_count)
{
	int i = blockIdx.x * gridDim.y + blockIdx.y,j=i*blockDim.x*blockDim.y+threadIdx.x*blockDim.y+threadIdx.y;
	int ver=d_frontier1[j];
	if(d_dist[ver]!=-1)
		return;
	//atomicAdd(d_count,1);
}

__global__ void findloc(int *d_front,int *d_vertex,int *d_loc)
{
	int i = blockIdx.x * gridDim.y + blockIdx.y,j=i*blockDim.x*blockDim.y+threadIdx.x*blockDim.y+threadIdx.y;
	int ver=d_front[j];
	d_loc[j]=d_vertex[ver+1]-d_vertex[ver];
}
*/
__global__ void parallel(int *d_vertex,int *d_verc,int *d_edge,int *d_dist,bool *d_over,int *d_weight)
{
	int i;
	int ver=threadIdx.x*blockDim.y+threadIdx.y+blockDim.x*blockDim.y*(blockIdx.x*gridDim.y+blockIdx.y);
	//threadIdx.x*blockIdx.y+threadIdx.y;
	//d_dist[1]=3;
	//int ver=(threadIdx.x+blockIdx.x*blockDim.x)*gridDim.y*blockDim.y+threadIdx.y+blockIdx.y*blockDim.y;
	if(ver>=d_verc[0])
		return ;
	int k=d_dist[ver];
	if(k!=-1)
	{
		int l=d_vertex[ver+1];
		for(i=d_vertex[ver];i<l;i++)
		{
			int m=d_edge[i],p=d_dist[m],q=k+d_weight[i];
			if(p==-1)
			{
				d_dist[m]=q;d_over[0]=true;
			}
			else if(p>q)
			{
				d_dist[m]=q;d_over[0]=true;
			}

		}
	}

}


	void dfs(int *d_vertex,int *d_verc,int *d_edge,int *d_dist,int *d_source)
{
	
/*
	int * d_frontier1,* d_frontier2,*d_loc,d_frontcount,done=0;
	int *temp,numBytes=sizeof(int)*d_vertex[*d_verc];
	cudaMalloc(&d_frontier1,numBytes);
	cudaMalloc(&d_loc,numBytes);
	cudaMalloc(&d_frontier2,numBytes);

	d_frontier1[0]=*d_source;
	d_frontcount=1;


	while(d_frontcount)
	{
		dim3 blocks(d_frontcount/16,d_frontcount/16);
		dim3 threads(16, 16);
		findloc <<< blocks,threads >>> (d_frontier1,d_vertex,d_loc);
		thrust::exclusive_scan(d_loc,d_loc+d_frontcount,d_loc);

		dfs_parallel<<<blocks,threads>>>(d_frontier1,d_vertex,d_loc,d_edge,d_frontier2);
		cudaMemset(d_frontier1,-1,numBytes);
		int l=d_frontier1[d_frontcount-1],k=d_loc[d_frontcount-1]+d_vertex[l+1]-d_vertex[l];
		d_frontcount=0;
		
		dim3 blocks1(k/16,k/16);
		dim3 threads1(16, 16);		
		cull<<<blocks1,threads1>>>(d_frontier2,d_dist,d_frontier1,&d_frontcount);
	}


	*/



}
int main(int argc, char **argv)
{
	
//	GpuTimer timer2;
	int i,j,k,numBytes;
	int *d_vertex,*d_verc,*d_edge, *d_dist,*d_source,*d_weight;
	int *h_vertex, h_verc,*h_edge, *h_dist, h_source,h_ec,*h_final,*h_weight;
	scanf("%d%d",&h_verc,&h_ec);
	numBytes=sizeof(int)*(h_verc+1);
	h_vertex=(int *)malloc(numBytes);
	h_dist=(int *)malloc(numBytes);
	h_final=(int *)malloc(numBytes);
	memset(h_dist,-1,numBytes);
	h_edge=(int *)malloc(sizeof(int)*h_ec);
	h_weight=(int *)malloc(sizeof(int)*h_ec);
	for(i=0;i<=h_verc;i++)
		scanf("%d",&h_vertex[i]);
	for(i=0;i<h_ec;i++)
		scanf("%d",&h_edge[i]);
	for(i=0;i<h_ec;i++)
		scanf("%d",&h_weight[i]);
	scanf("%d",&h_source);
	h_dist[h_source]=0;

//	timer2.Start();

	cudaMalloc(&d_vertex,numBytes);
	cudaMalloc(&d_edge,sizeof(int)*h_ec);
	cudaMalloc(&d_weight,sizeof(int)*h_ec);
	cudaMalloc(&d_dist,numBytes);
	cudaMalloc(&d_verc,sizeof(int));
	cudaMalloc(&d_source,sizeof(int));
	cudaMemcpy(d_vertex,h_vertex,numBytes,cudaMemcpyHostToDevice);
	cudaMemcpy(d_edge,h_edge,sizeof(int)*h_ec,cudaMemcpyHostToDevice);
	cudaMemcpy(d_weight,h_weight,sizeof(int)*h_ec,cudaMemcpyHostToDevice);
	cudaMemcpy(d_verc,&h_verc,sizeof(int),cudaMemcpyHostToDevice);
	//cudaMemcpy(d_source,&h_source,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d_dist,h_dist,numBytes,cudaMemcpyHostToDevice);
	
//	timer2.Stop();
//	printf("(Memory Transfer)Time Elapsed :%lfms\n",timer2.Elapsed());
//	timer2.Start();
	//dfs(d_vertex,d_verc,d_edge,d_dist,d_source);

	
	
	
	bool stop=true;
	bool *d_over;
	dim3 blocks(h_verc/1024+1,h_verc/1024+1);
	dim3 threads(32,32);
	cudaMalloc(&d_over,sizeof(bool));
	while(stop)
	{
		stop=false;
		cudaMemcpy(d_over, &stop, sizeof(bool), cudaMemcpyHostToDevice) ;
		parallel<<<blocks,threads>>>(d_vertex,d_verc,d_edge,d_dist,d_over,d_weight);
		cudaError_t err = cudaGetLastError();
		if(err!=cudaSuccess)
		{
			printf("Error: %s\n", cudaGetErrorString(err));
			printf("Not Ok");
		}
		cudaMemcpy(&stop, d_over, sizeof(bool), cudaMemcpyDeviceToHost) ;
		//printf("a");

	}
	
	
	
//	timer2.Stop();
//	printf("(Processing)Time Elapsed :%lfms\n",timer2.Elapsed());
//	timer2.Start();

	cudaMemcpy(h_final,d_dist,numBytes,cudaMemcpyDeviceToHost);
//	timer2.Stop();
//	printf("(Memory Transfer)Time Elapsed :%lfms\n",timer2.Elapsed());
	//cudaMemcpy(&h_verc,d_verc,sizeof(int),cudaMemcpyDeviceToHost);
	
	printf("Distance\n");
	for(i=0;i<h_verc;i++)
		printf("%d ",h_final[i]);
	printf("\n");
	getchar();
	cudaFree(d_vertex);
	cudaFree(d_verc);
	cudaFree(d_source);
	cudaFree(d_dist);
	cudaFree(d_edge);
	free(h_dist);
	free(h_edge);
	free(h_vertex);
}
