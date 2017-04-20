writing query graph
./graph > input.txt
give vertex_count avgerage_degree variance_degree alphabet_size
3 3 0 1 ->gives a 3 vertex complete graph with only 'A' value present in each node
3 3 0 2 ->gives a 3 vertex complete graph with only 'A' & 'B' value present in each node
10 5 0 1-> gives 10 vertex with 50 edges
writing data graph
./graph >> input.txt
Compilation :  /usr/local/cuda/bin/nvcc -arch=sm_35 -rdc=true  isokernelgpu.cu 
run command
time ./a.out < input.txt

