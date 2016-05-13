function Kernels=calcmutikernel(kerneltype,kernelpara,x1,x2)

[n1,d1]=size(x1);
if(nargin==4)
[n2,d2]=size(x2);
if(d1~=d2)error('data error');end;
elseif(nargin==3)
n2=n1;d2=d1;    
end;

count_kernel_type=length(kernelpara);
if(length(kerneltype)~=count_kernel_type)error('para error');end;
k_count=0;
for i=1:count_kernel_type,
k_count=k_count+length(kernelpara{i});
end;

Kernels=zeros(n1,n2,k_count);
iter_kernels=1;
for i=1:count_kernel_type,
    for j=1:length(kernelpara{i}),
    single_kernel_type=kerneltype{i};
    single_kernel_para=kernelpara{i}(j);
    
if(nargin==3)
 
    Kernels(:,:,iter_kernels)=calckernel(single_kernel_type,single_kernel_para,x1);
elseif(nargin==4)
 
    Kernels(:,:,iter_kernels)=calckernel(single_kernel_type,single_kernel_para,x1,x2);
end;

 

    iter_kernels=iter_kernels+1;
    end;
     
 
end;
 









