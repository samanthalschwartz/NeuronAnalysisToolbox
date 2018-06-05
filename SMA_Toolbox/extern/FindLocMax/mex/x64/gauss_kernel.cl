__kernel void gauss_kernel0 (__global float* in,
         int sz0,
         int sz1,
         int sz2,
         float b0,
         float b1,
         float b2,
         float b3,
         float B
        )
{
    const int idy=get_global_id(0); //y index
    const int idz=get_global_id(1); //z index
    float w0,w1,w2,w3;
    float temp;
    int ii=0;
    const int base=idz*sz0*sz1+idy*sz0;
    
    //forward
    w1=w2=w3=in[base];
    
    for (ii=0;ii<sz0;ii++)
    {
        w0=in[base+ii];
        temp=w0*B+(b1*w1+b2*w2+b3*w3)/b0;
        in[base+ii]=temp;
        w3=w2;
        w2=w1;
        w1=temp;
    }
   
    //backward
    w1=w2=w3=in[base+sz0-1];
    for (ii=sz0-1;ii>=0;ii--)
    {
        w0=in[base+ii];
        temp=w0*B+(b1*w1+b2*w2+b3*w3)/b0;
        in[base+ii]=temp;
        w3=w2;
        w2=w1;
        w1=temp;
    }
     
}

__kernel void gauss_kernel1 (__global float* in,
        const int sz0,
        const int sz1,
        const int sz2,
        const float b0,
        const float b1,
        const float b2,
        const float b3,
        const float B
        )
{
    const int idx=get_global_id(0);
    const int idz=get_global_id(1);
    float w0,w1,w2,w3;
    int ii;
    float temp;
    const int base=idz*sz0*sz1+idx;
    //forward
    w1=w2=w3=in[base];
    for (ii=0;ii<sz1;ii++)
    {
        w0=in[base+ii*sz0];
        temp=w0*B+(b1*w1+b2*w2+b3*w3)/b0;
        in[base+ii*sz0]=temp;
        w3=w2;
        w2=w1;
        w1=temp;
    }
    
    //backward
    w1=w2=w3=in[base+sz0*(sz1-1)];
    for (ii=sz1-1;ii>=0;ii--)
    {
        w0=in[base+ii*sz0];
        temp=w0*B+(b1*w1+b2*w2+b3*w3)/b0;
        in[base+ii*sz0]=temp;
        w3=w2;
        w2=w1;
        w1=temp;
    }
}

__kernel void subtract_kernel (__global float* A,
        __global float* B,
        const int sz0,
        const int sz1,
        const int sz2
        )
{
    const int id[3]={get_global_id(0),get_global_id(1),get_global_id(2)};
    const int imsz[3]={sz0,sz1,sz2};
    A[id[2]*imsz[0]*imsz[1]+id[1]*imsz[0]+id[0]]=
            A[id[2]*imsz[0]*imsz[1]+id[1]*imsz[0]+id[0]]
            -B[id[2]*imsz[0]*imsz[1]+id[1]*imsz[0]+id[0]];   
}

__kernel void maxdim0_kernel (__global const float* in,
        __global float* out,
        const int kernelsz,
        const int sz0,
        const int sz1,
        const int sz2,
        const float minval
        )
{
    
    const int id[3]={get_global_id(0),get_global_id(1),get_global_id(2)};
    
    int start =  max(0,id[0]-kernelsz);
    int end =  min(sz0-1,id[0]+kernelsz);
    
    float inpixel=in[id[2]*sz0*sz1+id[1]*sz0+id[0]];
    float maxval=minval;
    
    for (int ii=start;ii<end+1;ii++)
        maxval=max(maxval,in[id[2]*sz0*sz1+id[1]*sz0+ii]);
    
    out[id[2]*sz0*sz1+id[1]*sz0+id[0]]=
            (maxval>inpixel)*-maxval+
            (maxval==inpixel)*maxval;
    
}

__kernel void maxdim1_kernel (__global const float* in,
        __global float* out,
        const int kernelsz,
        const int sz0,
        const int sz1,
        const int sz2,
        const float minval
        )
{
    
    const int id[3]={get_global_id(0),get_global_id(1),get_global_id(2)};
    
    int start =  max(0,id[1]-kernelsz);
    int end =  min(sz1-1,id[1]+kernelsz);
    
    float inpixel=in[id[2]*sz0*sz1+id[1]*sz0+id[0]];
    float maxval=minval;
    
    for (int ii=start;ii<end+1;ii++)
        maxval=max(maxval,fabs(in[id[2]*sz0*sz1+ii*sz0+id[0]]));
    
    out[id[2]*sz0*sz1+id[1]*sz0+id[0]]=(maxval==inpixel);
    
}
