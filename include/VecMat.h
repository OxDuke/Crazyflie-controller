#ifndef _VECMAT_H_
#define _VECMAT_H_
#include <iostream>
#include <math.h>
#include <assert.h>
using namespace std;
//释放内存
template<class T> inline void V_Dele(T* B)
{
	delete[] B;
	B=NULL;
}

//将向量A的值赋给B
template<class T> inline void V_Copy(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++) B[I_]=A[I_];
}

//将三个值依次赋给B
template<class T> inline void V_Copy(T* B, T x, T y, T z)
{
	B[0]=x;
	B[1]=y;
	B[2]=z;
}

//将六个值依次赋给B
template<class T> inline void V_Copy(T* B, T x, T y, T z, T vx, T vy, T vz)
{
	B[0]=x;
	B[1]=y;
	B[2]=z;
	B[3]=vx;
	B[4]=vy;
	B[5]=vz;
}

//向量B与A每个元素都相等时返回真
template<class T> inline bool V_BoolEqua(const T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++)	if(B[I_]!=A[I_])return false;
	return true;
}

//B[i]=-A[i]
template<class T> inline void V_Opposite(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++) B[I_]=-A[I_];
}


//向量C[i]=A[i]+B[i]
template<class T> inline void V_Add(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]+B[I_];	
}

//向量C[i]=B[i]+A
template<class T> inline void V_Add(T* C, const T* B, T A, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A+B[I_];	
}
//C[i] += B[i]
template<class T> inline void V_Add(T* C, const T* B, int N) {
	for (int I_ = 0; I_<N; I_++) C[I_] += B[I_];
}
//向量C[i]=A[i]-B[i]
template<class T> inline void V_Minus(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]-B[I_];	
}

//向量C[i]=A[i]-B
template<class T> inline void V_Minus(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]-B;
}
//向量C[i]-=A[i]
template<class T> inline void V_Minus(T* C, const T* A, int N) {
	for (int I_ = 0; I_<N; I_++) C[I_] -= A[I_];
}

//向量C[i]=A[i]*B[i]
template<class T> inline void V_Multi(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]*B[I_];	
}

//向量C[i]=A[i]*B
template<class T> inline void V_Multi(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=B*A[I_];
}

//向量A[i]=A[i]*B
template<class T> inline void V_Multi(T* C, T B, int N) {
	for (int I_ = 0; I_<N; I_++) C[I_] *= B;
}

//向量C[i]=A[i]/B[i]
template<class T> inline void V_Divid(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]/B[I_];	
}

//向量C[i]=A[i]/B
template<class T> inline void V_Divid(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]/B;
}

//求内积
template<class T> inline T V_Dot(const T* A, const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++) result+=A[I_]*B[I_];
	return result;
}

//求外积C=AXB,不能用V_Cross(B,B,A)或V_Cross(B,A,B)
template<class T> inline void V_Cross(T* C, const T* A, const T* B)
{
	C[0]=A[1]*B[2]-A[2]*B[1];
	C[1]=A[2]*B[0]-A[0]*B[2];
	C[2]=A[0]*B[1]-A[1]*B[0];
}

//求A[i]=|B[i]|
template<class T> inline void V_Absol(T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++)
	{
		if(B[I_]>=0)
			A[I_]=B[I_];
		else
			A[I_]=-B[I_];
	}
}

//求1-范数
template<class T> inline T V_Norm1(const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++)
	{
		if(B[I_]>=0)
			result+=B[I_];
		else
			result-=B[I_];
	}
	return result;
}

//求2-范数
template<class T> inline T V_Norm2(const T* B, int N)
{
	T result=V_Dot(B,B,N);
	return sqrt(result);
}

//求无穷-范数
template<class T> inline T V_NormInf(const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++) 
	{
		if(B[I_]>=0) {if(B[I_]>result) result=B[I_];}
		else {if(-B[I_]>result) result=-B[I_];}
	}
	return result;
}

//求最大元素
template<class T> inline T V_Max(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;I_++) if(B[I_]>result) result=B[I_];
	return result;
}

//求最大元素
template<class T> inline T V_Max(int & index, const T* B, int N)
{
	T maximal=B[0];
	index=0;
	for(int I_=0;I_<N;I_++) if(B[I_]>maximal) {index=I_;maximal=B[I_];}
	return maximal;
}

//求最小元素
template<class T> inline T V_Min(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;I_++) if(B[I_]<result) result=B[I_];
	return result;
}

//求最小元素
template<class T> inline T V_Min(int& index, const T* B, int N)
{
	T minimal=B[0];
	index=0;
	for(int I_=0;I_<N;I_++) if(B[I_]<minimal) {index=I_;minimal=B[I_];}
	return minimal;
}

//输入流函数
template<class T> inline void V_Input(istream  &input, T* Vec, int N)
{
    for(int I_=0;I_<N;I_++)	input>>Vec[I_];
}

//输出流函数
template<class T> inline void V_Output(ostream &output, const T* Vec, int N)
{
    for(int I_=0;I_<N;I_++) output<<Vec[I_]<<endl;   
}
//输出流函数
template<class T> inline void V_Output(const T* Vec, int N) {
	for (int I_ = 0; I_ < N; I_++) printf("%22.14e\n", Vec[I_]);
}

/**********************************************************************************/
//矩阵，以一维数组表示

//释放内存，与矢量相同

//将九个值依次赋给B
template<class T> inline void M_Copy(T* B, T a11, T a12, T a13, T a21, T a22, T a23, T a31, T a32, T a33)
{
	T temp[9]={a11, a12, a13, a21, a22, a23, a31, a32, a33};
	for(int I_=0;I_<9;I_++)
		B[I_]=temp[I_];
}

template<class T> inline void M_Copy(T* B, T angle, int axis)
{
	assert(axis==1||axis==2||axis==3);
	for(int I_=0;I_<9;I_++) B[I_]=0;
	if(axis==1)
	{
		B[0]=1.0;
		B[4]=B[8]=cos(angle);
		B[5]=sin(angle);
		B[7]=-B[5];
	}
	if(axis==2)
	{
		B[4]=1.0;
		B[0]=B[8]=cos(angle);
		B[6]=sin(angle);
		B[2]=-B[6];
	}
	if(axis==3)
	{
		B[8]=1.0;
		B[0]=B[4]=cos(angle);
		B[1]=sin(angle);
		B[3]=-B[1];
	}
}
/*
//将矩阵B的值逐个赋给数组A
template<class T> inline void M_Copy(T* A, T** B, int N, int M)
{
	int i=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) A[i++]=B[I_][J_];
}


//矩阵B与A每个元素都相等时返回真
template<class T> inline bool M_BoolEqua(T* B, T** A, int N, int M)
{
	int index=0;
	for(int I_=0;I_<N;I_++)	for(int J_=0;J_<M;J_++) if(B[index++]!=A[I_][J_])return false;
	return true;
}
*/
//C[i][j]=A[i][k]*B[k][j],C:NXM,A:NXK,B:KXM,不能M_Multi(C,C,B...
template<class T> inline void M_Multi(T* C, const T* A, const T* B, int N, int M, int K)
{
	int s=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++)
	{
		s=I_*M+J_;
		C[s]=0;
		for(int i=0;i<K;i++) C[s]+=A[I_*K+i]*B[i*M+J_];
	}
}

//B[i][j]=A[j][i],B:NXM,A:MXN
template<class T> inline void M_Tranpose(T* B, const T* A, int N, int M)
{
	int index=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[index++]=A[J_*N+I_];
}

//B[i][j]=B[j][i],B:NXN
template<class T> inline void M_Tranpose(T* B, int N) {
	for (int I_ = 0; I_ < N; I_++) {
		for (int J_ = I_ + 1; J_ < N; J_++) {
			T tmp = B[I_ * N + J_];
			B[I_ * N + J_] = B[J_ * N + I_];
			B[J_ * N + I_] = tmp;
		}
	}
}

//求最大元素
template<class T> inline T M_Max(int& row, int& col, const T* B, int N, int M)
{
	T maximal=B[0];
	row=0;
	col=0;
	int s=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++)
	{
		s=I_*M+J_;
		if(B[s]>maximal) {row=I_; col=J_;maximal=B[s];}
	}
	return maximal;
}

//求最小元素
template<class T> inline T M_Min(int& row, int& col, const T* B, int N, int M)
{
	T minimal=B[0];
	row=0;
	col=0;
	int s=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++)
	{
		s=I_*M+J_;
		if(B[s]<minimal) {row=I_; col=J_;minimal=B[s];}
	}
	return minimal;
}

//输出流函数
template<class T> inline void M_Output(const T* Vec, int N, int M)//ostream &output, 
{
    for(int I_=0;I_<N;I_++)
	{
		//for(int J_=0;J_<M;J_++)	output<<setprecision(8)<<Vec[I_*M+J_]<<",";
		//output<<endl;
		for(int J_=0;J_<M;J_++)	printf("%.6f%s",Vec[I_*M+J_],",");
		printf("\n");
	}
}

//求逆
//	gjelim 
void M_Inverse(double* B, const double* A, int dim, double* wa); 

/***************************************************************/
/*
//二维指针数组表示
//释放内存
template<class T> inline void M_Dele(T** B, int N, int M)
{
	for(int I_=0;I_<N;I_++) {delete[] B[I_];B[I_]=NULL;}
	delete[] B;
	B=NULL;
}

//将A的值赋给B
template<class T> inline void M_Copy(T** B, T** A, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[I_][J_]=A[I_][J_];
}

//将九个值依次赋给B
template<class T> inline void M_Copy(T** B, T a11, T a12, T a13, T a21, T a22, T a23, T a31, T a32, T a33)
{
	B[0][0]=a11;
	B[0][1]=a12;
	B[0][2]=a13;
	B[1][0]=a21;
	B[1][1]=a22;
	B[1][2]=a23;
	B[2][0]=a31;
	B[2][1]=a32;
	B[2][2]=a33;
}

template<class T> inline void M_Copy(T**B, T angle, int axis)
{
	assert(axis==1||axis==2||axis==3);
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[I_][J_]=0;
	if(axis==1)
	{
		B[0][0]=1.0;
		B[1][1]=B[2][2]=cos(angle);
		B[1][2]=sin(angle);
		B[2][1]=-B[1][2];
	}
	if(axis==2)
	{
		B[1][1]=1.0;
		B[0][0]=B[2][2]=cos(angle);
		B[2][0]=sin(angle);
		B[0][2]=-B[2][0];
	}
	if(axis==3)
	{
		B[2][2]=1.0;
		B[0][0]=B[1][1]=cos(angle);
		B[0][1]=sin(angle);
		B[1][0]=-B[0][1];
	}
}
//将数组A的值逐个赋给B
template<class T> inline void M_Copy(T** B, T* A, int N, int M)
{
	int i=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[I_][J_]=A[i++];
}


//矩阵B与A每个元素都相等时返回真
template<class T> inline bool M_BoolEqua(T** B, T** A, int N, int M)
{
	for(int I_=0;I_<N;I_++)	for(int J_=0;J_<M;J_++) if(B[I_][J_]!=A[I_][J_])return false;
	return true;
}

//B[i][j]=-A[i][j]
template<class T> inline void M_Opposite(T** B, T** A, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[I_][J_]=-A[I_][J_];
}

//B[i][j]=A[j][i]
template<class T> inline void M_Tranpose(T** B, T** A, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[I_][J_]=A[J_][I_];
}

//C[i][j]=A[i][j]+B[i][j]
template<class T> inline void M_Add(T** C, T** A, T** B, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) C[I_][J_]=A[I_][J_]+B[I_][J_];	
}

//C[i][j]=B[i][j]+A
template<class T> inline void M_Add(T** C, T** B, T A, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) C[I_][J_]=A+B[I_][J_];	
}

//C[i][j]=A[i][j]-B[i][j]
template<class T> inline void M_Minus(T** C, T** A, T** B, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) C[I_][J_]=A[I_][J_]-B[I_][J_];	
}


//C[i][j]=A[i][k]*B[k][j]
template<class T> inline void M_Multi(T** C, T** A, T** B, int N, int M, int K)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++)
	{
		C[I_][J_]=0;
		for(int i=0;i<K;i++){ C[I_][J_]+=A[I_][i]*B[i][J_];	}
	}
}

//C[i]=A[i][j]*B[j]
template<class T> inline void M_Multi(T* C, T** A, T* B, int N, int M)
{
	for(int I_=0;I_<N;I_++)
	{
		C[I_]=0;
		for(int J_=0;J_<M;J_++)	C[I_]+=A[I_][J_];
	}
}

//C[i][j]=A[i][j]*B
template<class T> inline void M_Multi(T** C, T** A, T B, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) C[I_][J_]=B*A[I_][J_];
}


//C[i][j]=A[i][j]/B[i][j]
template<class T> inline void M_Divid(T** C, T** A, T** B, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) C[I_][J_]=A[I_][J_]/B[I_][J_];	
}

//求A[i][j]=|B[i][j]|
template<class T> inline void M_Absol(T** A, T** B, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++)
	{
		if(B[I_][J_]>=0)
			A[I_][J_]=B[I_][J_];
		else
			A[I_][J_]=-B[I_][J_];
	}
}

//求最大元素
template<class T> inline T M_Max(T** B, int N, int M)
{
	T result=B[0][0];
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) if(B[I_][J_]>result) result=B[I_][J_];
	return result;
}

//求最大元素
template<class T> inline T M_Max(T** B, int N, int M, int& row, int& col)
{
	T maximal=B[0][0];
	row=0;
	col=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) if(B[I_][J_]>result) {row=I_; col=J_;result=B[I_][J_];}
	return index;
}

//求最小元素
template<class T> inline T M_Min(T** B, int N, int M)
{
	T result=B[0][0];
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) if(B[I_][J_]<result) result=B[I_][J_];
	return result;
}

//求最小元素
template<class T> inline T M_Min(T** B, int N, int M, int& row, int& col)
{
	T minimal=B[0][0];
	row=0;
	col=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) if(B[I_][J_]<result) {row=I_; col=J_;result=B[I_][J_];}
	return index;
}

//输入流函数
template<class T> inline void M_Input(std::istream &input, T** Vec, int N, int M)
{
    for(int I_=0;I_<N;I_++)	for(int J_=0;J_<M;J_++) input>>Vec[I_][J_];
}

//输出流函数
template<class T> inline void M_Output(std::ostream &output, T** Vec, int N, int M)
{
    for(int I_=0;I_<N;I_++)
	{
		for(int J_=0;J_<M;J_++)	output<<Vec[I_][J_]<<",";
		output<<std::endl;
	}
}

//求逆
template<class T> inline void swaprows(T** A, int row0, int row1)
{
    T* temp=NULL;
    temp=A[row0];
    A[row0]=A[row1];
    A[row1]=temp;
}

//	gjelim 
void M_Inverse(double** B, double** A, int dim, double** wa); */

#endif