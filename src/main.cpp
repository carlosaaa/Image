#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num()  (0)
   #define omp_get_num_threads() (0)
   #define omp_get_max_threads() (0)  
   #define omp_set_num_threads(X) 
#endif

#include <time.h>
#include "EasyBMP.h"
#include "Matrix.h"

Matrix box_filter(const Matrix& imSrc, int r)
{	
	Matrix imCum = imSrc;
	imCum.cum_for_every_h();

	Matrix imDst = imCum;
	imDst.cum_diff_r_for_every_h(imCum, r);

	imCum = imDst;
	imCum.cum_for_every_w();

	imDst = imCum;
	imDst.cum_diff_r_for_every_w(imCum, r);

	return imDst;
}

//.导向滤波彩色图模糊
Matrix guided_filter(const Matrix& I, const Matrix& p, int r, double eps)
{
	Matrix m1(p);
	m1.set_all_value_1();
	m1.show("m1");

	Matrix N = box_filter(m1, r);	
	N.show("N");

	//mean_I = boxfilter(I, r) ./ N;
	Matrix mean_I = box_filter(I, r) / N;  
	mean_I.show("mean_I");

	// mean_p = boxfilter(p, r) ./ N;
	Matrix mean_p = box_filter(I, r) / N; 
	mean_p.show("mean_p");

	//mean_Ip = boxfilter(I.*p, r) ./ N;
	Matrix mean_Ip = box_filter(I*p, r) / N; 
	mean_Ip.show("mean_Ip");

	//cov_Ip = mean_Ip - mean_I .* mean_p; % this is the covariance of (I, p) in each local patch.
	Matrix cov_Ip = mean_Ip - mean_I * mean_p;
	cov_Ip.show("cov_Ip");
	
	//mean_II = boxfilter(I.*I, r) ./ N;
	Matrix mean_II = box_filter(I*I, r) / N;
	mean_II.show("mean_II");
	
	//var_I = mean_II - mean_I .* mean_I;
	Matrix var_I = mean_II - mean_I*mean_I;

	//a = cov_Ip ./ (var_I + eps); % Eqn. (5) in the paper;
	Matrix a = cov_Ip / (var_I + eps);
	a.show("a");
	
	//b = mean_p - a .* mean_I; % Eqn. (6) in the paper;
	Matrix b = mean_p - a * mean_I;
	b.show("b");

	//	mean_a = boxfilter(a, r) ./ N;
	Matrix mean_a = box_filter(a, r) / N;
	mean_a.show("mean_a");

	//mean_b = boxfilter(b, r) ./ N;
	Matrix mean_b = box_filter(b, r) / N;
	mean_b.show("mean_b");
	
	//q = mean_a .* I + mean_b; % Eqn. (8) in the paper;
	return mean_a * I + mean_b;

}


int main(int argc, char* argv[])
{
	if (argc == 1)
	{
		printf("a.exe x.bmp\n");
		return 1;
	}

	BMP bmp_input;
	if (!bmp_input.ReadFromFile(argv[1]))
	{
		printf("read %s failed\n", argv[1]);
		return 1;
	}

	printf("%d*%d bitdepth=%d numberofcoler=%d\n", 
		bmp_input.TellWidth(), bmp_input.TellHeight(), 
		bmp_input.TellBitDepth(), bmp_input.TellNumberOfColors());
	
	Matrix input(bmp_input);
	input.show("input");

	Matrix guided = input;
	guided.show("guided");

	int r = 4;
	double eps = 0.01;
	
	int begin = clock();
	Matrix output = guided_filter(guided, input, r, eps);
	int end = clock();

	printf("%d %d %d\n", begin , end, end-begin);
	output.show("output");

	BMP bmp_out;
	output.toBMP(bmp_out);

	bmp_out.WriteToFile("out.bmp");
		
	return 0;
}


