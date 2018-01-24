数据生成程序
#include nilibddc.h
#include stdlib.h
#include stdio.h
#includestring.h
#includemath.h

#define PI 3.1415926535897932384626433832795028841971
#define ddcChk(f) if (ddcError = (f), ddcError  0) goto Error; else
#ifdef nullChk
#undef nullChk
#endif
#define nullChk(p) if (!(p)) { ddcError = DDC_OutOfMemory; goto Error; } else


double tab_sin;
long FFT_N;

static int	GenerateData(char selected, FILE out, FILE outData, char FILE_PATH, long SelectGroup,char type);
int create_sin_tab(double tab_sin);
double sin_tab(double pi);
double cos_tab(double pi);

typedef struct Dcomplex
{
	double real;        
	double imag;        
}Dcomplex;

#define PI 3.1415926535897932384626433832795028841971  

void window_function(Dcomplex data,char type);
void c_plus(Dcomplex a, Dcomplex b, Dcomplex c);
void c_mul(Dcomplex a, Dcomplex b, Dcomplex c);
void c_sub(Dcomplex a, Dcomplex b, Dcomplex c); 
void fft(int N, Dcomplex f[]);

#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr  

void c_plus(Dcomplex a, Dcomplex b, Dcomplex c)
{
	c-real = a.real + b.real;
	c-imag = a.imag + b.imag;
}

void c_sub(Dcomplex a, Dcomplex b, Dcomplex c)
{
	c-real = a.real - b.real;
	c-imag = a.imag - b.imag;
}

void c_mul(Dcomplex a, Dcomplex b, Dcomplex c)
{
	c-real = a.real  b.real - a.imag  b.imag;
	c-imag = a.real  b.imag + a.imag  b.real;
}

void Wn_i(int n, int i, Dcomplex Wn, char flag)
{
	Wn-real = cos_tab(2  PIi  n);
	if (flag == 1)
		Wn-imag = -sin_tab(2  PIi  n);
	else if (flag == 0)
		Wn-imag = -sin_tab(2  PIi  n);
}


void fft(int N, Dcomplex f[])
{
	Dcomplex t, wn; 
	int i, j, k, m, n, l, r, M;
	int la, lb, lc;
	----计算分解的级数M=log2(N)----
	for (i = N, M = 1; (i = i  2) != 1; M++);
	----按照倒位序重新排列原信号----
	for (i = 1, j = N  2; i = N - 2; i++)
	{
		if (ij)
		{
			t = f[j];
			f[j] = f[i];
			f[i] = t;
		}
		k = N  2;
		while (k = j)
		{
			j = j - k;
			k = k  2;
		}
		j = j + k;
	}

	----FFT算法----
	for (m = 1; m = M; m++)
	{
		la = (int)pow(2, m); la=2^m代表第m级每个分组所含节点数       
		lb = la  2;    lb代表第m级每个分组所含碟形单元数  
						同时它也表示每个碟形单元上下节点之间的距离  
						----碟形运算----
		for (l = 1; l = lb; l++)
		{
			r = (l - 1)(int)pow(2, M - m);
			for (n = l - 1; nN - 1; n = n + la) 遍历每个分组，分组总数为Nla  
			{
				lc = n + lb;  n,lc分别代表一个碟形单元的上、下节点编号       
				Wn_i(N, r, &wn, 1);wn=Wnr  
				c_mul(f[lc], wn, &t);t = f[lc]  wn复数运算  
				c_sub(f[n], t, &(f[lc]));f[lc] = f[n] - f[lc]  Wnr  
				c_plus(f[n], t, &(f[n]));f[n] = f[n] + f[lc]  Wnr  
			}
		}
	}
}

int main(int argc, char argv[])
{
	int	ddcError = 0;
	FILE in,out,outData;
	char SelectName[100] = { 0 },FILE_PATH[100],id[100],idInfo[100],postfix[20],type[20];
	long SelectGroup;
	fopen_s(&in, input.txt, r);
	fscanf_s(in, %s%s%ld%s%s, FILE_PATH,100,SelectName, 100,&SelectGroup,id,100,&type,20);
	strcpy_s(idInfo, 100, id);

	strcpy_s(postfix, 20, .csv);
	strcat_s(id, 100, postfix);
	fopen_s(&outData, id, w);
	if (strcmp(SelectName,index)==0)
		fopen_s(&out, index.txt, w);
	else if(strcmp(SelectName,time)==0)
		fopen_s(&out, time.txt, w);
	else {
		strcpy_s(postfix, 20, .json);
		strcat_s(idInfo, 100, postfix);
		fopen_s(&out, idInfo, w);
	}
	printf(%sn, FILE_PATH);
	ddcChk(GenerateData(SelectName, out, outData, FILE_PATH, SelectGroup,type));
Error
	if (ddcError  0)
		printf(nError %sn, DDC_GetLibraryErrorDescription(ddcError));
	else
		printf(nNo errors.n);
	fclose(in);
	fclose(out);
	fclose(outData);
	return 0;
}

static int GenerateData(char SelectName,FILE out,FILE outData,char FILE_PATH,long SelectGroup,char type)
{
	int				ddcError = 0, length;
	DDCFileHandle	file = 0;
	DDCChannelGroupHandle	groups = 0;
	DDCChannelHandle		channels = 0;
	unsigned int			i, j,k, numChannels, numGroups;
	unsigned __int64		numDataValues,count=0;
	char					name = 0, time = 0,data1[1000][30];
	double					data = 0, avg = 0, dataMax, dataMin, FFTMax, FFTMin, data0[100000];;
	float                   dataP;
	Dcomplex		    	xin;
	long u;

	memset(data0,0,100000sizeof(float));
	ddcChk(DDC_OpenFileEx(FILE_PATH, TDMS, 1, &file));

	ddcChk(DDC_GetNumChannelGroups(file, &numGroups));
	nullChk(groups = calloc(numGroups, sizeof(DDCChannelGroupHandle)));
	ddcChk(DDC_GetChannelGroups(file, groups, numGroups));
	for (i = 0; i  numGroups; ++i)
	{
		if (SelectGroup != -1&&SelectGroup!=i)
			continue;
		ddcChk(DDC_GetChannelGroupStringPropertyLength(groups[i], DDC_CHANNELGROUP_NAME, &length));
		nullChk(time = malloc(length + 1));
		ddcChk(DDC_GetChannelGroupProperty(groups[i], DDC_CHANNELGROUP_NAME, time, length + 1));

		if (strcmp(SelectName, time) == 0)
		{
			fprintf(out,%sn,time);
			continue;
		}

		ddcChk(DDC_GetNumChannels(groups[i], &numChannels));
		nullChk(channels = calloc(numChannels, sizeof(DDCChannelHandle)));
		ddcChk(DDC_GetChannels(groups[i], channels, numChannels));


		for (j = 0; j  numChannels; ++j)
		{
			ddcChk(DDC_GetChannelStringPropertyLength(channels[j], DDC_CHANNEL_NAME, &length));
			nullChk(name = malloc(length + 1));
			ddcChk(DDC_GetChannelProperty(channels[j], DDC_CHANNEL_NAME, name, length + 1));

			if (strcmp(SelectName, index) == 0)
			{
				fprintf(out,%sn,name);
				continue;
			}

			if (strcmp(name, SelectName) != 0)
				continue;

			ddcChk(DDC_GetNumDataValues(channels[j], &numDataValues));
			nullChk(dataP = malloc(sizeof(float)  (unsigned int)numDataValues));
			data = malloc(sizeof(double)  (unsigned int)numDataValues);
			ddcChk(DDC_GetDataValues(channels[j], 0, (unsigned int)numDataValues, dataP));

			for (k = 0; k  numDataValues; k++)
			{
				data[k] = dataP[k];
				data0[count] = data[k];
				if(SelectGroup==-1)
					strcpy_s(data1[count],30, time);
				count++;
			}
			free(data);
			data = 0;
		}
		if (strcmp(SelectName,index) == 0)
			goto Error;
	}
	if (strcmp(SelectName, time) == 0)
		goto Error;
	dataMax = data0[0];
	dataMin = data0[0];
	for (k = 0; k  count; k++)
	{
		dataMax = dataMax  data0[k]  dataMax  data0[k];
		dataMin = dataMin  data0[k]  dataMin  data0[k];
		avg += data0[k];
	}
	avg = count;
	for (k = 0; k  numDataValues; k++)
		data0[k] = data0[k] - avg;
	dataMax -= avg;
	dataMin -= avg;

	if (SelectGroup != -1)
	{
		for (FFT_N = 1; FFT_N  numDataValues; FFT_N = 2);

		tab_sin = calloc(FFT_N,sizeof(double));
		create_sin_tab(tab_sin);

		xin = calloc(FFT_N, sizeof(Dcomplex));
		for (u = 0; u  FFT_N; u++)
		{
			xin[u].real = data0[u];
			xin[u].imag = 0;
		}
		window_function(xin,type);
		fft(FFT_N,xin);
		for (u = 0; uFFT_N; u++)                           求变换后结果的模值，存入复数的实部部分  
			xin[u].real=sqrt(xin[u].realxin[u].real+xin[u].imagxin[u].imag)FFT_N2;

		FFTMax = xin[0].real;
		FFTMin = xin[0].real;
		for (u = 0; u FFT_N; u++)
		{
			FFTMax = FFTMax  xin[u].real  FFTMax  xin[u].real;
			FFTMin = FFTMin  xin[u].real  FFTMin  xin[u].real;
		}
		fprintf(outData, data,data_fftn);
		for (u = 0; u  numDataValues; u++)
			if(uFFT_N2)
				fprintf(outData, %lf,%lfn, data0[u], xin[u].real);
			else fprintf(outData, %lf,-1n, data0[u]);
		free(tab_sin);
		free(xin);
	}
	else{
		fprintf(outData, data,timen);
		for (k = 0; k  count; k++)
			fprintf(outData, %lf,%sn, data0[k],data1[k]);
	}
	fprintf(out, {n	max%lf,n	min%lf,n	FFTmax%lf,n	FFTmin%lf,n	name%s,n	rate%I64d,n	time%sn}n,
		dataMax, dataMin,FFTMax,FFTMin,SelectName,numDataValues,time);

Error
	if (groups)
		free(groups);
	if (name)
		free(name);
	if (file)
		DDC_CloseFile(file);
	if (data)
		free(data);
	if (channels)
		free(channels);
	if (time)
		free(time);
	return ddcError;
}
void window_function(Dcomplex data,char type)
{
	long i;
	if (strcmp(type,hanning)==0)
		for (i = 0; i  FFT_N; i++)
			data[i].real = (data[i].real  0.5(1 - cos(2  PI(i+1)  (FFT_N+1))))2;
	else if(strcmp(type,hamming)==0)
		for (i = 0; i  FFT_N; i++)
			data[i].real = (data[i].real  (0.54 - 0.46cos(2  PIi  (FFT_N-1))))2;
}
int create_sin_tab(double tab_sin) 
{ 
	int i;   
	for (i = 0; i = FFT_N  4; i++)  
		tab_sin[i] = sin(2  PIi  FFT_N);
	return 0;
}
double sin_tab(double pi)
{
	int n;
	double result;
	n = (int)(piFFT_N  2  PI);
	if (n = 0 && n = FFT_N  4)
		result = tab_sin[n];
	else if (n  FFT_N  4 && n  FFT_N  2)
	{
		n -= FFT_N  4;
		result = tab_sin[FFT_N  4 - n];
	}
	else if (n = FFT_N  2 && n  3  FFT_N  4)
	{
		n -= FFT_N  2;
		result = -tab_sin[n];
	}
	else if (n = 3  FFT_N  4 && n  3  FFT_N)
	{
		n = FFT_N - n;
		result = -tab_sin[n];
	}
	return result;
}
double cos_tab(double pi) 
{ 
	double pi2;   
	pi2 = pi + PI  2;   
	if (pi22  PI)     
		pi2 -= 2  PI;   
	return sin_tab(pi2);    
}
