#include <iostream>
#include <conio.h>
#include <math.h>
#include <iomanip>
using namespace std;

double dispersion(int N, int m, double* x, double* y,double* solutions) 
{
	double sum = 0;
	double* save = new double[N];

	for (int i = 0; i < N; i++)
	{
		//считаю разность в скобке
		for (int j = 0; j <= m; j++)
		{
			y[i] -= solutions[j] * pow(x[i], j);
		}

		//возвожу результат разности значения функции в точке с полиномом m-ой степени в квадрат
		save[i] = y[i] * y[i];
	}

	//сумма квадрав розности функции и полинома m степени
	for (int i = 0; i < N; i++)
	{
		sum += save[i];
	}

	double dispersion = sum / (N - m - 1);

	return dispersion;
}

void Gauss(double** WorkMatrix, double solutions[],  int n)
{
	for (int i = 0; i < n; i++)			 //Прямой ход метода Гаусса
	{										 //поиск максимального элемента в первом столбце
		double maxElem = fabs(WorkMatrix[i][i]);		 //определим главный элемент i - ого столбца
		int lineNumb = i;

		for (int j = i; j < n; j++)
		{
			if (fabs(WorkMatrix[j][i]) > maxElem)
			{
				maxElem = fabs(WorkMatrix[j][i]);		//нашли макс элемент в 'i' строке и закунилуи его на верх
				lineNumb = j;
			}
		}

		if (lineNumb != i)					//перестановка строк
		{
			double* buf = WorkMatrix[i];
			WorkMatrix[i] = WorkMatrix[lineNumb];
			WorkMatrix[lineNumb] = buf;
		}

		double leadElem = WorkMatrix[i][i];             //определение ведущего элемента

		for (int z = i; z < n + 1; z++)			//делю строку на первый свобоный элемент
		{
			WorkMatrix[i][z] = WorkMatrix[i][z] / leadElem;
		}

		for (int j = i + 1; j < n; j++)			//зануляю нижнюю строку
		{
			double down_number = WorkMatrix[j][i];
			for (int m = i; m < n + 1; m++)
			{
				WorkMatrix[j][m] = WorkMatrix[j][m] - WorkMatrix[i][m] * down_number;
			}
		}
	}

	for (int i = n - 1; i > 0; i--)			//обратный ход Гаусса
	{
		for (int j = i - 1; j >= 0; j--)
		{
			WorkMatrix[j][n] = WorkMatrix[j][n] - WorkMatrix[j][i] * WorkMatrix[i][n];
		}
	}
	
	for (int i = n - 1, j = n-1; i >= 0; i--,j--)			//записываю корни в массив
	{
		solutions[j] = WorkMatrix[i][n];
		
	}
}

int main()
{
	setlocale(LC_ALL, "ru");
	int N = 6, m = 4;	

	double* x = new double[N];
	double* y = new double[N];

	x[0] = 0; y[0] = 29.5;
	x[1] = 20; y[1] = 18.4;
	x[2] = 40; y[2] = 11.9;
	x[3] = 60; y[3] = 8.6;
	x[4] = 80; y[4] = 5.0;
	x[5] = 100; y[5] = 3.3;

	double* POWERX = new double[2 * m];
	for (int i = 0; i < 2 * m; i++)
	{
		POWERX[i] = 0;
		for (int j = 0; j < N; j++)
			POWERX[i] += pow(x[j], i + 1);
	}

	double** SUMX = new double* [m + 1];
	for (int i = 0; i < m + 1; i++)
		SUMX[i] = new double[m + 1];

	for (int l = 0; l < m + 1; l++)
	{
		for (int j = 0; j < m + 1; j++)
		{
			SUMX[0][0] = N;
			SUMX[l][j] = POWERX[l + j - 1];
		}
	}

	double* PRAW = new double[m + 1];
	for (int l = 0; l < m + 1; l++)
	{
		PRAW[l] = 0;
		for (int i = 0; i < N; i++)
			PRAW[l] += y[i] * pow(x[i], l);
	}

	int n = m + 1;
	double** A = new double* [n];
	for (int i = 0; i < n; i++)
		A[i] = new double[n + 1];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A[i][j] = SUMX[i][j];
		}
		A[i][n] = PRAW[i];
	}

	//cout << endl << "A:" << endl;
	//for (int i = 0; i < n; i++) {
	//	for (int j = 0; j < n + 1; j++) {
	//		cout << setw(15) << A[i][j];
	//	}
	//	cout << endl;
	//}

	double* solutions = new double[m + 1];
	Gauss(A, solutions, m + 1);

	double s = dispersion(N,m,x,y,solutions);

	cout << "P(x) = ";
	for (int i = 0; i < m + 1; i++)
	{
		if (i == 0)
		{
			cout << " " << solutions[i] << " + ";
		}
		else if( i != m)
		{
			cout << " (" << solutions[i] << "*x^" << i << ") + ";
		}
		else
		{
			cout << " (" << solutions[i] << "*x^" << i << ")\n\n";
		}
	}

	cout << "Остаточноя дисперсия: " << s << endl;
	cout << "Среднеквадротическое отклонение: " << pow(s, 0.5) << endl << endl;

	system("pause");
}