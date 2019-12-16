#include "mpi.h"
#include <Windows.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <cstddef>
#include <float.h>
using namespace std;

int ProcNum = 0, ProcRank = 0;
int* pReceiveNum;  // Количество элементов, посылаемых процессом
int* pReceiveInd;  // Индекс элемента данных в результирующем векторе

struct crsMatrix
{
	int N;//размер матрицы N*N
	int NZ;//количество ненулевых элементов
	//массив значений (размеров NZ)	
	double* Value;
	//массив номеров столбцов (размером NZ)
	int* Col;
	//массив индексов строк	(размером N + 1)
	int* RowIndex;
};

void InitializeMatrix(int n, int NZ, crsMatrix* mtx)
{
	mtx[0].N = n;
	mtx[0].NZ = NZ;
	mtx[0].Value = new double[NZ];
	mtx[0].Col = new int[NZ];
	mtx[0].RowIndex = new int[n + 1];
}

void GenerateRegularCRS(int n, int cntInRow, crsMatrix* mtx)
{
	int i, j, k, f, tmp, notNull, c;
	notNull = cntInRow * n;
	for (i = 0; i < n; i++) {
		// Формируем номера столбцов в строке i 
		for (j = 0; j < cntInRow; j++)
		{
			do
			{
				mtx[0].Col[i * cntInRow + j] = rand() % n;
				f = 0;
				for (k = 0; k < j; k++)
					if (mtx[0].Col[i * cntInRow + j] == mtx[0].Col[i * cntInRow + k])
						f = 1;
			} while (f == 1);
		}
		// Сортируем номера столбцов в строке i 
		for (j = 0; j < cntInRow - 1; j++)
			for (k = 0; k < cntInRow - 1; k++)
				if (mtx[0].Col[i * cntInRow + k] > mtx[0].Col[i * cntInRow + k + 1])
				{
					tmp = mtx[0].Col[i * cntInRow + k];
					mtx[0].Col[i * cntInRow + k] = mtx[0].Col[i * cntInRow + k + 1];
					mtx[0].Col[i * cntInRow + k + 1] = tmp;
				}
	}
	// Заполняем массив значений 
	for (i = 0; i < cntInRow * n; i++)
		mtx[0].Value[i] = rand() % 99 + 1; // Чтобы не было нулевых значений
	// Заполняем массив индексов строк 
	c = 0;
	for (i = 0; i <= n; i++)
	{
		mtx[0].RowIndex[i] = c;
		c += cntInRow;
	}
}

void Print(int n, crsMatrix* mtx)
{
	printf("Matrix in format CRS: \n");
	int k, i;
	k = mtx[0].NZ;
	printf("\n Value: ");
	for (i = 0; i < k; i++) {
		printf(" %.1f", mtx[0].Value[i]);
	}
	printf("\n Col: ");
	for (i = 0; i < k; i++) {
		printf(" %d", mtx[0].Col[i]);
	}
	printf("\n RowIndex: ");
	for (i = 0; i < n + 1; i++) {
		printf(" %d", mtx[0].RowIndex[i]);
	}
	printf("\n");
	fflush(stdout);
}

// Транспонируем матрицу алгоритмом Густавсона для дальнейшего перемножения матриц
crsMatrix Transpose(int N, int NZ, int cntInRow, crsMatrix* A)
{
	int i;
	crsMatrix AT; // Транспонированная матрица А
	InitializeMatrix(N, NZ, &AT);
	memset(AT.RowIndex, 0, (N + 1) * sizeof(int));
	for (i = 0; i < A->NZ; i++)
		AT.RowIndex[A->Col[i] + 1]++;

	int S = 0;
	for (i = 1; i <= A->N; i++)
	{
		int tmp = AT.RowIndex[i];
		AT.RowIndex[i] = S;
		S = S + tmp;
	}
	for (i = 0; i < A->N; i++)
	{
		int j1 = A->RowIndex[i]; int j2 = A->RowIndex[i + 1];
		int Col = i; // Столбец в AT - строка в А
		for (int j = j1; j < j2; j++)
		{
			double V = A->Value[j]; // Значение
			int RIndex = A->Col[j]; // Строка в AT
			int IIndex = AT.RowIndex[RIndex + 1];
			AT.Value[IIndex] = V;
			AT.Col[IIndex] = Col;
			AT.RowIndex[RIndex + 1]++;
		}
	}
	//Print(N, &AT);
	return AT;
}

void Multiplicate(crsMatrix A, crsMatrix B, crsMatrix& C)
{
	int N = A.N;
	vector<int> columns; vector<double> values; vector<int> row_index;
	int rowNZ; row_index.push_back(0);
	for (int i = 0; i < N; i++) {
		rowNZ = 0;
		for (int j = 0; j < N; j++) {
			double sum = 0;
			// Считаем скалярное произведение строк А и BT
			for (int k = A.RowIndex[i]; k < A.RowIndex[i + 1]; k++)
				for (int l = B.RowIndex[j]; l < B.RowIndex[j + 1]; l++)
					if (A.Col[k] == B.Col[l])
					{
						if (A.Col[k] == B.Col[l])
							sum += A.Value[k] * B.Value[l];
						break;
					}

			if (fabs(sum) > 0) {
				columns.push_back(j); values.push_back(sum);  rowNZ++;
			}
		}
		row_index.push_back(rowNZ + row_index[i]);
	}
	InitializeMatrix(N, columns.size(), &C);
	for (int j = 0; j < columns.size(); j++) {
		C.Col[j] = columns[j];
		C.Value[j] = values[j];
	}
	for (int i = 0; i <= N; i++) C.RowIndex[i] = row_index[i];
	//printf("\nResult of multiplication: \n");
	//Print(N, &C);

	//return columns.size();
}

// Для какого процесса какую часть матрицы
void DataDistribution(int N, crsMatrix mtx1, crsMatrix* mtx2, int cntInRow)
{
	int* pSendNum;    	// Количество элементов, посылаемых процессу
	int* pSendInd;    	// Индекс первого элемента данных, посылаемого процессу
	int RestRows = N;
	int RowNum = (N / ProcNum);
	pSendNum = new int[ProcNum];
	pSendInd = new int[ProcNum];
	pSendNum[0] = RowNum * cntInRow;//столько элементов 0-му процессу
	pSendInd[0] = 0; // 0-му с 0-ого индекса
	for (int i = 1; i < ProcNum; i++)//раскидываем оставшиеся строки
	{
		RestRows -= RowNum;
		RowNum = RestRows / (ProcNum - i);
		pSendNum[i] = RowNum * cntInRow;
		pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
	}

	// Рассылка массива Val (значений матрицы)
	MPI_Scatterv(mtx1.Value, pSendNum, pSendInd, MPI_DOUBLE, mtx2[0].Value, pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// Рассылка массива Col (столбцов матрицы)
	MPI_Scatterv(mtx1.Col, pSendNum, pSendInd, MPI_INT, mtx2[0].Col, pSendNum[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);

}


void RowIndexProc(int N, crsMatrix* mtx, int cntInRow)
{
	mtx[0].RowIndex[0] = 0;
	for (int i = 1; i < N + 1; i++)
		mtx[0].RowIndex[i] = mtx[0].RowIndex[i - 1] + cntInRow;
}


int main(int argc, char* argv[])
{
	srand(time(NULL));
	 int N = 2000;
	//N = atoi(argv[1]);
	//количество элементов в строке, для равномерной нагрузки на процессы
	int cntInRow = 2;	//atoi(argv[2]);
	int NZ = cntInRow * N;
	int i;
	crsMatrix mtx;
	crsMatrix mtx2;
	crsMatrix mtx_tmp, mtx_tmp2;
	crsMatrix Answer;
	InitializeMatrix(N, N, &Answer);
	double startwtime, endwtime, dur_paral, dur_consequen;
	int RestRows = N;
	//int* pSendNum;    	// Количество элементов, посылаемых процессу
	//int* pSendInd;    	// Индекс первого элемента данных, посылаемого процессу*/
	//pSendNum = new int[ProcNum];	pSendInd = new int[ProcNum];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	pReceiveNum = new int[ProcNum];	
	pReceiveInd = new int[ProcNum];
	pReceiveInd[0] = 0;	
	pReceiveNum[0] = N / ProcNum;

	for (i = 1; i < ProcNum; i++) {
		RestRows -= pReceiveNum[i - 1];
		pReceiveNum[i] = RestRows / (ProcNum - i);
		pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
	}
	InitializeMatrix(N, NZ, &mtx);	
	InitializeMatrix(N, NZ, &mtx2);
	if (ProcRank == 0)
	{
		GenerateRegularCRS(N, cntInRow, &mtx);
		GenerateRegularCRS(N, cntInRow, &mtx2);
		Print(N, &mtx);
		Print(N, &mtx2);
		MPI_Datatype StructType;
	}
	startwtime = MPI_Wtime();
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&cntInRow, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//fflush(stdout);

	InitializeMatrix(pReceiveNum[ProcRank], pReceiveNum[ProcRank] * cntInRow, &mtx_tmp);
	DataDistribution(N, mtx, &mtx_tmp, cntInRow);
	InitializeMatrix(pReceiveNum[ProcRank], pReceiveNum[ProcRank] * cntInRow, &mtx_tmp2);
	DataDistribution(N, mtx2, &mtx_tmp2, cntInRow);

	RowIndexProc(pReceiveNum[ProcRank], &mtx_tmp, cntInRow);
	RowIndexProc(pReceiveNum[ProcRank], &mtx_tmp2, cntInRow);
	crsMatrix result;
	//int NotNull = 
	Multiplicate(mtx_tmp, mtx_tmp2, result);

	MPI_Gatherv(result.Value, result.N, MPI_DOUBLE, Answer.Value, pReceiveNum, pReceiveInd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(result.Col, result.N, MPI_INT, Answer.Col, pReceiveNum, pReceiveInd, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gatherv(result.RowIndex, result.N, MPI_INT, Answer.RowIndex, pReceiveNum, pReceiveInd, MPI_INT, 0, MPI_COMM_WORLD);
	endwtime = MPI_Wtime();
	dur_paral = endwtime - startwtime;

	if (ProcRank == 0)
	{
		cout << "\nResult of multiplication: \n\n";
		int max = 0;
		for (i = 0; i < ProcNum; i++) {
			max += pReceiveNum[i];
		}
		printf("Values:\n ");
		for (i = 0; i < max; i++)
		{
			//printf(" %.1f, %d, %d ", Answer.Value[i], Answer.Col[i], Answer.RowIndex[i]);
			//if ((Answer.Value[i] != NULL) && (Answer.Value[i] > 0) && (Answer.Value[i] < DBL_MAX - 1))
			printf("%.1f  ", Answer.Value[i]);
		}
		printf("\nColumns:\n ");
		for (i = 0; i < max; i++)
		{
			printf("%d  ", Answer.Col[i]);
		}
		printf("\nRowIndex:\n ");
		int S = 0;
		int tmp;
		//Answer.RowIndex[0] = 0;
		for (i = 0; i < max; i++)
		{
			tmp = Answer.RowIndex[i];
			Answer.RowIndex[i] = S;
			S = S + tmp;
			printf("%d  ", Answer.RowIndex[i]);
		}
		printf("\n\nTime (parallel computing) = %f \n", dur_paral);
		crsMatrix ToCompare;
		InitializeMatrix(N, N, &ToCompare);
		startwtime = MPI_Wtime();
		Multiplicate(mtx, mtx2, ToCompare);
		endwtime = MPI_Wtime();
		//Print(N, &ToCompare);
		dur_consequen = endwtime - startwtime;
		printf("\nTime (consequent computing) = %f \n\n", dur_consequen);
		int res;
		if (dur_consequen >= dur_paral)
		{
			printf("Parallel computing is %f faster\n", dur_consequen - dur_paral);
			printf("The difference is %f times\n", dur_consequen / dur_paral);
		}
		else
		{
			printf("Consequent computing is %f faster\n ", dur_paral - dur_consequen);
			printf("The difference is %f times\n", dur_paral / dur_consequen);
		}
	}
	MPI_Finalize();
	return 0;
}