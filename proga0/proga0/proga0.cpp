#include <iostream>
#include <string>
#include <fstream>
using namespace std;

ifstream f1("input.txt", ios::in);
ofstream f2("output.txt", ios::out);

class Matrix
{
	int** M; //матрица
	int N;   //размерность матрицы
public:
	Matrix() //default-конструктор
	{
		N = 0; //нужно что-то ещё?
		M = new int*[0];
	}
	Matrix(int value) 
	{
		N = 1;
		M = (int**) new int*;
		M[0] = (int*) new int;
		M[0][0] = value;
	}
	Matrix(unsigned int k, int* value) // инициализация диагональной матрицы
	{
		N = k;
		M = new int* [N];
		for (int i = 0; i < N; i++)
			M[i] = new int[N];

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				if (i == j)
					M[i][j] = value[i];
				else
					M[i][j] = 0;
			}
	}
	Matrix(unsigned int k, bool value) //инициализация квадратной матрицы
	{
		N = k;
		M = (int**) new int* [N];
		for (int i = 0; i < N; i++)
			M[i] = (int*)new int[N];

		for (int i = 0; i < N; i++) 
			for (int j = 0; j < N; j++)
				f1>>M[i][j];
		
	}
	Matrix(unsigned int V, int value) // инициализировали К диагональными элементами
	{
		this->N = V;
		M = new int* [N];
		for (int i = 0; i < N; i++)
			M[i] = new int[N];

		for (int i = 0; i < N; i++) 
			for (int j = 0; j < N; j++)
				M[i][j] =0;
		for (int i = 0; i < N; i++)
			M[i][i] = value;
	}
	
	/*~Matrix()
	{
		for (int i = 0; i < this->N; i++)
			delete[] M[i];
		delete[]M;
	}*/
	int get_element(int i, int j) const
	{
		 return M[i][j];
	}
	void change_element(int i, int j, int value) 
	{
		M[i][j] = value;
	}
	Matrix operator + (const Matrix& second) const
	{
		if (this->N != second.N)
		{
			f2 << "ERROR" << "\n";
			return 0;
		}
		else {
			Matrix Sum = Matrix(second.N, 0);
			for (int i = 0; i < this->N; i++)
				for (int j = 0; j < this->N; j++)
					Sum.change_element(i, j, get_element(i, j) + second.get_element(i, j));
			return Sum;
		}
	}
	Matrix operator - (const Matrix& second)const
	{
		if (this->N != second.N)
		{
			f2 << "ERROR" << "\n";
			return 0;
		}
		else {
			Matrix Sum = Matrix(second.N, 0);
			for (int i = 0; i < this->N; i++)
				for (int j = 0; j < this->N; j++)
					Sum.change_element(i, j, get_element(i, j) - second.get_element(i, j));
			return Sum;
		}
	}
	Matrix operator * (const Matrix& second)const
	{
		if (this->N != second.N)
		{
			f2 << "ERROR" << "\n";
			return 0;
		}
		else {
			int P;
			Matrix Sum = Matrix(second.N, 0);
			for (int i = 0; i < this->N; i++)
				for (int j = 0; j < this->N; j++) {
					P = 0;
					for (int k = 0; k < this->N; k++)
						P += get_element(i, k) * second.get_element(k, j);
					Sum.change_element(i, j, P);
				}
			return Sum;
		}
	}
	bool operator == (const Matrix& second)const
	{
		if (this->N != second.N)
		{
			f2 << "ERROR" << "\n";
			return 0;
		}
		else {
			int P = 0;
			for (int i = 0; i < this->N; i++)
				for (int j = 0; j < this->N; j++) {
					if (get_element(i, j) != second.get_element(i, j))
					{
						P = 1;
						break;
					}

				}
			if (P == 1)
				return false;
			else
				return true;
		}
		
	}
	bool operator != (const Matrix& second)const
	{
		if (this->N != second.N)
		{
			f2 << "ERROR" << "\n";
			return 0;
		}
		else {
			int P = 0;
			for (int i = 0; i < this->N; i++)
				for (int j = 0; j < this->N; j++) {
					if (get_element(i, j) != second.get_element(i, j))
					{
						P = 1;
						break;
					}

				}
			if (P == 1)
				return true;
			else
				return false;
		}

	}
	Matrix operator !( ) const //транспонирование
	{
		Matrix Sum = Matrix(this->N, 0);
		for (int i = 0; i < this->N; i++)
			for (int j = 0; j < this->N; j++)
				Sum.change_element(i, j, this->get_element(j, i));
		return Sum;
	}
	void print() const
	{
		for (int i = 0; i < this->N; i++) {
			for (int j = 0; j < this->N; j++)
				f2 << this->get_element(i, j) << " ";
			f2 << "\n";
		}
	}
	Matrix operator()(int p, int m) const //Минор
	{
		Matrix Sum = Matrix(this->N - 1, 0);
		for (int i = 0; i < p-1; i++)
			for (int j = 0; j < m-1; j++)
				Sum.change_element(i, j, this->get_element(i, j));

		for (int i = p-1; i < this->N-1; i++)
			for (int j = 0; j < m-1; j++)
				Sum.change_element(i, j, this->get_element(i+1, j));

		for (int i = 0; i < p - 1; i++)
			for (int j = m - 1; j < this->N - 1; j++)
				Sum.change_element(i, j, this->get_element(i, j+1));
		
		for (int i = p - 1; i < this->N - 1; i++)
			for (int j = m - 1; j < this->N - 1; j++)
				Sum.change_element(i, j, this->get_element(i + 1, j + 1));
		return Sum;
	}
	int* operator[](int p) const //обращение к строке
	{
		return M[p - 1];
	}
	int* operator ()(int p) const //обращение к столбцу
	{
		int* E = new int[this->N];
		for (int i = 0; i < this->N; i++)
			E[i] = M[i][p-1];
		return E;
	}
	//как создать оператор к типу int* 
};
int main()
{
	int N;
	int k;
	f1 >> N;
	f1 >> k;
	Matrix K(N,k);
	bool y = true;
	Matrix A(N, y);
	Matrix B(N, y);
	Matrix C(N, y);
	Matrix D(N, y);
	Matrix U = (B * (!C) + A + K) * (!D);
	U.print();
	//f2 << "\n";
	//Matrix Z = U(1, 2);
	//Z.print();
	/*K.~Matrix();
	A.~Matrix();
	B.~Matrix();
	C.~Matrix();
	D.~Matrix();
	U.~Matrix();*/
	return 0;
}