#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <vector>

using namespace std;

template <class T>
class Matrix {
    	vector<vector<T>> matrix;
    	int rows;
    	int cols;
public:
	Matrix();
	Matrix(int);
	Matrix(int, int);
	Matrix(const Matrix&);
	template<class F>friend ostream& operator<<(ostream&, const Matrix<F>&);
	template<class F>friend Matrix<F> operator+(const Matrix<F>&, const Matrix<F>&);
	template<class F>friend Matrix<F> operator-(const Matrix<F>&, const Matrix<F>&);
	template<class F>friend vector<F> operator*(const Matrix<F>&, const vector<F>&);
	template<class F>friend Matrix<F> operator*(const Matrix<F>&, const F&);
	Matrix<T> operator=(const Matrix<T>&);
	Matrix<T> operator+=(const Matrix<T>&);
	Matrix<T> operator-=(const Matrix<T>&);
	T& operator()(int, int);
	const T& operator()(int, int) const;
	T scalar(const Matrix<T>&, const Matrix<T>&); 
	vector<T>& operator[](int);
	const vector<T>& operator[](int) const;
	int row() const;
	int col() const;

};

//конструкторы

template<class T>Matrix<T>::Matrix(): rows(0), cols(0) {} 
template<class T>Matrix<T>::Matrix(int m, int n): rows(m), cols(n) {
	vector<T> v(cols); 
	for(int i = 0; i < rows; i++) matrix.push_back(v);
}  
template<class T>Matrix<T>::Matrix(int n): rows(n), cols(n) {
	vector<T> v(cols); 
	for(int i = 0; i < rows; i++) matrix.push_back(v); 
}
template<class T>Matrix<T>::Matrix(const Matrix& M): rows(M.rows), cols(M.cols) {
	vector<T> v(cols);
	for(int i = 0; i < rows; i++) matrix.push_back(v);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) matrix[i][j] = M(i, j);
	}
}

//операторы перегрузки

template<class F> ostream& operator<<(ostream& os, const Matrix<F>& M) {
	for (int i = 0; i < M.rows; i++) {
		for (int j = 0; j < M.cols; j++) os << setprecision(2) << setw(12) << M(i, j);
		os << endl;
	}
    	return os;
}
template<class F> Matrix<F> operator+(const Matrix<F>& A, const Matrix<F>& B) {
	try {
		if (A.rows != B.rows || A.cols != B.cols) throw "Error!";
	}
	catch (const char* str) {cout << str << endl;}
	Matrix<F> Res(A.rows, A.cols);
	for (int i = 0; i < Res.rows; i++) {
		for (int j = 0; j < Res.cols; j++) Res(i, j) = A(i, j) + B(i, j);
	}
    	return Res;
}
template<class F> Matrix<F> operator-(const Matrix<F>& A, const Matrix<F>& B) {
    	try {
		if (A.rows != B.rows || A.cols != B.cols) throw "Error!";
	}
	catch (const char* str) {cout << str << endl;}
	Matrix<F> Res(A.rows, A.cols);
	for (int i = 0; i < Res.rows; i++) {
		for (int j = 0; j < Res.cols; j++) Res(i, j) = A(i, j) - B(i, j);
	}
    	return Res;
}
template<class F> vector<F> operator*(const Matrix<F>& A, const vector<F>& v) {
	try {
		if (A.cols != v.size()) throw "Error!";
	}
	catch (const char* str) {cout << str << endl;}
	vector<F> res(A.rows);
	for (int i = 0; i < A.rows; i++) {
			res[i] = 0;
			for (int j = 0; j < A.cols; j++) res[i] += A(i, j) * v[j];	
	}
	return res;
}	
template<class F> Matrix<F>operator*(const Matrix<F>& A, const F& a) {
	Matrix<F> Res(A.rows, A.cols);
	for (int i = 0; i < Res.rows; i++) {
		for (int j = 0; j < Res.cols; j++) Res(i, j) = A(i, j) * a;
	}
	return Res;
}
template<class T> Matrix<T> Matrix<T>::operator=(const Matrix<T>& B) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) matrix[i][j] = B(i, j);
	}
	return *this;
}
template<class T>Matrix<T>Matrix<T>::operator+=(const Matrix& B) {
	try {
		if (rows != B.rows || cols != B.cols) throw "Error!";
	}
	catch (const char* str) {cout << str << endl;}
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) matrix[i][j] += B(i, j);
	}
    	return *this;
}
template<class T>Matrix<T>Matrix<T>::operator-=(const Matrix& B) {
	try {
		if (rows != B.rows || cols != B.cols) throw "Error!";
	}
	catch (const char* str) {cout << str << endl;}
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) matrix[i][j] -= B(i, j);
	}
    	return *this;
}
template<class T> T& Matrix<T>::operator()(int i, int j) {return matrix[i][j];}
template<class T> const T& Matrix<T>::operator()(int i, int j) const {return matrix[i][j];}
template<class T> vector<T>& Matrix<T>::operator[](int i) {return matrix[i];}
template<class T> const vector<T>& Matrix<T>::operator[](int i) const {return matrix[i];}

//функции

template<class T> int Matrix<T>::row() const {return rows;}
template<class T> int Matrix<T>::col() const {return cols;}
template<class T> T Matrix<T>::scalar(const Matrix<T>& A, const Matrix<T>& B) {
	try {
		if (A.rows != B.rows || A.cols != B.cols) throw "Error!";
	}
	catch (const char* str) {cout << str << endl;}
	T res = 0;
	for (int i = 0; i < A.rows; i++) {
		for (int j =0; j < A.cols; j++) res += A(i, j) * B(i, j);
	}
	return res;
}
//для векторов
template<class T> vector<T> operator+(const vector<T>& a, const vector<T>& b) {
	int n = a.size();
	vector<T> res(n);
	for (int i = 0; i < n; i++) res[i] = a[i] + b[i];
	return res;
}
template<class T> ostream& operator<<(ostream& os, const vector<T>& v) {
	int n = v.size();
	for (int i = 0; i < n; i++) os << setprecision(2) << setw(10) << v[i];
    	return os;
}
vector<double> linspace(double start, double end, int num) {
	vector<double> linspaced;
  	if (!num) {return linspaced;}
  	if (num == 1) {
      		linspaced.push_back(start);
      		return linspaced;
    	}
  	double delta = (end - start) / (num - 1);
  	for(int i = 0; i < num - 1; i++) linspaced.push_back(start + delta * i);
  	linspaced.push_back(end);                     
  	return linspaced;
}
