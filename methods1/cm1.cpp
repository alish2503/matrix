#include "Matrix1.h"
#include "Poisson-equation/func.h"
#include <fstream>
#include <unistd.h>
#include <chrono>

using namespace std::chrono;

typedef Matrix<double> dm;
typedef vector<double> dv;

void Jacobi(const dm&, dm&, const dm&, double);
void Seidel(dm&, const dm&, double);
void modified_Seidel(dm&, const dm&, double); 
void Sor(dm&, const dm&, double, double);
double c_norm(const dm&);
double matr_max(const dm&);
double matr_min(const dm&);
double approx(const dm&, const dm&, double);

int main(void) {
	int n = 16;
	int np = n + 1;
	int p = 3;
	double h = 1.0 / n;
	double h2 = h * h;
	dv x = linspace(0.0, 1.0, np);
	dv y = x;
	dm f(n), ua(np), u1(np), u2(np);
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < n; j++) f(i, j) = fun1(x[i], y[j]);
	}
	for (int i = 0; i < np; i++) {
    		for (int j = 0; j < np; j++) ua(i, j) = u_a1(x[i], y[j]);   
    	}
	for (int i = 0; i < np; i++) {
		u1(i, 0) = u2(i, 0) = south(x[i]); u1(0, i) = u2(0, i) = west(y[i]);
		u1(i, n) = u2(i, n) = north(x[i]); u1(n, i) = u2(n, i) = east(y[i]);
	}
	double t1 = approx(u1, f, h2);
	double t4 = c_norm(u1 - ua);
	double rs, r, tmp, omega; 
	int opt, m;
	cout << "Номер метода ";
	cin >> opt;
	switch (opt) {
		case 1:
			m = ceil(2 * p * log(10) / (M_PI * M_PI * h2));
			rs = cos(M_PI * h);
			break;
		case 2:
			m = ceil(p * log(10) / (M_PI * M_PI * h2));
			rs = 1.0 - M_PI * M_PI * h2;
			break;
		case 3:
			m = ceil(2 * p * log(10) / (M_PI * h));
			r = cos(M_PI * h);
			omega = 2.0 / (1.0 + sqrt(1.0 - r * r));
			rs = omega - 1.0;
			break;
		case 4:
			m = ceil(p * log(10) / (M_PI * M_PI * h2));
			rs = 1.0 - M_PI * M_PI * h2;
			break;
		default:
			cout << "Ошибка!" << endl;
			return 1;
	}
	cout << endl << "Число итераций " << m << endl;
	cout << "Спектральный радиус " << rs << endl << endl;
	cout << "k \t   |F-AUk|    |F-AUk|/|F-AU0|\t  |Uk-aU|    |Uk-aU|/|U0-aU|\t|Uk-U(k-1)|   Оц-ка пог-сти\t   rs_exp" << endl;  
	auto start = high_resolution_clock::now();
	for (int k = 1; k <= m; k++) {	
		switch (opt) {
			case 1:
				Jacobi(u1, u2, f, h2);
				break;
			case 2:
				Seidel(u2, f, h2);
				break;
			case 3:
				Sor(u2, f, h2, omega);
				break;
			case 4:
				modified_Seidel(u2, f, h2);
				break;
		}
		double t5 = c_norm(u2 - u1);
		double t2 = approx(u2, f, h2);
		double t3 = c_norm(u2 - ua); 
		if (opt == 3) {
			if ((k % 10 == 0) || k == m) {
				printf("%d", k);
				if (k >= 100) printf("%14.7lf", t2); 
				else printf("%15.7lf", t2);
				printf("%20.7lf %12.7lf %12.7lf %22.8lf %14.8lf %13.5lf\n", t2 / t1, t3, t3 / t4, t5, rs * t5 / (1.0 - rs), t5 / tmp); 
			}
		}
		else if (k % 100 == 0 || k == m) {
			printf("%d", k);
			if (k >= 1000) printf(" %14.4lf", t2); 
			else printf(" %15.4lf", t2);
			printf("%18.7lf %12.5lf %15.5lf %18.7lf %14.6lf %13.5lf\n", t2 / t1, t3, t3 / t4, t5, rs * t5 / (1.0 - rs), t5 / tmp); 
		}
		tmp = t5;
		u1 = u2;
	}
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start) * 1.0e-6;
	cout << endl << duration.count() << " секунд" << endl;
	cout << endl << "Umax = " << matr_max(u1) << endl << "Umin = " << matr_min(u1) << endl;
	ofstream ofs("data.txt");
        for (int i = 0; i < np; i++) {
		for (int j = 0; j < np; j++) ofs << x[i] << ' ' << y[j] << ' ' << u1(i, j) << endl;
	}
	if (execlp("./plot", "./plot", NULL) < 0) {
	 	perror( "failed");
	 	return 1;
	}
	ofs.close();
	return 0;
}

void Jacobi(const dm& u1, dm& u2, const dm& f, double h2) {
	int n = u2.row() - 1;
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < n; j++) {
			u2(i, j) = (u1(i - 1, j) + u1(i + 1, j) + u1(i, j - 1) + u1(i, j + 1) + h2 * f(i, j)) / 4;
		}
	}	
}
void Seidel(dm& u, const dm& f, double h2) {
	int n = u.row() - 1;
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < n; j++) {
			u(i, j) = (u(i - 1, j) + u(i + 1, j) + u(i, j - 1) + u(i, j + 1) + h2 * f(i, j)) / 4;
		}
	}	
}
void modified_Seidel(dm& u, const dm& f, double h2)  {
	int n = u.row() - 1;
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < n; j++) {
			if ((i + j) % 2 == 0) {
				u(i, j) = (u(i - 1, j) + u(i + 1, j) + u(i, j - 1) + u(i, j + 1) + h2 * f(i, j)) / 4;
			}
		}
	}
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < n; j++) {
			if ((i + j) % 2 == 1) {
				u(i, j) = (u(i - 1, j) + u(i + 1, j) + u(i, j - 1) + u(i, j + 1) + h2 * f(i, j)) / 4;
			}
		}
	}		
}
void Sor(dm& u, const dm& f, double h2, double omega) {
	int n = u.row() - 1;
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < n; j++) {
			double tmp = (u(i - 1, j) + u(i + 1, j) + u(i, j - 1) + u(i, j + 1) + h2 * f(i, j)) / 4 - u(i, j);
			tmp *= omega;
			u(i, j) += tmp;
		}
	}	
}
double c_norm(const dm& a) {
	int n = a.row();
	double norm = 0.0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double max = fabs(a(i, j));
			if (norm < max) norm = max;
		}
	}
	return norm;
}
double matr_max(const dm& a) {
	int n = a.row();
	double max = a(0, 0);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (max < a(i, j)) max = a(i, j);
		}
	}
	return max;
}
double matr_min(const dm& a) {
	int n = a.row();
	double min = a(0, 0);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (min > a(i, j)) min = a(i, j);
		}
	}
	return min;
}
double approx(const dm& u, const dm& f, double h2) {
	int n = u.row() - 1;
	double norm = 0.0;
	for (int i = 1; i < n; i++) {	
		for (int j = 1; j < n; j++) {
            		double max = fabs((u(i - 1, j) + u(i + 1, j) + u(i, j - 1) + u(i, j + 1) - 4 * u(i, j)) / h2 + f(i, j));
            		if (norm < max) norm = max;
            	}
        }        
	return norm;
}
 

