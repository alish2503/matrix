#include "Matrix2.h"
#include "Poisson-equation/func.h"
#include <fstream>
#include <unistd.h>
#include <tuple>
#include <chrono>

using namespace std::chrono;

typedef Matrix<double> dm;
typedef vector<double> dv;

dm inv(dm);
dm matr_sweep(const dm&, const dm&);
dm Ax(const dm&); 
tuple<int, double, int, dm> CG(const dm&, double, int);

int main(void) {
	int n = 64;
	int nm1 = n - 1;
	int nm2 = n - 2;
	int np = n + 1;
	double h = 1.0 / n;
	double h2 = h * h;
	dv x = linspace(0.0, 1.0, np);
	dv y = x;
	dm f(nm1), u(nm1), ua(nm1), c(nm1);
	for (int j = 0; j < nm1; j++) {
		for (int i = 0; i < nm1; i++) {
			f(j, i) += fun1(x[i + 1], y[j + 1]) * h2;
			if (!j) f(0, i) += south(x[i + 1]);
			else if (j == nm2) f(nm2, i) += north(x[i + 1]);
		}
		f(j, 0) += west(y[j + 1]);
		f(j, nm2) += east(y[j + 1]);
	}
	int opt, flag, iter, max_iter = 100;
	double eps = 1.e-5;
	cout << "Номер метода ";
	cin >> opt;
	auto start = high_resolution_clock::now();
	switch (opt) {
		case 1:
			for (int i = 0; i < nm1; i++) {
				c(i, i) = 4.0;
				if (i > 0) c(i - 1, i) = -1.0;
				if (i < nm1 - 1) c(i + 1, i) = -1.0;
			}	
			u = matr_sweep(c, f);
			break;
		case 2:
			tie(flag, eps, iter, u) = CG(f, eps, max_iter);
			break;
		default:
			cout << "Ошибка!" << endl;
			return 1;
	}
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start) * 1.0e-6;
	cout << duration.count() << " секунд" << endl;
	double dlt = 0.0;
	double I, J;
	for (int j = 0; j < nm1; j++) {
    		for (int i = 0; i < nm1; i++) { 
    			ua(j, i) = u_a1(x[i + 1], y[j + 1]);     
        		double t = fabs(u(j, i) - ua(j, i));
        		if (t > dlt) {
            			I = i; J = j; dlt = t;
            		}
            	}
        }
        if (opt == 2) cout << "flag = " << flag << ", eps = " << eps << ", число итераций = " << iter << endl << endl;
        cout << "dlt = " << dlt << ", x[i]= " << x[I] << ", y[j]= " << y[J] << endl << endl;
        ofstream ofs("data.txt");
	for (int j = 0; j < nm1; j++) {
		for (int i = 0; i < nm1; i++) ofs << x[i + 1] << ' ' << y[j + 1] << ' ' << u(j, i) << endl;
	}
	if (execlp("Poisson-equation/plot","Poisson-equation/plot","./plot", "./plot", NULL) < 0) {
	 	perror( "failed");
	 	return 1;
	}
	ofs.close();
	return 0;
}

dm inv(dm a) {
	int n = a.row();
	dm e(n);
	for (int i = 0; i < n; i++) e(i, i) = 1;
	double eps = 1.0e-100;
	for (int k = 0; k < n - 1; k++) {
		double max = fabs(a(k, k));
		int pos = k;
		for (int t = k + 1; t < n; t++) {
			if (fabs(a(t, k)) > max) {
				max = fabs(a(t, k));
				pos = t;
			}
		}
		try {
			if (max < eps) throw "Error!";		
		}
		catch (const char* str) {cout << str << endl;}
		if (pos != k) {
			for (int j = k; j < n; j++) swap(a(pos, j), a(k, j));
			for (int s = 0; s < n; s++) swap(e(pos, s), e(k, s));	
			
		}
		for (int i = k + 1; i < n; i++) {
			double tmp = a(i, k) / a(k, k);
			a(i, k) = 0;
			for (int j = k + 1; j < n; j++) a(i, j) -= tmp * a(k, j);
			for (int s = 0; s < n; s++) e(i, s) -= tmp * e(k, s);
		}
	}
	dm x(n);
	for (int k = 0; k < n; k++) {
		for (int i = n - 1; i >= 0; i--) {
			x(i, k) = e(i, k);
			double sum = 0 ;
			for (int j = n - 1; j > i; j--) sum += a(i, j) * x(j, k);
			x(i, k) = (x(i, k) - sum) / a(i, i);
		}
	}
	return x;
}
dm matr_sweep(const dm& c, const dm& f) {
	int n = c.row();
	vector<dm> alph;
	dm x(n), bet(n);
	alph.push_back(inv(c));
	bet[0] = alph[0] * f[0];
	for (int i = 0; i < n - 1; i++) {
		alph.push_back(inv(c - alph[i]));
		bet[i + 1] = alph[i + 1] * (bet[i] + f[i + 1]);
	}
	x[n - 1] = bet[n - 1];
	for (int i = n - 2; i >= 0; i--) x[i] = alph[i] * x[i + 1] + bet[i];
	return x;
}
dm Ax(const dm& x) {
	int n = x.row();
	dm matr(n);
	for (int i = 0; i < n; i++) {
		matr(0, i) = 4 * x(0, i) - x(1, i);
		if (i > 0) matr(0, i) -= x(0, i - 1);
		if (i < n - 1) matr(0, i) -= x(0, i + 1);
	}
	for (int j = 1; j < n - 1; j++) {
		for (int i = 0; i < n; i++) {
			matr(j, i) = 4 * x(j, i) - x(j - 1, i) - x(j + 1, i);
			if (i > 0) matr(j, i) -= x(j, i - 1);
			if (i < n - 1) matr(j, i) -= x(j, i + 1);
		}
	}
	for (int i = 0; i < n; i++) {
		matr(n - 1, i) = 4 * x(n - 1, i) - x(n - 2, i);
		if (i > 0) matr(n - 1, i) -= x(n - 1, i - 1);
		if (i < n - 1) matr(n - 1, i) -= x(n - 1, i + 1);
	}
	return matr;
}	
tuple<int, double, int, dm> CG(const dm& f, double eps, int max_iter) { 
	int n = f.row();
	dm z(n), x(n), p(n), r(n), matr;
	p = r = f;
	double sr1 = matr.scalar(r, r);
	double nf = sqrt(sr1);
	if (!nf) nf = 1.0;
	double res = 1.0; //  res = |r| / nf
	if (res <= eps) {
		eps = res;
		max_iter = 0;
		return make_tuple(0, eps, max_iter, x);
	}
	for (int k = 1; k < max_iter; k++) {
		z = Ax(p);
		double v = sr1 / matr.scalar(p, z);
		x += p * v;
		r -= z * v;
		double sr2 = matr.scalar(r, r);
		res = sqrt(sr2) / nf;
		if (res <= eps) {
			eps = res;
			max_iter = k;
			return 
			make_tuple(0, eps, max_iter, x);
		}
		double m = sr2 / sr1;		
		p = r + p * m;
		sr1 = sr2;
	}
	eps = res;
	return make_tuple(1, eps, max_iter, x);
}

