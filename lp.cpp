#include <iostream>
#include <vector>
#include <glpk.h>
using namespace std;

void solve_lp_relaxation(const vector<vector<int>> &a, const vector<int> &b, const vector<int> &c) {
	int m = a.size();
	int n = a[0].size();
	
	int* ia = new int[1 + n * m];
	int* ja = new int[1 + n * m];
	double* ar = new double[1 + n * m];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			ia[i * n + j + 1] = i + 1;
			ja[i * n + j + 1] = j + 1;
			ar[i * n + j + 1] = a[i][j];
		}
	}

	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);

	glp_add_rows(lp, m);
	for (int i = 0; i < m; i++) {
		glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, b[i]);
	}

	glp_add_cols(lp, n);
	for (int j = 0; j < n; j++) {
		glp_set_col_bnds(lp, j + 1, GLP_DB, 0.0, 1.0);
		glp_set_obj_coef(lp, j + 1, c[j]);
	}

	glp_load_matrix(lp, n * m, ia, ja, ar);
	glp_simplex(lp, NULL);

	vector<double> res(n);

	cout << "solution:" << endl;
	cout << "z: " << glp_get_obj_val(lp) << endl;
	cout << "j,\tx,\treduced cost" << endl;
	for (int j = 0; j < n; j++) {
		res[j] = glp_get_col_prim(lp, j + 1);
		cout << j + 1 << "\t" << res[j] << "\t" << glp_get_col_dual(lp, j + 1) << endl;
	}
}


int main(void) {     
    vector<vector<int>> a{
    	vector<int>{7, 19, 30, 22, 30, 44, 11, 21, 35, 14, 29, 18, 3, 36, 42, 87},
    	vector<int>{3, 5, 7, 35, 24, 31, 25, 37, 35, 25, 40, 21, 7, 17, 22, 75},
    	vector<int>{20, 33, 17, 45, 12, 21, 20, 2, 7, 17, 21, 11, 11, 9, 21, 65},
    	vector<int>{15, 17, 9, 11, 5, 5, 12, 21, 17, 10, 5, 13, 9, 7, 13, 55},
    };

    vector<int> c{36, 83, 59, 71, 43, 67, 23, 52, 93, 25, 67, 89, 60, 47, 64};
    vector<int> b{87, 75, 65, 55};

    solve_lp_relaxation(a, b, c);

    return 0;
}
