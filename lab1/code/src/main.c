#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef double (*MathFunc)(double);

void nodes(double a, double b, int n, MathFunc f, double **x_h, double **y_h) {
  *x_h = (double *)malloc(n * sizeof(double));
  *y_h = (double *)malloc(n * sizeof(double));

  for (int i = 0; i < n; i++) {
    double cos_arg = (2 * i + 1) * M_PI / (2 * n);
    (*x_h)[i] = (a + b) / 2 + (b - a) / 2 * cos(cos_arg);
  }

  for (int i = 0; i < n; i++) {
    (*y_h)[i] = f((*x_h)[i]);
  }
}

double *divided_differences(double *x, double *y, int n) {
  double *dd = (double *)malloc(n * sizeof(double));
  for (int i = 0; i < n; i++)
    dd[i] = y[i];

  for (int j = 1; j < n; j++) {
    for (int i = n - 1; i >= j; i--) {
      dd[i] = (dd[i] - dd[i - 1]) / (x[i] - x[i - j]);
    }
  }
  return dd;
}

double newton_eval(double x_val, double *x, double *coeffs, int n) {
  double result = coeffs[n - 1];
  for (int i = n - 2; i >= 0; i--) {
    result = result * (x_val - x[i]) + coeffs[i];
  }
  return result;
}

void check_accuracy(double *x, double *y, double *coeffs, int n) {
  printf("Проверка в узлах:\n");
  double max_err = 0.0;
  for (int i = 0; i < n; i++) {
    double p_val = newton_eval(x[i], x, coeffs, n);
    double err = fabs(p_val - y[i]);
    if (err > max_err)
      max_err = err;
    printf("Узел %d: f=%.6f P=%.6f Ошибка=%.2e\n", i, y[i], p_val, err);
  }
  printf("Максимальная ошибка: %.2e\n", max_err);
}

double max_error(double a, double b, double *x, double *coeffs, int n,
                 MathFunc f) {
  double max_err = 0.0;
  int test_points = 1000;
  double step = (b - a) / test_points;

  for (double t = a; t <= b; t += step) {
    // Пропускаем узлы интерполяции
    int is_node = 0;
    for (int i = 0; i < n; i++) {
      if (fabs(t - x[i]) < 1e-9) {
        is_node = 1;
        break;
      }
    }
    if (is_node)
      continue;

    double err = fabs(f(t) - newton_eval(t, x, coeffs, n));
    if (err > max_err)
      max_err = err;
  }
  return max_err;
}

double sample_func(double x) { return x * x + 1 - acos(x); }

int main() {
  double a = -0.6, b = 0.2;
  int n = 8;

  double *x, *y;
  nodes(a, b, n, sample_func, &x, &y);

  double *coeffs = divided_differences(x, y, n);

  check_accuracy(x, y, coeffs, n);

  FILE *fout = fopen("graph_data.txt", "w");
  for (double i = a; i <= b + 0.02; i += 0.02) {
    double p_val = newton_eval(i, x, coeffs, n);
    fprintf(fout, "%lf %.8lf\n", i,
            sample_func(i) - p_val);
  }
  fclose(fout);

  free(x);
  free(y);
  free(coeffs);

  int max_nodes = 80;

  FILE *ferr = fopen("error_vs_nodes.txt", "w");
  fprintf(ferr, "Nodes MaxError\n");

  for (int n = 2; n <= max_nodes; n++) {
    double *x, *y;
    nodes(a, b, n, sample_func, &x, &y);

    double *coeffs = divided_differences(x, y, n);

    double err = max_error(a, b, x, coeffs, n, sample_func);
    printf("n = %d, max error = %.6e\n", n, err);
    fprintf(ferr, "%d %.6e\n", n, err);

    free(x);
    free(y);
    free(coeffs);
  }

  fclose(ferr);

  return 0;
}
