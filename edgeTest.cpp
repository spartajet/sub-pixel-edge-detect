#include "edgeTest.h"
using namespace std;

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif
#ifndef M_PI
#endif
void error(char *msg) {
    fprintf(stderr, "error: %s\n", msg);
    exit(EXIT_FAILURE);
}

void *xmalloc(size_t size) {
    void *p;
    if (size == 0) error("xmalloc: zero size");
    p = malloc(size);
    if (p == NULL) error("xmalloc: out of memory");
    return p;
}

int greater_round(double a, double b) {
    if (a <= b) return FALSE;
    if ((a - b) < 1000 * DBL_EPSILON) return FALSE;
    return TRUE;
}

double dist(double x1, double y1, double x2, double y2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

void gaussian_kernel(double *pfGaussFilter, int n, double sigma, double mean) {
    double sum = 0.0;
    double val;
    int i;
    if (pfGaussFilter == NULL) error("gaussian_kernel: kernel not allocated");
    if (sigma <= 0.0) error("gaussian_kernel: sigma must be positive");

    // compute Gaussian kernel
    for (i = 0; i < n; i++) {
        val = ((double) i - mean) / sigma;
        pfGaussFilter[i] = exp(-0.5 * val * val);
        sum += pfGaussFilter[i];
    }

    // normalization
    if (sum > 0.0)
        for (i = 0; i < n; i++)
            pfGaussFilter[i] /= sum;
}

void gaussian_filter(uchar *image, uchar *out, int X, int Y, double sigma) {
    int x, y, offset, i, j, nx2, ny2, n;
    double *kernel;
    double *tmp;
    double val, prec;

    if (sigma <= 0.0) error("gaussian_filter: sigma must be positive");
    if (image == NULL || X < 1 || Y < 1) error("gaussian_filter: invalid image");


    tmp = (double *) xmalloc(X * Y * sizeof(double));
    prec = 3.0;
    offset = (int) ceil(sigma * sqrt(2.0 * prec * log(10.0)));
    n = 1 + 2 * offset; //kernel size
    kernel = (double *) xmalloc(n * sizeof(double));
    gaussian_kernel(kernel, n, sigma, (double) offset);
    // for (int i = 0; i < n; i++) {
    //     // printf("%d--%f\n", i, kernel[i]);
    // }
    // auxiliary variables for the double of the image size
    nx2 = 2 * X;
    ny2 = 2 * Y;

    // x axis convolution
    for (x = 0; x < X; x++)
        for (y = 0; y < Y; y++) {
            val = 0.0;
            for (i = 0; i < n; i++) {
                j = x - offset + i;

                // symmetry boundary condition
                while (j < 0)
                    j += nx2;
                while (j >= nx2)
                    j -= nx2;
                if (j >= X)
                    j = nx2 - 1 - j;

                val += (double) image[j + y * X] * kernel[i];
            }
            tmp[x + y * X] = val;
        }

    // y axis convolution
    for (x = 0; x < X; x++)
        for (y = 0; y < Y; y++) {
            val = 0.0;
            for (i = 0; i < n; i++) {
                j = y - offset + i;

                // symmetry boundary condition
                while (j < 0)
                    j += ny2;
                while (j >= ny2)
                    j -= ny2;
                if (j >= Y)
                    j = ny2 - 1 - j;

                val += tmp[x + j * X] * kernel[i];
            }
            out[x + y * X] = (uchar) val;
        }

    free((void *) kernel);
    free((void *) tmp);
}

double chain(int from, int to, double *Ex, double *Ey, double *Gx, double *Gy, int X, int Y) {
    double dx, dy;
    if (Ex == NULL || Ey == NULL || Gx == NULL || Gy == NULL)
        error("chain: invalid input");
    if (from < 0 || to < 0 || from >= X * Y || to >= X * Y)
        error("chain: one of the points is out the image");

    //check that the points are different and valid edge points,otherwise return invalid chaining
    if (from == to)
        return 0.0; // same pixel, not a valid chaining
    if (Ex[from] < 0.0 || Ey[from] < 0.0 || Ex[to] < 0.0 || Ey[to] < 0.0)
        return 0.0; // one of them is not an edge point, not a valid chaining

    dx = Ex[to] - Ex[from];
    dy = Ey[to] - Ey[from];
    if ((Gy[from] * dx - Gx[from] * dy) * (Gy[to] * dx - Gx[to] * dy) <= 0.0)
        return 0.0;
    if ((Gy[from] * dx - Gx[from] * dy) >= 0.0)
        return 1.0 / dist(Ex[from], Ey[from], Ex[to], Ey[to]);
    else
        return -1.0 / dist(Ex[from], Ey[from], Ex[to], Ey[to]);
}

void compute_gradient(double *Gx, double *Gy, double *modG, uchar *image, int X, int Y) {
    int x, y;

    if (Gx == NULL || Gy == NULL || modG == NULL || image == NULL)
        error("compute_gradient: invalid input");

    // approximate image gradient using centered differences
    for (x = 1; x < (X - 1); x++)
        for (y = 1; y < (Y - 1); y++) {
            Gx[x + y * X] = (double) image[(x + 1) + y * X] - (double) image[(x - 1) + y * X];
            Gy[x + y * X] = (double) image[x + (y + 1) * X] - (double) image[x + (y - 1) * X];
            modG[x + y * X] = sqrt(Gx[x + y * X] * Gx[x + y * X] + Gy[x + y * X] * Gy[x + y * X]);
        }
}

void compute_edge_points(double *Ex, double *Ey, double *modG,
                         double *Gx, double *Gy, int X, int Y) {
    int x, y, i;

    if (Ex == NULL || Ey == NULL || modG == NULL || Gx == NULL || Gy == NULL)
        error("compute_edge_points: invalid input");
    for (i = 0; i < X * Y; i++) Ex[i] = Ey[i] = -1.0;
    for (x = 2; x < (X - 2); x++)
        for (y = 2; y < (Y - 2); y++) {
            int Dx = 0;
            int Dy = 0;
            double mod = modG[x + y * X];
            double L = modG[x - 1 + y * X];
            double R = modG[x + 1 + y * X];
            double U = modG[x + (y + 1) * X];
            double D = modG[x + (y - 1) * X];
            double gx = fabs(Gx[x + y * X]);
            double gy = fabs(Gy[x + y * X]);
            if (greater_round(mod, L) && !greater_round(R, mod) && gx >= gy) {
                Dx = 1;
            } else if (greater_round(mod, D) && !greater_round(U, mod) && gx <= gy) {
                Dy = 1;
            }
            if (Dx > 0 || Dy > 0) {
                double a = modG[x - Dx + (y - Dy) * X];
                double b = modG[x + y * X];
                double c = modG[x + Dx + (y + Dy) * X];
                double offset = 0.5 * (a - c) / (a - b - b + c);

                Ex[x + y * X] = x + offset * Dx;
                Ey[x + y * X] = y + offset * Dy;
            }
        }
}

void chain_edge_points(int *next, int *prev, double *Ex, double *Ey, double *Gx, double *Gy, int X, int Y) {
    int x, y, i, j, alt;

    if (next == NULL || prev == NULL || Ex == NULL || Ey == NULL || Gx == NULL || Gy == NULL)
        error("chain_edge_points: invalid input");

    for (i = 0; i < X * Y; i++) next[i] = prev[i] = -1;

    for (x = 2; x < (X - 2); x++)
        for (y = 2; y < (Y - 2); y++)
            if (Ex[x + y * X] >= 0.0 && Ey[x + y * X] >= 0.0) {
                int from = x + y * X;
                double fwd_s = 0.0;
                double bck_s = 0.0;
                int fwd = -1;
                int bck = -1;
                for (i = -2; i <= 2; i++)
                    for (j = -2; j <= 2; j++) {
                        int to = x + i + (y + j) * X;
                        double s = chain(from, to, Ex, Ey, Gx, Gy, X, Y);
                        if (s > fwd_s) {
                            fwd_s = s;
                            fwd = to;
                        }
                        if (s < bck_s) {
                            bck_s = s;
                            bck = to;
                        }
                    }

                if (fwd >= 0 && next[from] != fwd &&
                    ((alt = prev[fwd]) < 0 || chain(alt, fwd, Ex, Ey, Gx, Gy, X, Y) < fwd_s)) {
                    if (next[from] >= 0) prev[next[from]] = -1;
                    next[from] = fwd;
                    if (alt >= 0) next[alt] = -1;
                    prev[fwd] = from;
                }
                if (bck >= 0 && prev[from] != bck &&
                    ((alt = next[bck]) < 0 || chain(alt, bck, Ex, Ey, Gx, Gy, X, Y) > bck_s)) {
                    if (alt >= 0) prev[alt] = -1;
                    next[bck] = from;
                    if (prev[from] >= 0) next[prev[from]] = -1;
                    prev[from] = bck;
                }
            }
}

void thresholds_with_hysteresis(int *next, int *prev,
                                double *modG, int X, int Y, double th_h, double th_l) {
    int *valid;
    int i, j, k;

    if (next == NULL || prev == NULL || modG == NULL)
        error("thresholds_with_hysteresis: invalid input");

    valid = (int *) xmalloc(X * Y * sizeof(int));
    for (i = 0; i < X * Y; i++) valid[i] = FALSE;

    for (i = 0; i < X * Y; i++)
        if ((prev[i] >= 0 || next[i] >= 0) && !valid[i] && modG[i] >= th_h) {
            valid[i] = TRUE;
            for (j = i; j >= 0 && (k = next[j]) >= 0 && !valid[k]; j = next[j])
                if (modG[k] < th_l) {
                    next[j] = -1;
                    prev[k] = -1;
                } else
                    valid[k] = TRUE;
            for (j = i; j >= 0 && (k = prev[j]) >= 0 && !valid[k]; j = prev[j])
                if (modG[k] < th_l) {
                    prev[j] = -1;
                    next[k] = -1;
                } else
                    valid[k] = TRUE;
        }

    for (i = 0; i < X * Y; i++)
        if ((prev[i] >= 0 || next[i] >= 0) && !valid[i])
            prev[i] = next[i] = -1;

    free((void *) valid);
}

void list_chained_edge_points(double **x, double **y, int *N, int **curve_limits,
                              int *M, int *next, int *prev, double *Ex, double *Ey, int X, int Y) {
    int i, k, n;

    *x = (double *) xmalloc(X * Y * sizeof(double));
    *y = (double *) xmalloc(X * Y * sizeof(double));
    *curve_limits = (int *) xmalloc(X * Y * sizeof(int));
    *N = 0;
    *M = 0;

    for (i = 0; i < X * Y; i++)
        if (prev[i] >= 0 || next[i] >= 0) {
            (*curve_limits)[*M] = *N;
            ++(*M);

            for (k = i; (n = prev[k]) >= 0 && n != i; k = n);

            do {
                (*x)[*N] = Ex[k];
                (*y)[*N] = Ey[k];
                ++(*N);

                n = next[k];
                next[k] = -1;
                prev[k] = -1;


                k = n;
            } while (k >= 0);
        }
    (*curve_limits)[*M] = *N;
}

void devernay(double **x, double **y, int *N, int **curve_limits, int *M,
              uchar *image, uchar *gauss, int X, int Y, double sigma, double th_h, double th_l) {
    double *Gx = (double *) xmalloc(X * Y * sizeof(double));
    double *Gy = (double *) xmalloc(X * Y * sizeof(double));
    double *modG = (double *) xmalloc(X * Y * sizeof(double));
    double *Ex = (double *) xmalloc(X * Y * sizeof(double));
    double *Ey = (double *) xmalloc(X * Y * sizeof(double));
    int *next = (int *) xmalloc(X * Y * sizeof(int));
    int *prev = (int *) xmalloc(X * Y * sizeof(int));
    if (sigma == 0.0) compute_gradient(Gx, Gy, modG, image, X, Y);
    else {
        gaussian_filter(image, gauss, X, Y, sigma);
        compute_gradient(Gx, Gy, modG, gauss, X, Y);
    }
    // for (int i = 0; i < X * Y; i++) {
    //     if (gauss[i] < 255) {
    //         // printf("%d %d\n", i, gauss[i]);
    //     }
    // }
    for (int i = 0; i < X * Y; i++) {
        if (Gx[i] == 0 || Gy[i] == 0 || modG[i] == 0) {
            continue;
        }
    }
    compute_edge_points(Ex, Ey, modG, Gx, Gy, X, Y);

    chain_edge_points(next, prev, Ex, Ey, Gx, Gy, X, Y);

    thresholds_with_hysteresis(next, prev, modG, X, Y, th_h, th_l);

    list_chained_edge_points(x, y, N, curve_limits, M, next, prev, Ex, Ey, X, Y);

    free((void *) Gx);
    free((void *) Gy);
    free((void *) modG);
    free((void *) Ex);
    free((void *) Ey);
    free((void *) next);
    free((void *) prev);
}
