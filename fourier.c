#include <math.h>

#include "fourier.h"

void devolve_linha(double complex matrix[MAX_SIZE][MAX_SIZE], int width, double complex vector[], int linha){

    //percorre a largura da matrix
    for (int i=0; i<width; i++){
        vector[i] = matrix[linha][i];
    }
}

void devolve_coluna(double complex matrix[MAX_SIZE][MAX_SIZE], int height, double complex vector[], int coluna){

    //percorre a altura da matrix
    for (int i=0; i<height; i++){
        vector[i] = matrix[i][coluna];
    }
}

void nft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
    for (int k = 0; k < n; k++) {
        t[k] = 0;

        for (int j = 0; j < n; j++) {
            t[k] += s[j] * cexp(sign * 2 * PI * k * j * I / n);
        }
    }
}

void nft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
    nft(s, t, n, -1);
}

void nft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
    nft(t, s, n, 1);

    for (int k = 0; k < n; k++) {
        s[k] /= n;
    }
}

void fft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {

    //se n for igua a 1, k está entre 0 e n-1, ou seja k = 0, desta forma 
    //o expoente de euler é 0.
    if(n == 1){
        t[0] = s[0];
        return;
    }
    
    //declara vetores pares e ímpares
    double complex sp[n/2];
    double complex si[n/2];

    //declara as transformadas dos respectivos vetores
    double complex tp[n/2];
    double complex ti[n/2];

    int contador_sp = 0; //índice de sp
    int contador_si = 0; //indice de si

    //cria o vetor dos índices pares
    for(int i=0; i<n; i+=2){
        sp[contador_sp] = s[i];
        contador_sp += 1; 
    }

    //cria o vetor dos índices ímpares
    for(int i=1; i<n; i+=2){
        si[contador_si] = s[i];
        contador_si += 1; 
    }

    //transoformada de Fourier de sp
    //ao invés de chamar a função nft_foward, chama essa mesma função, 
    //e tem fé de que ela fará o calculo do nft_foward da maneira correta 
    fft(sp, tp, n/2, sign);

    //transformada de Fourier de si
    fft(si, ti, n/2, sign);

    for (int k = 0; k < n/2; k++) {
        t[k] = tp[k] + ti[k]*cexp(sign * 2 * PI * k * I / n);
        t[k + n/2] = tp[k] - ti[k]*cexp(sign * 2 * PI * k * I / n);
    }
}

void fft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
    fft(s, t, n, -1);
}

void fft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
    fft(t, s, n, 1);

    for (int k = 0; k < n; k++) {
        s[k] /= n;
    }
}

void fft_forward_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {

    //para cada linha da matrix, aplicaremos
    //a transformada unidimensional sobre l
    for (int l = 0; l < height; l++){

        double complex t_linha[width]; //cria o vetor t
        double complex vector_linha[width]; //cria o vector que receberá a linha

        //altera o vetor para receber a linha da matrix desejada
        devolve_linha(matrix, width, vector_linha, l);

        // t = transformada de Fourier matrix[l]
        fft_forward(vector_linha, t_linha, width);

        //para cada coluna da matrix_t, a coluna recebe o mesmo valor de t
        for (int coluna = 0; coluna < width; coluna++){
            matrix[l][coluna] = t_linha[coluna];
        }
    }

    //para cada coluna da matrix, aplicaremos
    //a transformada unidimensional sobre c
    for (int c = 0; c < width; c++){

        double complex t_coluna[height]; //cria o vetor t
        double complex vector_coluna[height]; //cria o vector que receberá a coluna

        //altera o vetor para receber a coluna da matrix desejada
        devolve_coluna(matrix, height, vector_coluna, c);

        // t = transformada de Fourier matrix[c]
        fft_forward(vector_coluna, t_coluna, height);

        //para cada linha da matrix_t, a linha recebe o mesmo valor de t
        for (int linha = 0; linha < height; linha++){
            matrix[linha][c] = t_coluna[linha];
        }
    }
}

void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    
    //para cada linha da matrix, aplicaremos
    //a transformada unidimensional sobre l
    for (int l = 0; l < height; l++){
        
        double complex t_linha[width]; //cria o vetor t
        double complex vector_linha[width]; //cria o vector que receberá a linha

        //altera o vetor para receber a linha da matrix desejada
        devolve_linha(matrix, width, vector_linha, l);

        // t = transformada de Fourier matrix[l]
        fft_inverse(vector_linha, t_linha, width);

        //para cada coluna da matrix_t, a coluna recebe o mesmo valor de t
        for (int coluna = 0; coluna < width; coluna++){
            matrix[l][coluna] = t_linha[coluna];
        }
    }

    //para cada coluna da matrix, aplicaremos
    //a transformada unidimensional sobre c
    for (int c = 0; c < width; c++){

        double complex t_coluna[height]; //cria o vetor t        
        double complex vector_coluna[height]; //cria o vector que receberá a coluna

        //altera o vetor para receber a coluna da matrix desejada
        devolve_coluna(matrix, height, vector_coluna, c);

        // t = transformada de Fourier matrix[c]
        fft_inverse(vector_coluna, t_coluna, height);

        //para cada linha da matrix_t, a linha recebe o mesmo valor de t
        for (int linha = 0; linha < height; linha++){
            matrix[linha][c] = t_coluna[linha];
        }
    }
}

void filter(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;

            double g = exp(d / variance);

            if (flip) {
                g = 1 - g;
            }

            output[y][x] = g * input[y][x];
        }
    }
}

void filter_lp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 0);
}

void filter_hp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 1);
}
