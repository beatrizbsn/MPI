// MPI.cpp 

#include <iostream>
#include <vector>
#include <mpi.h>
#include <iomanip>

#define Tamanho 50
#define Ngeracao 4


float getNeighbords_media(const std::vector<std::vector<float>>& grid, int i, int j, const std::vector<float>& Down , const std::vector<float>& Top) {
    float vizinhosVivos = 0.0;
    int jmenos1, jmais1, imenos1, imais1;

    jmais1 = j + 1;
    jmenos1 = j - 1;
    imais1 = i + 1;
    imenos1 = i - 1;

    // A  | B  | C
    // D  |i,j | E
    // F  |G   | H

    if (jmais1 > Tamanho - 1) {
        jmais1 = 0;
    }
    if (jmenos1 < 0) {
        jmenos1 = Tamanho - 1;
    }


    if (imais1 > Tamanho - 1) {
        vizinhosVivos = Down [jmenos1] + Down [j] + Down [jmais1];
    }
    else {
        vizinhosVivos = grid[imais1][jmenos1] + grid[imais1][j] + grid[imais1][jmais1];
    }


    if (imenos1 < 0) {
        vizinhosVivos = vizinhosVivos + Top[jmenos1] + Top[j] + Top[jmais1];
    }
    else {
        vizinhosVivos = vizinhosVivos + grid[imenos1][jmenos1] + grid[imenos1][j] + grid[imenos1][jmais1];
    }

    vizinhosVivos = vizinhosVivos + grid[i][jmenos1] + grid[i][jmais1];

    return vizinhosVivos / 2.00;
}


int getNeighbords(const std::vector<std::vector<float>>& grid, int i, int j,const std::vector<float>& Down , const std::vector<float>& Top) {
    int vizinhosVivos = 0;
    int jmenos1, jmais1, imenos1, imais1;

    jmais1 = j + 1;
    jmenos1 = j - 1;
    imais1 = i + 1;
    imenos1 = i - 1;

    // A  | B  | C
    // D  |i,j | E
    // F  |G   | H

    if (jmais1 > Tamanho - 1) {
        jmais1 = 0;
    }
    if (jmenos1 < 0) {
        jmenos1 = Tamanho - 1;
    }


    if (imais1 > Tamanho - 1) {
        if (Down [jmenos1] > 0) { vizinhosVivos = vizinhosVivos + 1; }
        if (Down [j] > 0) { vizinhosVivos = vizinhosVivos + 1; }
        if (Down [jmais1] > 0) { vizinhosVivos = vizinhosVivos + 1; }
    }
    else {
        if (grid[imais1][jmenos1] > 0) { vizinhosVivos = vizinhosVivos + 1; }
        if (grid[imais1][j] > 0) { vizinhosVivos = vizinhosVivos + 1; }
        if (grid[imais1][jmais1] > 0) { vizinhosVivos = vizinhosVivos + 1; }
    }


    if (imenos1 < 0) {
        if (Top[jmenos1] > 0) { vizinhosVivos = vizinhosVivos + 1; }
        if (Top[j] > 0) { vizinhosVivos = vizinhosVivos + 1; }
        if (Top[jmais1] > 0) { vizinhosVivos = vizinhosVivos + 1; }
    }
    else {
        if (grid[imenos1][jmenos1] > 0) { vizinhosVivos = vizinhosVivos + 1; } 
        if (grid[imenos1][j] > 0) { vizinhosVivos = vizinhosVivos + 1; } 
        if (grid[imenos1][jmais1] > 0) { vizinhosVivos = vizinhosVivos + 1; }

    }

    
    if (grid[i][jmenos1] > 0) { vizinhosVivos = vizinhosVivos + 1; }
    if (grid[i][jmais1] > 0) { vizinhosVivos = vizinhosVivos + 1; }
   

    return vizinhosVivos;
}



int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int Id, Proc;
    int soma = 0;

    double start_time = MPI_Wtime();


    MPI_Comm_rank(MPI_COMM_WORLD, &Id);
    MPI_Comm_size(MPI_COMM_WORLD, &Proc);

    int Linha = (Tamanho / Proc) + (Tamanho % Proc > Id);
    std::vector<std::vector<float>> grid(Linha, std::vector<float>(Tamanho, 0.0));
    std::vector<std::vector<float>> newgrid(Linha, std::vector<float>(Tamanho, 0.0));
    std::vector<float> Top(Tamanho, 0.0);
    std::vector<float> Down (Tamanho, 0.0);

    if (Id == 0) {

        //GLIDER
        grid[1][2] = 1.0;
        grid[2][3] = 1.0;
        grid[3][1] = 1.0;
        grid[3][2] = 1.0;
        grid[3][3] = 1.0;

        //R-pentomino
        grid[10][31] = 1.0;
        grid[10][32] = 1.0;
        grid[11][30] = 1.0;
        grid[11][31] = 1.0;
        grid[12][31] = 1.0;
    }
   

    MPI_Barrier(MPI_COMM_WORLD);
    
    
    for (int geracao = 0; geracao < Ngeracao; geracao++) {

     
        int top = (Id + 1) % Proc;
        int down = (Id - 1 + Proc) % Proc;
       

        MPI_Status status;
        MPI_Sendrecv(&grid[0][Tamanho - 1], Tamanho, MPI_FLOAT,top, 0,
                         &Top, Tamanho, MPI_FLOAT, top, 0, MPI_COMM_WORLD, &status);
       
        MPI_Status status1;
        MPI_Sendrecv(&grid[Linha - 1][Tamanho - 1], Tamanho, MPI_FLOAT, down, 0,
                              &Down , Tamanho, MPI_FLOAT, down, 0, MPI_COMM_WORLD, &status1);
       
       

       
        for (int i = 1; i < Linha - 1; i++) {
            for (int j = 0; j < Tamanho; j++) {

                int vizinhosVivos = getNeighbords(grid, i,j, Down  , Top);

                // vivas
                if (grid[i][j] > 0) {
                    if (vizinhosVivos == 2 || vizinhosVivos == 3) {
                        newgrid[i][j] = grid[i][j];
                        
                    }
                    else {
                        newgrid[i][j] = 0.0;
                    }
                    // morta
                }
                else {

                    if (vizinhosVivos == 3) {
                        newgrid[i][j] = getNeighbords_media(grid, i, j, Down , Top);
                       
                    }
                    else {
                        newgrid[i][j] = 0.0;
                    }
                }
            }
        }

        
        std::swap(grid, newgrid);

        MPI_Barrier(MPI_COMM_WORLD);      
    }

  

    MPI_Barrier(MPI_COMM_WORLD);
    int sum = 0;
    for (int n = 1; n < Linha - 1; n++) {
        for (int m = 0; m < Tamanho ; m++) {
            sum = sum + grid[n][m];
        }
    }

    MPI_Request request;
    if (Id == 0) {
        int sumAux;

        for (int k = 1; k < Proc ; k++) {
            MPI_Recv(&sumAux, 1, MPI_INT, k, 10, MPI_COMM_WORLD, NULL);
            sum += sumAux;
        }

        std::cout << "Final: " << sum;
        double end_time = MPI_Wtime();
        std::cout << "Tempo decorrido: " << end_time - start_time << " segundos" << std::endl;
    }
    else {
        MPI_Isend(&sum, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &request);
    }

    

    MPI_Finalize();


    
    
    return 0;
}