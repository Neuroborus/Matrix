enum WhatToDo{SumRws, SwapRws};
struct HistB{    // история преобразований для столбца свободных членов
    int* firstRow;
    int* secondRow;
    double* multiplier;
    WhatToDo* wtd;
    int counter;

    void SumRows(int fR, int sR, double m = 1){ //Строка(firstRow), умноженная на multiplier, прибавляется к secondRow, если wtd = SumRows
        firstRow[counter] = fR;
        secondRow[counter] = sR;
        multiplier[counter] = m;
        wtd[counter] = SumRws;
        counter++;
    }
    void SwapRows(int fR, int sR){  //если wtd = SwapRows, свапает строки соответственно
        firstRow[counter] = fR;
        secondRow[counter] = sR;
        wtd[counter] = SwapRws;
        counter++;
    }

    HistB(int n){
        firstRow = new int[n]{};
        secondRow = new int[n]{};
        multiplier = new double[n]{};
        wtd = new WhatToDo[n]{};
        counter = 0;
    }
    ~HistB(){
        delete[] firstRow;
        delete[] secondRow;
        delete[] multiplier;
        delete[] wtd;
    }
};