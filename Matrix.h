#pragma once
#include "MSize.h"
#include "HistB.h"

#include <iostream>
#include <cmath>
using std::cout;
using std::cerr;
using std::swap;
using std::runtime_error;

template<typename T>
class Matrix
{
private:
	T **arr;    //Матрица эл-тов
	MSize size; //Структура размерности матрицы
	bool signTumbler;   //Флаг, нужно ли менять знак у определителя матрицы (учитывает перестановку строк)
	bool ReCreateMatrix(MSize); //Приватный метод для создания/пересоздания матрицы

public:
	Matrix();
    explicit Matrix(MSize);
	Matrix(unsigned int, unsigned int);
    Matrix(const Matrix&);
	virtual ~Matrix();

    Matrix UTriangular();   //Метод приведения матрицы к верхнетругольному виду. Не изменяет саму матрицу, возвращает измененную
    Matrix Transpose();   //Метод транспонирования марицы. Не изменяет саму матрицу, возвращает измененную
	T Determinant();    //Метод возвращает детерминант матрицы

    Matrix InvertSOLE();    //Метод нахождения обратной матрицы через решение систем линейных уравнений. Саму матрицу не изменяет, возвращает измененную.
	Matrix InvertNaive();   //Метод нахождения обратной матрицы с помощью матрицы алгебраических дополнений. Саму матрицу не изменяет, возвращает измененную.

	void Fill();    //Метод заполнения матрицы вручную
	void Show();
	MSize GetSize() {return size;}
	////////////////////////////////////////////////
	Matrix& operator=(const Matrix<T>&);
	Matrix& operator+=(const Matrix&);
	Matrix& operator+=(const T&);
	Matrix& operator-=(const Matrix&);
	Matrix& operator-=(const T&);
	//
	Matrix operator*(const Matrix);
	Matrix operator*(const T);
	Matrix operator/(const Matrix);
	Matrix operator/(const T t);
};
/////////////////////////////////////////
/////////////////////////////////////////
template<typename T>
Matrix<T> operator+(const Matrix<T>& lm, const Matrix<T>& rm);
template<typename T>
Matrix<T> operator+(const Matrix<T>& lm, const T& rt);
template<typename T>
Matrix<T> operator-(const Matrix<T>& lm, const Matrix<T>& rm);
template<typename T>
Matrix<T> operator-(const Matrix<T>& lm, const T &rt);
//
template<typename T>
Matrix<T>::Matrix(): signTumbler(false) {
}
template<typename T>
Matrix<T>::Matrix(MSize ms): signTumbler(false) {
    if (ReCreateMatrix(ms))
        cout << "Matrix recreated successfully!\n";
    else
        cout << "Matrix created successfully!\n";
}
template<typename T>
Matrix<T>::Matrix(unsigned int nr, unsigned int nc): signTumbler(false) {
    MSize tms(nr, nc);
    if (ReCreateMatrix(tms))
        cout << "Matrix recreated successfully!\n";
    else
        cout << "Matrix created successfully!\n";
}
template<typename T>
Matrix<T>::Matrix(const Matrix &other): signTumbler(false){
    if (this->size != other.size) {
        if (ReCreateMatrix(other.size))
            cout << "Matrix recreated successfully!\n";
        else
            cout << "Matrix created successfully!\n";
    }
    for (int i = 0; i < size.nrow; i++)
        for (int j = 0; j < size.ncol; j++)
            this->arr[i][j] = other.arr[i][j];
    this->signTumbler=other.signTumbler;
    cout << "Copying ctor done!\n";
}
template<typename T>
Matrix<T>::~Matrix(){
    if (size) {
        for (int i = 0; i < size.nrow; i++)
            if (arr[i])
                delete[] arr[i];
        if (arr)
            delete[] arr;
    }
}
//
template<typename T>
bool Matrix<T>::ReCreateMatrix(MSize tms) {
    bool flag = false;
    if (size) {
        flag = true;
        for (int i = 0; i < size.nrow; i++)
            if (arr[i])
                delete[] arr[i];
        if (arr)
            delete[] arr;
    }
    arr = new T *[tms.nrow];
    for (int i = 0; i < tms.nrow; i++)
        arr[i] = new T[tms.ncol];
    size = tms;

    return flag;
}
template<typename T>
Matrix<T> Matrix<T>::UTriangular(){
    if (!size)
        throw runtime_error("Matrix does not exists!");
    if (!size.IsSquare())
        throw runtime_error("Error! Non-square matrix!");
    int i, j;
    Matrix<double> tmp(size);
    for(i = 0; i < size.nrow; i++)
        for(j = 0; j < size.ncol; j++)
            tmp.arr[i][j] = arr[i][j];
    int im = 0, jm = 0;
    int loopCounter = 0;
    bool flag;
    double var;
    while(im < size.nrow-1 && jm < size.ncol-1) {
        //1. В первом столбце выбрать элемент, отличный от нуля
        // ведущий элемент). Строку с ведущим элементом
        // (ведущая строка), если она не первая, переставить на место
        // первой строки (преобразование I типа). Если в первом
        // столбце нет ведущего (все элементы равны нулю), то
        // исключаем этот столбец, и продолжаем поиск ведущего
        // элемента в оставшейся части матрицы. Преобразования
        // заканчиваются, если исключены все столбцы или в оставшейся
        // части матрицы все элементы нулевые.
        flag = false;
        for (j = jm; j < tmp.size.ncol; j++) {      //Найдем отличный от нуля элемент
            for (i = im; i < tmp.size.nrow; i++)
                if (tmp.arr[i][j]) {
                    im = i;
                    jm = j;
                    flag = true;
                    break;
                }
            if (flag)
                break;
        }
        if(im!=loopCounter) {   //Переместим строку вверх, если она не вверху
            swap(tmp.arr[im], tmp.arr[loopCounter]);
            im = loopCounter;
            tmp.signTumbler=!tmp.signTumbler;   //У определителя меняется знак
        }
        if (!flag) {  //[под]Матрица полностью состоит из нулей
            cout << "[sub-]Matrix consists of zeroes!\n";
            return tmp;
        }
        //2. Разделить все элементы ведущей строки
        // на ведущий элемент (преобразование II типа). Если ведущая строка последняя,
        // то на этом преобразования следует закончить.
        Matrix<double> st(1, size.ncol-jm);
        for(j = 0; j < size.ncol - jm; j++)
            st.arr[0][j] = tmp.arr[im][j + jm];
        for (j = 0; j < tmp.size.ncol - jm; j++)
            st.arr[0][j] /= tmp.arr[im][jm];
        //3. К каждой строке, расположенной ниже ведущей, прибавить ведущую строку,
        // умноженную соответственно на такое число, чтобы элементы, стоящие под ведущим
        // оказались равными нулю (преобразование III типа).
        for(i = im+1; i<tmp.size.nrow; i++) {
            var = 0;
            for (j = 0; j < st.size.ncol; j++) {
                if (!var)
                    var = (-tmp.arr[i][j + jm]) / st.arr[0][j];
                tmp.arr[i][j + jm] += var * st.arr[0][j];
            }
        }
        //4. Исключив из рассмотрения строку и столбец, на пересечении которых стоит ведущий элемент,
        // перейти к пункту 1, в котором все описанные действия применяются к оставшейся части матрицы.
        im++;
        jm++;
        loopCounter++;
    }
    return tmp;
}
template<typename T>
Matrix<T> Matrix<T>::Transpose() {
    if (!size)
        throw runtime_error("Matrix does not exists!");
    Matrix tmp(size.ncol,size.nrow);
    for (int i = 0; i < tmp.size.nrow; i++)
        for (int j = 0; j < tmp.size.ncol; j++)
            tmp.arr[i][j] = arr[j][i];
    return tmp;
}
template<typename T>
T Matrix<T>::Determinant() {
    if (!size)
        throw runtime_error("Matrix does not exists!");
    if (!size.IsSquare())
        throw runtime_error("Error! Non-square matrix!");
    T det;

    if(size.ncol==1 && size.nrow==1)
        return arr[0][0];
    if(size.ncol==2 && size.nrow==2)
        return arr[0][0]*arr[1][1] - arr[0][1]*arr[1][0];

    Matrix tmp = this->UTriangular();
    det = arr[0][0];
    for(int i = 1; i<size.nrow; i++) {
        if(tmp.arr[i][i])
            det *= tmp.arr[i][i];
        else
            return 0;
    }
    if(signTumbler && det)
        det=-det;
    return det;
}
template<typename T>
Matrix<T> Matrix<T>::InvertSOLE() {    //Нахождение обратной матрицы с помощью решения систем линейных алгебраических уравнений
    if (!size)
        throw runtime_error("Matrix does not exists!");
    if (!size.IsSquare())
        throw runtime_error("Error! Non-square matrix!");
    int i, j, k;
    int im, jm;
    int loopCounter;
    bool flag;
    double var;
    HistB hist(size.nrow*size.nrow);  //Объект для записи истории преобразований, чтобы применить потом к столбцу свободных членов
    Matrix<double> tmp(size);
    for (i = 0; i < size.nrow; i++)
        for (j = 0; j < size.ncol; j++)
            tmp.arr[i][j] = arr[i][j];
    //////////////////////////////
    im = 0, jm = 0;
    loopCounter = 0;
    while (im < size.nrow - 1 && jm < size.ncol - 1) {
        flag = false;
        for (j = jm; j < tmp.size.ncol; j++) {      //Найдем отличный от нуля элемент
            for (i = im; i < tmp.size.nrow; i++)
                if (tmp.arr[i][j]) {
                    im = i;
                    jm = j;
                    flag = true;
                    break;
                }
            if (flag)
                break;
        }
        if (im != loopCounter) {   //Переместим строку вверх, если она не вверху
            hist.SwapRows(im, loopCounter); //Записываем изменения в "журнал"
            swap(tmp.arr[im], tmp.arr[loopCounter]);
            im = loopCounter;
            tmp.signTumbler = !tmp.signTumbler;   //У определителя меняется знак
        }
        if (!flag) {  //[под]Матрица полностью состоит из нулей
            cout << "[sub-]Matrix consists of zeroes!\n";
            return tmp;
        }
        Matrix<double> st(1, size.ncol - jm);
        for (j = 0; j < size.ncol - jm; j++)
            st.arr[0][j] = tmp.arr[im][j + jm];
        for (j = 0; j < tmp.size.ncol - jm; j++)
            st.arr[0][j] /= tmp.arr[im][jm];
        for (i = im + 1; i < tmp.size.nrow; i++) {
            var = 0;
            for (j = 0; j < st.size.ncol; j++) {
                if (!var)
                    var = (-tmp.arr[i][j + jm]) / st.arr[0][j];
                tmp.arr[i][j + jm] += var * st.arr[0][j];
            }
            hist.SumRows(im, i, 1/tmp.arr[im][jm]*var);
        }
        im++;
        jm++;
        loopCounter++;
    }
    ///////////////////////////////
    Matrix<double> sol(size);   //Для записи решений слау
    double *b = new double[size.nrow];    //столбец свободных членов
    T sum;
    for (k = 0; k < size.nrow; k++) {
        //
        for (i = 0; i < size.nrow; i++)
            if (i == k)
                b[i] = 1.0;
            else
                b[i] = 0.0;
        for(i = 0; i < hist.counter; i++)
            switch(hist.wtd[i])
            {
                case SumRws:
                    b[hist.secondRow[i]]+=b[hist.firstRow[i]]*hist.multiplier[i];
                    break;
                case SwapRws:
                    swap(b[hist.firstRow[i]], b[hist.secondRow[i]]);
                    break;
                default:
                    throw runtime_error("Error! Unknown key!");
                    break;
            }
        for(i = size.nrow-1; i>=0; i--){
            sum = 0;
            for(j = size.nrow-1; j > i; j--)
                sum += sol.arr[j][k] * tmp.arr[i][j];
            if(tmp.arr[i][i])
                sol.arr[i][k] = (b[i]-sum)/tmp.arr[i][i];
            else
                throw runtime_error("Error! SOLE is incompatible!");    //Система несовместна
        }
    }
    ///////////////////////
    return sol;
}
template<typename T>
Matrix<T> Matrix<T>::InvertNaive() {    //Нахождение обратной матрицы с помощью матрицы алгебраических дополнений
    if (!size)
        throw runtime_error("Matrix does not exists!");
    if (!size.IsSquare())
        throw runtime_error("Error! Non-square matrix!");
    T t;
    Matrix tmp(size);
    Matrix fMinor(size.nrow - 1, size.ncol - 1);
    unsigned int skipRow, skipCol;
    for (int i = 0; i < size.nrow; i++) {
        for (int j = 0; j < size.ncol; j++) {
            skipRow = 0;
            for (int r = 0; r < size.nrow; r++) {
                if (r == i) {
                    skipRow++;
                    continue;
                }
                skipCol = 0;
                for (int c = 0; c < size.ncol; c++) {
                    if (c == j) {
                        skipCol++;
                        continue;
                    }
                    fMinor.arr[r - skipRow][c - skipCol] = this->arr[r][c];
                }
            }
            tmp.arr[i][j] = fMinor.Determinant() * pow(-1, i + j);    //Получим матрицу алгебраических дополнений
        }
    }
    t = Determinant();
    if(!t)
        throw runtime_error("Error! Det == 0!");
    else {
        tmp = tmp.Transpose();
        tmp = tmp / t;
    }
    return tmp;
}
template<typename T>
void Matrix<T>::Fill()
{
    if (!size)
        throw runtime_error("Matrix does not exists!");
    cout<<"M["<<size.nrow<<"]["<<size.ncol<<"]("<<(int)size<<")->";
    for (int i = 0; i < size.nrow; i++)
        for (int j = 0; j < size.ncol; j++)
            std::cin >> arr[i][j];
}
template<typename T>
void Matrix<T>::Show(){
    if (!size)
        throw runtime_error("Matrix does not exists!");
    cout << "M{\n";
    for (int i = 0; i < size.nrow; i++) {
        cout << '\t';
        for (int j = 0; j < size.ncol; j++) {
            cout << arr[i][j] << '\t';
        }
        cout << '\n';
    }
    cout << "}\n";
}

///////////////////////////
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) {
    if (this->size != other.size){
        if (ReCreateMatrix(other.size))
            cout << "Matrix recreated successfully!\n";
        else
            cout << "Matrix created successfully!\n";
    }
    for (int i = 0; i < size.nrow; i++)
        for (int j = 0; j < size.ncol; j++)
            this->arr[i][j] = other.arr[i][j];
    cout << "Assignment done!\n";
    return *this;
}
template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix& other)
{
    if (this->size == other.size)
        for (int i = 0; i < size.nrow; i++)
            for (int j = 0; j < size.ncol; j++)
                arr[i][j] += other.arr[i][j];
    else
        cerr << "\tError! Different sizes!\n";
    return *this;
}
template<typename T>
Matrix<T>& Matrix<T>::operator+=(const T& t)
{
    for (int i = 0; i < size.nrow; i++)
        for (int j = 0; j < size.ncol; j++)
            arr[i][j] += t;
    return *this;
}
template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix& other)
{
    if (this->size == other.size)
        for (int i = 0; i < size.nrow; i++)
            for (int j = 0; j < size.ncol; j++)
                arr[i][j] -= other.arr[i][j];
    else
        cerr << "\tError! Different sizes!\n";
    return *this;
}
template<typename T>
Matrix<T>& Matrix<T>::operator-=(const T& t)
{
    for (int i = 0; i < size.nrow; i++)
        for (int j = 0; j < size.ncol; j++)
            arr[i][j] -= t;
    return *this;
}
//
template<typename T>
Matrix<T> operator+(const Matrix<T>& lm, const Matrix<T>& rm)
{
    Matrix<T> tmp(lm);
    tmp += rm;
    return tmp;
}
template<typename T>
Matrix<T> operator+(const Matrix<T>& lm, const T& rt)
{
    Matrix<T> tmp(lm);
    tmp += rt;
    return tmp;
}
template<typename T>
Matrix<T> operator-(const Matrix<T>& lm, const Matrix<T>& rm)
{
    Matrix<T> tmp(lm);
    tmp -= rm;
    return tmp;
}
template<typename T>
Matrix<T> operator-(const Matrix<T>& lm, const T &rt) {
    Matrix<T> tmp(lm);
    tmp -= rt;
    return tmp;
}

//
template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> rm)
{
    if (size.ncol != rm.size.nrow)
        throw runtime_error("\tError! Incompatible sizes!\n");
    T t;
    Matrix tmp(size.nrow, rm.size.ncol);
    for (int r1 = 0; r1 < size.nrow; r1++)	//строка первой матрицы
        for (int c2 = 0; c2 < rm.size.ncol; c2++)	//столбец второй матрицы
        {
            t = 0;
            for (int c1r2 = 0; c1r2 < size.ncol; c1r2++)    //столбец первой и строка второй матрицы
                t += arr[r1][c1r2] * rm.arr[c1r2][c2];
            tmp.arr[r1][c2] = t;
        }
    return tmp;
}
template<typename T>
Matrix<T> Matrix<T>::operator*(const T t)
{
    Matrix tmp = *this;
    for (int i = 0; i < size.nrow; i++)
        for (int j = 0; j < size.ncol; j++)
            tmp.arr[i][j] *= t;
    return tmp;
}
template<typename T>
Matrix<T> Matrix<T>::operator/(const Matrix<T> rm)  //По сути, деление в данном случае - умножение делимого на обратную матрицу делителя
{
    Matrix invRm = rm.InvertSOLE();
    if (size.ncol != invRm.size.nrow)
        throw runtime_error("\tError! Incompatible sizes!\n");

    T t;
    Matrix tmp(size.nrow, invRm.size.ncol);
    for (int r1 = 0; r1 < size.nrow; r1++)
        for (int c2 = 0; c2 < invRm.size.ncol; c2++)
        {
            t = 0;
            for (int c1r2 = 0; c1r2 < size.ncol; c1r2++)
                t += arr[r1][c1r2] * arr[c1r2][c2];
            tmp.arr[r1][c2] = t;
        }
    return tmp;
}
template<typename T>
Matrix<T> Matrix<T>::operator/(const T t)
{
    Matrix tmp = *this;
    for (int i = 0; i < size.nrow; i++)
        for (int j = 0; j < size.ncol; j++)
            tmp.arr[i][j] /= t;
    return tmp;
}
