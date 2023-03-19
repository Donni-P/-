#include <iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include <string>

using namespace std;

const double pi = 4 * atan(1);
const short N = 32;
void SIGNAL(double*);
void HINDRANCE(double*);
void DPF(double d[][2], double*);
void ODPF(double d[][2], double*);

int main()
{
    setlocale(LC_ALL, "Russian");
    ofstream fileOut;
    fileOut.open("out.txt");
    double dpf[N][2], odpf[N], dpfh[N][2];
    double signal[N], hindrance[N], hs[N];
    SIGNAL(signal);
    DPF(dpf, signal);
    fileOut << "---------------------------------------------------------------------------------------------------------\n";
    fileOut << "u |                 ДПФ(Re + j*Im)                      |      ДПФ Амплитуда    |        ДПФ Фаза       |\n";
    fileOut << "---------------------------------------------------------------------------------------------------------\n";
    for (short u = 0; u < N; u++) {
        fileOut << setw(2) << u << "|";
        fileOut << fixed << setw(23) << setprecision(20) << dpf[u][0];
        fileOut << "  +  j*" << setw(23) << dpf[u][1] << "|" << setw(23) << sqrt(dpf[u][0] * dpf[u][0] + dpf[u][1] * dpf[u][1]);
        fileOut << "|" << setw(23) << atan2(dpf[u][1], dpf[u][0]) << "|\n";
    }
    ODPF(dpf, odpf);
    fileOut << "---------------------------------------------------------------------------------------------------------\n";
    fileOut << "k |     Входной сигнал    |         ОДПФ          | Ошибка восстановления |\n";
    fileOut << "---------------------------------------------------------------------------\n";
    for (short k = 0; k < N; k++) {
        fileOut << setw(2) << k << "|";
        fileOut << setw(23) << signal[k] << "|" << setw(23) << odpf[k];
        fileOut << "|" << setw(23) << fabs(odpf[k] - signal[k]) << "|\n";
    }
    HINDRANCE(hindrance);
    DPF(dpfh, hindrance);
    for (short k = 0; k < N; k++)
        hs[k] = signal[k] + hindrance[k];
    DPF(dpf, hs);
    for (short u = 0; u < N; u++)
    {
        dpf[u][0] -= dpfh[u][0];
        dpf[u][1] -= dpfh[u][1];
    }
    ODPF(dpf, odpf);
    fileOut << "---------------------------------------------------------------------------------------------------------------------------\n";
    fileOut << "k |         Cигнал        |         Помеха        |      Сигнал+помеха    |     Фильтрация ОДПФ   |         Ошибка        |\n";
    fileOut << "---------------------------------------------------------------------------------------------------------------------------\n";
    for (short k = 0; k < N; k++) {
        fileOut << setw(2) << k << "|" << setw(23) << signal[k] << "|" << setw(23) << hindrance[k] << "|";
        fileOut << setw(23) << hs[k] << "|" << setw(23) << odpf[k] << "|";
        fileOut << setw(23) << fabs(odpf[k] - signal[k]) << "|" << endl;
    }
    fileOut << "---------------------------------------------------------------------------------------------------------------";
    fileOut << "------------------------------------------------------\n";
    fileOut << "u |              ДПФ помехи(Re + j*Im)                  |          ДПФ сигнал+помеха(Re + j*Im)               |";
    fileOut << "         отфильтрованный ДПФ(Re + j*Im)              |\n";
    fileOut << "---------------------------------------------------------------------------------------------------------------";
    fileOut << "------------------------------------------------------\n";
    for (short u = 0; u < N; u++) {
        fileOut << setw(2) << u << "|" << setw(23) << dpfh[u][0] << "  +  j*" << setw(23) << dpfh[u][1] << "|";
        fileOut << setw(23) << dpf[u][0] + dpfh[u][0] << "  +  j*" << setw(23) << dpf[u][1] + dpfh[u][1] << "|";
        fileOut << setw(23) << dpf[u][0] << "  +  j*" << setw(23) << dpf[u][1] << "|\n";
    }
    fileOut << "---------------------------------------------------------------------------------------------------------------";
    fileOut << "------------------------------------------------------\n";
    fileOut.close();
    ifstream file("out.txt");
    string line;
    while (getline(file, line))
        cout << line << endl;
    file.close();
}
void HINDRANCE(double* h) {
    for (short k = 0; k < N; k++)
        h[k] = cos(20.0 * pi * k / N - pi / 4) + sin(22.0 * pi * k / N + pi / 4);
}
double sign(double x) { return (x > 0) ? 1.0 : (x < 0) ? -1.0 : 0.0; }
void SIGNAL(double* s) {
    for (short k = 0; k < N; k++)
        s[k] = sign(cos(2.0 * pi * k / N));
}
void DPF(double d[][2], double* sig)
{
    double sumRe, sumIm;
    for (short u = 0; u < N; u++)
    {
        sumRe = 0;
        sumIm = 0;
        for (short k = 0; k < N; k++)
        {
            sumRe += sig[k] * cos(2.0 * pi * k * u / N);
            sumIm += sig[k] * sin(2.0 * pi * k * u / N);
        }
        d[u][0] = 1.0 / sqrt(N) * sumRe;
        d[u][1] = 1.0 / sqrt(N) * sumIm;
    }
}
void ODPF(double d[][2], double* s)
{
    double sum;
    for (short k = 0; k < N; k++)
    {
        sum = 0;
        s[k] = 0.0;
        for (short u = 0; u < N; u++)
        {
            sum += (d[u][0] * cos(2.0 * pi * k * u / N) + d[u][1] * sin(2.0 * pi * k * u / N));
        }
        s[k] = 1.0 / sqrt(N) * sum;
    }
}
