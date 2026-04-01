#define _USE_MATH_DEFINES // Garante que o M_PI seja habilitado em todos os compiladores
#include <stdio.h>
#include <math.h>
#include "pbPlots.h"
#include "supportLib.h"

typedef struct {
    double B;      // Bottom width (m)
    double m;      // Side slope (H:V)
    double n;      // Manning's roughness
    double S0;     // Slope (m/m)
    double Q;      // Flow rate (m3/s)
    double Y;      // Flow depth (m)
    double H_max;  // Total wall height (m)
} Section;

// Funções auxiliares
double calcularArea(double y, Section s) {
    return (s.B + s.m * y) * y;
}

double calcularPerimetro(double y, Section s) {
    return s.B + 2.0 * y * sqrt(1.0 + s.m * s.m);
}

double calcularLarguraTopo(double y, Section s) {
    return s.B + 2.0 * s.m * y;
}



void ManningSolver(int Case, Section *Sect) {

    // Manning: Q = (1/n) * (A^5/3)/(P^2/3) * (S0^1/2) 

    /*
    Case Known      Quantities      Unknown Quantities
        1           n,b,Y,m,S               Q
        2           Q,b,Y,m,S               n
        3           Q,n,Y,m,S               b
        4           Q,n,b,m,S               Y
        5           Q,n,b,Y,S               m
        6           Q,n,b,Y,m               S
    */

    /*
        F(r) = n * Q * P^(2/3) - A^(5/3) * S0^(1/2)= 0

        dF/dr = 2/3 * n * Q * P^(-1/3) * dP/dr - 5/3 * S0^(1/2) * A^(2/3) * dA/dr
    */

    /*
        Y unknown: DA/DY = b +2mY = T       DP/DY = 2(m^2 + 1)^1/2
        m unknown: DA/Dm = Y^2               DP/Dm = 2mY/(m^2 + 1)^1/2
        b unknown: DA/Db = Y                DP/Db = 1
    */


    if (Case == 1 || Case == 2 || Case == 6) {
        double A = calcularArea(Sect->Y, *Sect);
        double P = calcularPerimetro(Sect->Y, *Sect);
        if (Case == 1) Sect->Q = (1.0 / Sect->n) * pow(A, 5.0/3.0) * pow(P, -2.0/3.0) * sqrt(Sect->S0);
        if (Case == 2) Sect->n = (1.0 / Sect->Q) * pow(A, 5.0/3.0) * pow(P, -2.0/3.0) * sqrt(Sect->S0);
        if (Case == 6) Sect->S0 = pow((Sect->Q * Sect->n * pow(P, 2.0/3.0)) / pow(A, 5.0/3.0), 2.0);
        return;
    }

    double y_inic = Sect->H_max / 2.0; // Chute inicial (metade da altura)
    double b_calc = Sect->B;
    double m_calc = Sect->m;
    double erro = 1e-5;
    int max_iter = 100, iter = 0;
    double F, dF, dA, dP, A, P;
    double y_calc = (Case == 4) ? y_inic : Sect->Y;

    do {
        Section temp = *Sect;
        temp.B = b_calc;
        temp.m = m_calc;
        
        A = calcularArea(y_calc, temp);
        P = calcularPerimetro(y_calc, temp);
        

        if (Case == 3)
        {
            dP = 1.0;
            dA = y_calc;
        }

        if (Case == 4)
        {
            dP = 2.0 * pow((pow(m_calc, 2.0) + 1.0), 0.5);
            dA = b_calc + 2.0 * m_calc * y_calc;
        }

        if (Case == 5)
        {
            dP = (2.0 * m_calc * y_calc)/ pow((pow(m_calc, 2.0) + 1.0), 0.5);
            dA = pow(y_calc, 2.0);
        }
        
        
        F = Sect->n * Sect->Q * pow(P, 2.0/3.0) - pow(A, 5.0/3.0) * sqrt(Sect->S0);
        
        dF = ((2.0/3.0) * Sect->n * Sect->Q * pow(P, -1.0/3.0) * dP) - 
            ((5.0/3.0) * sqrt(Sect->S0) * pow(A, 2.0/3.0) * dA);

        double delta = F / dF;

        
        if (Case == 3) b_calc -= delta;
        else if (Case == 4) y_calc -= delta;
        else if (Case == 5) m_calc -= delta;

        iter++;

    } while (fabs(F) > erro && iter < max_iter);

    Sect->B = b_calc;
    Sect->Y = y_calc;
    Sect->m = m_calc;

    if (Case == 3) printf("Largura da base (b) calculada: %.4f m\n", b_calc);
    if (Case == 4) printf("Profundidade (Y) calculada: %.4f m\n", y_calc);
    if (Case == 5) printf("Inclinação (m) calculada: %.4f\n", m_calc);
    printf("Number of Iterations: %d", iter);

}


void EfficientSection(Section *Sect) {
    Sect->m = sqrt(3.0) / 3.0; // m = 0.57735
    
    double constante_SI = 0.70205;
    Sect->B = pow( (Sect->n * Sect->Q) / (constante_SI * sqrt(Sect->S0)), 3.0/8.0 );
    
    Sect->Y = Sect->B * (sqrt(3.0) / 2.0); 

    Sect->H_max = Sect->Y * 1.2;
}


void PlotChannel(Section *Sect) {
    double DisBase = Sect->B/2;
    double DisTopo = Sect->m * Sect->H_max;

    double xs[] = {
        -(DisBase + DisTopo),
        -(DisBase),
        (DisBase),
        (DisBase + DisTopo)
    };

    double ys[] = {
        Sect->H_max,
        0.0,
        0.0,
        Sect->H_max,
    };

    double xs_water[] = {
        -(DisBase + (Sect->m * Sect->Y)),
        (DisBase + (Sect->m * Sect->Y))
    };

    double ys_water[] = {
        Sect->Y,
        Sect->Y,
    };


    _Bool success;
    StartArenaAllocator();


    ScatterPlotSeries *seriesCanal = GetDefaultScatterPlotSeriesSettings();
    seriesCanal->xs = xs;
    seriesCanal->xsLength = 4;
    seriesCanal->ys = ys;
    seriesCanal->ysLength = 4;
    seriesCanal->linearInterpolation = true;
    seriesCanal->lineType = L"solid";
    seriesCanal->lineThickness = 2;
    seriesCanal->color = CreateRGBColor(0, 0, 0); // Preto

    
    ScatterPlotSeries *seriesAgua = GetDefaultScatterPlotSeriesSettings();
    seriesAgua->xs = xs_water;
    seriesAgua->xsLength = 2;
    seriesAgua->ys = ys_water;
    seriesAgua->ysLength = 2;
    seriesAgua->linearInterpolation = true;
    seriesAgua->lineType = L"solid";
    seriesAgua->lineThickness = 3;
    seriesAgua->color = CreateRGBColor(0, 0, 1); // Azul

    
    ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
    settings->width = 800;
    settings->height = 600;
    settings->autoBoundaries = true;
    settings->autoPadding = true;
    settings->title = L"Channel profile";
    settings->titleLength = wcslen(settings->title);
    settings->xLabel = L"X(m)";
    settings->xLabelLength = wcslen(settings->xLabel);
    settings->yLabel = L"Y(m)";
    settings->yLabelLength = wcslen(settings->yLabel);


    ScatterPlotSeries *seriesList[] = {seriesCanal, seriesAgua};
    settings->scatterPlotSeries = seriesList;
    settings->scatterPlotSeriesLength = 2;

    RGBABitmapImageReference *canvasReference = CreateRGBABitmapImageReference();
    StringReference *errorMessage = CreateStringReference(L"", 0);

    success = DrawScatterPlotFromSettings(canvasReference, settings, errorMessage);

    if (success) {
        ByteArray *pngdata = ConvertToPNG(canvasReference->image);
        WriteToFile(pngdata, "canal_trapezoidal.png");
        printf("\n[Plot] Grafico gerado com sucesso: canal_trapezoidal.png\n");
        DeleteImage(canvasReference->image);
    }
    else {
        fprintf(stderr, "\n[Plot] Erro ao gerar grafico.\n");
    }

    FreeAllocations();

}

int main() {
    Section mySection = {0};
    int solvingCase;
    double slopePercent;

    printf("--- Manning's Equation Solver ---\n");
    printf("Select the unknown variable:\n");
    printf("1 - Flow Rate (Q)\n2 - Manning's n\n3 - Bottom Width (b)\n");
    printf("4 - Flow Depth (Y)\n5 - Side Slope (m)\n6 - Slope (S0)\n");
    printf("Choice: ");
    scanf("%d", &solvingCase);

    // Input flow parameters (always needed except for Q or n)
    if (solvingCase != 1) { printf("Flow Rate Q (m3/s): "); scanf("%lf", &mySection.Q); }
    if (solvingCase != 2) { printf("Manning's n: "); scanf("%lf", &mySection.n); }
    if (solvingCase != 6) { 
        printf("Slope in %% (e.g., 0.5 for 0.5%%): "); 
        scanf("%lf", &slopePercent); 
        mySection.S0 = slopePercent / 100.0;
    }

    // Input geometry parameters
    if (solvingCase != 3) { printf("Bottom Width B (m): "); scanf("%lf", &mySection.B); }
    if (solvingCase != 5) { printf("Side Slope m (H:V): "); scanf("%lf", &mySection.m); }
    if (solvingCase != 4) { printf("Flow Depth Y (m): "); scanf("%lf", &mySection.Y); }
    
    // H_max is always useful for the Y-initial guess
    printf("Total Wall Height H_max (m): ");
    scanf("%lf", &mySection.H_max);

    printf("\n--- Solving ---\n");
    ManningSolver(solvingCase, &mySection);

    printf("\n--- Results ---\n");
    printf("Q: %.3f m3/s | n: %.4f | S0: %.6f\n", mySection.Q, mySection.n, mySection.S0);
    printf("B: %.3f m | m: %.3f | Y: %.4f m\n", mySection.B, mySection.m, mySection.Y);

    PlotChannel(&mySection);


    EfficientSection(&mySection);

    printf("\n--- Most Efficient Results ---\n");
    printf("Calculated B: %.4f m\n", mySection.B);
    printf("Calculated Y: %.4f m\n", mySection.Y);
    printf("Fixed m: %.4f\n", mySection.m);

    system("PAUSE");
    return 0;
}