#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
void print_menu()
{
    printf("\nWhat do you want to watch?\n");
    printf("0. exit\n");
    printf("1. 9.5\n");
    printf("2. 9.6\n");

}
void preobr(int p, int n, double **mas, int *ans)
{
    int i, j;
    int imax = p, jmax = p;

    for (i = p; i < n; i++)
        for (j = p; j < n; j++)
            if (fabs(mas[imax][jmax]) < fabs(mas[i][j]))
              {
                  imax = i;

                  jmax = j;

              }

    for (j = p; j < n + 1; j++)
     {
        double t = mas[p][j];
        mas[p][j] = mas[imax][j];   //Меняем местами строки
        mas[imax][j] = t;
     }

    for (i = 0; i < n; i++)
     {
        double t = mas[i][p];
        mas[i][p] = mas[i][jmax];   //Меняем местами столбцы
        mas[i][jmax] = t;
     }
    i = ans[p];

    ans[p] = ans[jmax];

    ans[jmax] = i;
}

void program95()
{
    setlocale(LC_ALL, "Rus");

    double **mas;
    double *x;
    int i,j,p;
    int *ans;
    int n;

    printf("Введите количество уравнений в системе: ");
    scanf("%d", &n);

    mas = (double**) malloc(n*sizeof(double*));
    x = (double*) malloc(n*sizeof(double));
    ans = (int*) malloc(n*sizeof(int));

    for (int i = 0; i < n; ++i)
     {
        mas[i] = (double*) malloc((n+1)*sizeof(double));
     }

    printf("Введите строки матрицы (коэффициенты при переменных в каждом уравнении и также свободный член в нём)\n");
    for (i = 0; i < n; ++i)
        for (j = 0; j <= n; ++j)
            scanf("%lf", &mas[i][j]);



    for (i = 0; i < n; ++i)
        ans[i] = i;

    for (p = 0; p < n; ++p)
    {
        preobr(p, n, mas, ans);

        if (fabs(mas[p][p]) < 0.00001 )
         {
            printf("");
            return 0;
         }

        for (j = n; j >= p; --j)
            mas[p][j]= mas[p][j]/mas[p][p];


        for (i = p + 1; i < n; ++i)
            for (j = n; j >= p; --j)
                mas[i][j] = mas[i][j] - (mas[p][j] * mas[i][p]);

    }

    for ( i = 0; i < n; ++i)
        x[i] = mas[i][n];

    for (i = n - 2; i >= 0; --i)
        for (j = i + 1; j < n; ++j)
            x[i]= x[i] - (x[j] * mas[i][j]);

    printf("\n");

    printf("Ответ:\n");
    for (i = 0; i < n; ++i)
        for (j = 0; j <= n; ++j)

            if (i == ans[j])
            {
                printf("%f\n", x[j]);
                break;
            }

}

void gauss(double **mas, double **masr, int n)
{
    int k = 0, ind;
    while (k < n)
    {
        double t = mas[k][k];
        if(t == 0)
        {
            for(int i = 0; i < n; i++)
                if(mas[i][k] != 0)
                {
                    ind = i;
                    break;
                }
            for(int i = 0; i < n; i++)
            {
                double t = mas[k][i], t1 = masr[k][i];
                mas[k][i] = mas[ind][i], masr[k][i] = masr[ind][i];
                mas[ind][i] = t, masr[ind][i] = t1;
            }
            t = mas[k][k];
        }
        for (int i = 0; i < n; i++)
        {
            mas[k][i] /= t;
            masr[k][i] /= t;
        }
        for (int i = k + 1; i < n; i++)
        {
            double t = mas[i][k];
            for (int j = 0; j < n; j++)
            {
                mas[i][j] -= mas[k][j] * t;
                masr[i][j] -= masr[k][j] * t;
            }
        }
        k++;
    }
    for (k = n - 1; k > 0; k--)
    {
        for(int j = k - 1; j >= 0; j--)
        {
            double t = mas[j][k];
            mas[j][k] -= mas[k][k] * t;
            for(int i = n-1; i >= 0; i--)
            {
                masr[j][i] -= masr[k][i] * t;
            }
        }
    }
}
void program96()
{
    setlocale(LC_ALL, "Rus");
    int n;
    double **mas, **masr, **copy_m;

    printf("Введите размерность матрицы m*m: ");
    scanf("%d", &n);

    mas = (double**) malloc(n*sizeof(double*));
    masr = (double**) malloc(n*sizeof(double*));
    copy_m = (double**) malloc(n*sizeof(double*));


    for (int i = 0; i < n; ++i){
        mas[i] = (double*) malloc(n*sizeof(double));
        masr[i] = (double*) malloc(n*sizeof(double));
        copy_m[i] = (double*) malloc(n*sizeof(double));
    }
    printf("Введите саму матрицу: \n");

    for (int i = 0; i < n; i++)

        for (int j = 0; j < n; j++)

        {
            if(i == j)
                masr[i][j] = 1;
                  else masr[i][j] = 0;

            scanf("%lf", &mas[i][j]);
            copy_m[i][j] = mas[i][j];
        }

    gauss(mas, masr, n);
    printf("\nОбратная матрица: \n");

     for(int i = 0; i < n; i++)
      {
        for(int j = 0; j < n; j++)
            printf("%lf ", masr[i][j]);
        printf("\n");
      }

    double **c = (double**) malloc(n*sizeof(double*));

    for (int i = 0; i < n; i++)
    {
        c[i] = (double*) malloc(n*sizeof(double));

        for (int j = 0; j < n; j++)

        {
            c[i][j] = 0;

            for (int k = 0; k < n; k++)
                c[i][j] += copy_m[i][k] * masr[k][j];

            if (fabs(c[i][j]) < 0.000001)
                c[i][j] = 0;
        }

    }
    printf("\nРезультат перемножения матриц: \n");
    for(int i = 0; i < n; i++)
    {
        for(int j=0; j<n; j++)
        {
            printf("%lf ", c[i][j]);
        }
        printf("\n");
    }

}

int main()
{
    setlocale(LC_ALL,"Rus");
    int k=1, variant;

    while (k)
    {
        print_menu();
        printf("\nchoose the number = ");
        scanf("%d", &variant);

        switch (variant)
        {
        case 0:
            k=0;
            break;
        case 1:
            program95();
            break;

        case 2:
            program96();
            break;

        }

        variant=0;
    }


    return 0;
}


