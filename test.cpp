#include <cmath> 
#include <vector>
#include <iostream>

using namespace std;


double F(double x1, double x2)
{
    return pow((x1 + 6 * x2), 2) + pow((x1 + 2), 2);  //заданная функция
}
double Pr1(double x1, double x2)      //1-я производная по х1
{
    return 4 * x1 + 12 * x2 + 4;
}
double Pr2(double x1, double x2)        //1-я производная по х2
{
    return 12 * x1 + 72 * x2;
}




void gradient_metod() {
    double x1_prev, x2_prev;
    double x1 = -2, x2 = 0.4, eps = 0.01, alpha = 0.1;//alpha - шаг
    double dx1;
    double dx2;
    double Y1, Y2;
    int max_iter = 100;

    int k = 0, l = 0;

    do
    {
        x1_prev = x1;
        x2_prev = x2;

        dx1 = Pr1(x1, x2);//вычисляем 1-ю производную по х1
        dx2 = Pr2(x1, x2);//вычисляем 1-ю производную по х2    

        x1 = x1 - alpha * dx1;
        x2 = x2 - alpha * dx2;

        Y1 = F(x1_prev, x2_prev);
        Y2 = F(x1, x2);

        k++;
        alpha = alpha / 2;
        std::cout << "\nIteration " << k << "   x1 = " << x1 << "    x2 = " << x2 << "  dx1 = " << dx1 << "  dx2 = " << dx2 << "   F(x1,x2) = " << F(x1, x2);


        if (k > max_iter) break;

    } while ((abs(x1 - x1_prev) >= eps) && abs(x2 - x2_prev) > eps);


    std::cout << "\nSolution: x* = " << x1 << "    x2* = " << x2 << "   F(x1*,x2*) = " << F(x1, x2);
}


double func(double x1, double x2) {
    return pow((x1 + 6 * x2), 2) + pow((x1 + 2), 2);
}

// Вычисление градиента функции
vector<double> computeGradient(double x1, double x2) {
    vector<double> gradient;
    double df_dx1 = 2 * (x1 + 6 * x2) + 2 * (x1 + 2);
    double df_dx2 = 12 * (x1 + 6 * x2);
    gradient.push_back(df_dx1);
    gradient.push_back(df_dx2);
    return gradient;
}

// Вычисление матрицы Гессе функции
vector<vector<double>> computeHessian(double x1, double x2) {
    vector<vector<double>> hessian;
    vector<double> row1, row2;
    double d2f_dx1_dx1 = 2 + 2;
    double d2f_dx1_dx2 = 2 * 6;
    double d2f_dx2_dx1 = 2 * 6;
    double d2f_dx2_dx2 = 12 * 12;
    row1.push_back(d2f_dx1_dx1);
    row1.push_back(d2f_dx1_dx2);
    row2.push_back(d2f_dx2_dx1);
    row2.push_back(d2f_dx2_dx2);
    hessian.push_back(row1);
    hessian.push_back(row2);
    return hessian;
}

void marquardtMethod(double initialX1, double initialX2) {
    double x1 = initialX1;
    double x2 = initialX2;
    double lambda = 0.1;    // Инициализация параметра lambda
    double tolerance = 1e-3; // Допустимая погрешность
    int maxIterations = 100; // Максимальное количество итераций

    int sch = 0;

    for (int i = 0; i < maxIterations; i++) {
        vector<double> gradient = computeGradient(x1, x2);
        vector<vector<double>> hessian = computeHessian(x1, x2);

        vector<vector<double>> modifiedHessian = hessian;
        modifiedHessian[0][0] += lambda;
        modifiedHessian[1][1] += lambda;

        double det = modifiedHessian[0][0] * modifiedHessian[1][1] - modifiedHessian[0][1] * modifiedHessian[1][0];

        // Проверка на невырожденность матрицы
        if (det != 0) {
            double deltaX1 = det * (modifiedHessian[1][1] * gradient[0] - modifiedHessian[0][1] * gradient[1]) / (det * det + lambda);
            double deltaX2 = det * (modifiedHessian[0][0] * gradient[1] - modifiedHessian[1][0] * gradient[0]) / (det * det + lambda);

            // Обновление значений переменных x1 и x2
            x1 -= deltaX1;
            x2 -= deltaX2;

            // Проверка на достижение допустимой погрешности
            if (abs(deltaX1) < tolerance && abs(deltaX2) < tolerance) {
                break;
            }
        }

        if (lambda < 10) {
            lambda *= 10;
        }
        sch = i;
    }

    std::cout << "Function points: (" << x1 << ", " << x2 << ")" << endl;
    std::cout << "Function value: " << func(x1, x2) << endl;
    std::cout << "Iteration number: " << sch << endl;

}

double f(double x1, double x2) {
    return pow((x1 + 6 * x2), 2) + pow((x1 + 2), 2);
}

double dfdx1(double x1, double x2) {
    return 2 * (x1 + 6 * x2) + 2 * (x1 + 2);
}

double dfdx2(double x1, double x2) {
    return 12 * (x1 + 6 * x2);
}

void gradient_descent(double& x1, double& x2, double learning_rate, int iterations) {
    for (int i = 0; i < iterations; ++i) {
        double gradient_x1 = dfdx1(x1, x2);
        double gradient_x2 = dfdx2(x1, x2);

        x1 -= learning_rate * gradient_x1;
        x2 -= learning_rate * gradient_x2;
    }
}


void gradient() {

    double x1 = 0.0; // starting point for x1
    double x2 = 0.0; // starting point for x2
    double learning_rate = 0.01; // learning rate
    int iterations = 100; // number of iterations

    gradient_descent(x1, x2, learning_rate, iterations);

    std::cout << "Point of extremum: (" << x1 << ", " << x2 << ")" << std::endl;
    std::cout << "Minimum value of f(x1, x2): " << f(x1, x2) << std::endl;

}

int main() {


    std::cout << "---------------Gradient metod with fractional step---------------" << endl;
    gradient_metod();
    std::cout << endl;

    std::cout << "---------------Marquard metod---------------" << endl;
    double initialX1 = 0;
    double initialX2 = 0;
    marquardtMethod(initialX1, initialX2);
    std::cout << endl;

    std::cout << "---------------Fastest gradient metod---------------" << endl;
    gradient();
    return 0;
}
