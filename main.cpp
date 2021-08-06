#include <iostream>
#include <cmath>

class Equation {
private:
    double epsilon;
    double delta;
    double a;
    double b;
    double c;

public:
    double roots[3]{};
    int number_of_roots = 1;

    Equation(double epsilon, double delta, double a, double b, double c) {
        this->epsilon = epsilon;
        this->delta = delta;
        this->a = a;
        this->b = b;
        this->c = c;
    }

    [[nodiscard]] std::pair<bool, std::pair<double, double>> get_derivative_roots() const {
        double discriminant = 4 * a * a - 12 * b;

        if (discriminant >= epsilon * epsilon) {
            auto x1 = (-2 * a + sqrt(discriminant)) / (6);
            auto x2 = (-2 * a - sqrt(discriminant)) / (6);
            return {true, {x1, x2}};
        }

        return {false, {0.0, 0.0}};
    }

    [[nodiscard]] double function(double x) const {
        return pow(x, 3) + a * pow(x, 2) + b * x + c;
    }

    [[nodiscard]] double dichotomy(std::pair<double, double> interval) const {
        double expected_root;

        if (fabs(function(interval.first)) <= epsilon) {
            return interval.first;
        }

        if (fabs(function(interval.second)) <= epsilon) {
            return interval.second;
        }

        while (true) {
            expected_root = (interval.first + interval.second) / 2;
            double f_mid = function(expected_root);
            if (fabs(f_mid) < epsilon || expected_root == interval.first
                || expected_root == interval.second)
                return expected_root;

            auto f1 = function(interval.first);
            auto f2 = function(interval.second);

            if (f_mid * f1 < 0) {
                interval.second = expected_root;
            }

            if (f_mid * f2 < 0) {
                interval.first = expected_root;
            }
        }
    }

    void roots_calculate() {
        auto[hasRoots, D_roots] = get_derivative_roots();
        std::pair<double, double> interval = {0, 0};

        if (!hasRoots) {
            double f_zero = function(0);
            if (fabs(f_zero) < epsilon) {
                roots[0] = 0;
                number_of_roots = 1;

                return;
            }

            if (f_zero < -epsilon) {
                int i = 0;
                while (true) {
                    auto f_alpha_plus_delta_1 = function(delta * i);
                    auto f_alpha_plus_delta_2 = function(delta * (i + 1));
                    i++;
                    if ((f_alpha_plus_delta_1 * f_alpha_plus_delta_2) <= 0) break;
                }
                interval.first = function(delta * i);
                interval.second = function(delta * (i + 1));

                roots[0] = dichotomy(interval);
                number_of_roots = 1;

                return;
            }
            if (f_zero > epsilon) {
                int i = 0;
                while (true) {
                    auto f_alpha_plus_delta_1 = function(delta * i);
                    auto f_alpha_plus_delta_2 = function(delta * (i - 1));
                    i--;
                    if ((f_alpha_plus_delta_1 * f_alpha_plus_delta_2) <= 0) break;
                }
                interval.first = delta * i;
                interval.second = delta * (i - 1);

                roots[0] = dichotomy(interval);
                number_of_roots = 1;

                return;
            }
        }

        auto[x1, x2] = D_roots;

        double alpha = x1 < x2 ? x1 : x2;
        double beta = x1 > x2 ? x1 : x2;

        if (function(alpha) > epsilon && function(beta) > epsilon) {
            int i = 0;
            while (true) {
                auto f_alpha_plus_delta_1 = function(alpha + delta * i);
                auto f_alpha_plus_delta_2 = function(alpha + delta * (i - 1));
                if ((f_alpha_plus_delta_1 * f_alpha_plus_delta_2) < 0) break;
                i--;
            }
            interval.first = alpha + delta * i;
            interval.second = alpha + delta * (i - 1);

            roots[0] = dichotomy(interval);
            number_of_roots = 1;

            return;
        }
        if (function(alpha) < -epsilon && function(beta) < -epsilon) {
            int i = 0;
            while (true) {
                auto f_beta_plus_delta_1 = function(beta + delta * i);
                auto f_beta_plus_delta_2 = function(beta + delta * (i + 1));
                if ((f_beta_plus_delta_1 * f_beta_plus_delta_2) < 0) break;
                i++;
            }
            interval.first = beta + delta * i;
            interval.second = beta + delta * (i + 1);

            roots[0] = dichotomy(interval);
            number_of_roots = 1;

            return;
        }
        if (function(alpha) > epsilon && fabs(function(beta)) < epsilon) {
            int i = 0;
            while (true) {
                auto f_alpha_plus_delta_1 = function(alpha + delta * i);
                auto f_alpha_plus_delta_2 = function(alpha + delta * (i - 1));
                if ((f_alpha_plus_delta_1 * f_alpha_plus_delta_2) < 0) break;
                i--;
            }
            interval.first = alpha + delta * i;
            interval.second = alpha + delta * (i - 1);

            roots[0] = beta;
            roots[1] = dichotomy(interval);
            number_of_roots = 2;

            return;
        }
        if (fabs(function(alpha)) < epsilon && function(beta) < -epsilon) {
            int i = 0;
            while (true) {
                auto f_beta_plus_delta_1 = function(beta + delta * i);
                auto f_beta_plus_delta_2 = function(beta + delta * (i + 1));
                if ((f_beta_plus_delta_1 * f_beta_plus_delta_2) < 0) break;
                i++;
            }

            interval.first = beta + delta * i;
            interval.second = beta + delta * (i + 1);

            roots[0] = alpha;
            roots[1] = dichotomy(interval);
            number_of_roots = 2;

            return;
        }
        if (function(alpha) > epsilon && function(beta) < -epsilon) {
            int i = 0;
            while (true) {
                auto f_alpha_plus_delta_1 = function(alpha + delta * i);
                auto f_alpha_plus_delta_2 = function(alpha + delta * (i - 1));
                if ((f_alpha_plus_delta_1 * f_alpha_plus_delta_2) < 0) break;
                i--;
            }
            interval.first = alpha + delta * i;
            interval.second = alpha + delta * (i - 1);

            roots[0] = dichotomy(interval);

            interval.first = alpha;
            interval.second = beta;

            roots[1] = dichotomy(interval);

            i = 0;
            while (true) {
                auto f_beta_plus_delta_1 = function(beta + delta * i);
                auto f_beta_plus_delta_2 = function(beta + delta * (i + 1));
                if ((f_beta_plus_delta_1 * f_beta_plus_delta_2) < 0) break;
                i++;
            }
            interval.first = beta * delta * i;
            interval.second = beta + delta * (i + 1);

            roots[2] = dichotomy(interval);
            number_of_roots = 3;

            return;
        }
        if (fabs(function(alpha)) < epsilon && fabs(function(beta)) < epsilon) {
            roots[0] = (alpha + beta) / 2;
            number_of_roots = 1;

            return;
        }
    }
};

int main(int argc, char *argv[]) {
    double epsilon = 10E-10;
    double delta = 0.1;

    if (4 != argc) {
        std::cout << "Bad arguments, expected equation coefficients:\n\t x^3 + a*x^2 + b*x + c\n";
    }

    double a = strtod(argv[1], nullptr);
    double b = strtod(argv[2], nullptr);
    double c = strtod(argv[3], nullptr);

    Equation eq(epsilon, delta, a, b, c);
    eq.roots_calculate();

    for (int i = 0; i < eq.number_of_roots; i++) {
        std::cout << "\t[" << i+1 << "] " << eq.roots[i] << std::endl;
    }

    return EXIT_SUCCESS;
}