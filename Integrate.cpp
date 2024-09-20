#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

// Function f(α, β, z) as defined in the problem
double f(double alpha, double beta, double z) {
    double term = sin(10 * (z - alpha) * (z - beta)) * cos(10 * (z - alpha) * (z - beta));
    return exp(-20 * pow(term, 2));
}

// Simpson's rule for numerical integration
double simpsonsRule(double alpha, double beta, double z_start, double z_end, int n) {
    double h = (z_end - z_start) / n;
    double integral = f(alpha, beta, z_start) + f(alpha, beta, z_end);

    for (int i = 1; i < n; i += 2) {
        integral += 4 * f(alpha, beta, z_start + i * h);
    }
    for (int i = 2; i < n - 1; i += 2) {
        integral += 2 * f(alpha, beta, z_start + i * h);
    }

    return (h / 3) * integral;
}

int main() {
    // Define the range of α and β
    double alpha_start = 0.0, alpha_end = 2 * M_PI;
    double beta_start = 0.0, beta_end = 2 * M_PI;
    double step = 0.01;
    
    int alpha_steps = (alpha_end - alpha_start) / step;
    int beta_steps = (beta_end - beta_start) / step;

    // Open the output file
    std::ofstream outfile("output.dat");

    // Set the precision for output
    outfile << std::fixed << std::setprecision(6);

    // Iterate over α and β
    for (double alpha = alpha_start; alpha <= alpha_end; alpha += step) {
        for (double beta = beta_start; beta <= beta_end; beta += step) {
            // Compute the integral for each pair (α, β)
            double integral_value = simpsonsRule(alpha, beta, 0, M_PI, 1000);  // Using 1000 subdivisions
            
            // Write to the file
            outfile << "α = " << alpha << ", β = " << beta << ", I(α, β) = " << integral_value << std::endl;
        }
    }

    // Close the file
    outfile.close();

    std::cout << "Integration complete. Results saved to output.dat." << std::endl;

    return 0;
}
