#include <iostream>
#include <cmath>
#include <vector>
#include <sciplot/sciplot.hpp>
#include <random>

/**
 * Generates a vector of uniformly distributed random values in the range [-1, 1].
 *
 * @param n the number of random values to generate
 *
 * @return a vector of random values
 *
 * @throws None
 */
std::vector<double> random_values(int n) {
    //initialize random number generator
    std::random_device rd;
    std::default_random_engine engine(rd());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    std::vector<double> values(n);
    for (int i = 0; i < n; i++){
        values[i] = dist(engine);
    }

    return values;
}

/**
 * Generates a 2D vector of random values for the x and y coordinates of points.
 *
 * @param n the number of random values to generate for each coordinate
 *
 * @return a 2D vector of size 2, where the first element is a vector of random x values and the second element is a vector of random y values
 *
 * @throws None
 */
std::vector<std::vector<double>> points(int n) {
    std::vector<std::vector<double>> points(2, std::vector<double>(n));
    points[0] = random_values(n);
    points[1] = random_values(n);
    return points;
}

/**
 * Calculates the Euclidean norm of a 2D point.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 *
 * @return The Euclidean norm of the point.
 *
 * @throws None.
 */
double euclidian_norm(double x, double y) {
    return sqrt(pow(x, 2) + pow(y, 2));
}

/**
 * Determines whether the given point is inside the unit circle.
 *
 * @param x the x-coordinate of the point to check
 * @param y the y-coordinate of the point to check
 *
 * @return true if the point is inside the circle, false otherwise
 *
 * @throws None
 */
bool in_circle(double x, double y) {
    return euclidian_norm(x, y) <= 1;
}



/**
 * Generates a 2D vector of points on the unit circle.
 *
 * @return A 2D vector of size 2, where the first element is a vector of x-coordinates and the second element is a vector of y-coordinates.
 *
 * @throws None
 */
std::vector<std::vector<double>> unit_circle() {
    //define pi
    const double pi = 4 * atan(1);
    int n = 1000;
    std::vector<std::vector<double>> points(2, std::vector<double>(n));
    for (int i = 0; i < n; i++) {
        points[0][i] = cos(2 * pi * i / n);
        points[1][i] = sin(2 * pi * i / n);
    }
    return points;
}

/**
 * Reads a vector of points and returns a vector containing the partitioned points inside and outside a circle.
 *
 * @param points The input vector of points.
 *
 * @return A 2D vector containing two vectors: one for points inside the circle and one for points outside the circle.
 *
 * @throws None
 */
std::vector<std::vector<std::vector<double>>> inside_outside(std::vector<std::vector<double>> points) {
    int n = points[0].size();
    std::vector<std::vector<double>> inside(2);
    std::vector<std::vector<double>> outside(2);
    for (int i = 0; i < n; i++) {
        if (in_circle(points[0][i], points[1][i])) {
            inside[0].push_back(points[0][i]);
            inside[1].push_back(points[1][i]);
        } else {
            outside[0].push_back(points[0][i]);
            outside[1].push_back(points[1][i]);
        }
    }
    return {inside, outside};
}

/**
 * Creates a 2D scatter plot of generated points.
 *
 * @param n The number of points to generate.
 *
 * @return void
 *
 * @throws None
 */
void plot_points(int n){
    //initializes vector of n random points
    std::vector<std::vector<double>> point_vector = points(n);

    //partitions the points into those inside and those outside the unit circle 
    std::vector<std::vector<std::vector<double>>> in_out = inside_outside(point_vector);
    std::vector<std::vector<double>> points_inside = in_out[0];
    std::vector<std::vector<double>> points_outside = in_out[1];

    int inside = points_inside[0].size();
    int outside = points_outside[0].size();

    //calculates approximation of pi
    double ratio = double(inside) / double(inside + outside);
    double pi = ratio * 4;
    
    //creates points on the unit circle
    std::vector<std::vector<double>> unit = unit_circle();

    //plotting the unit circle and the points that are inside and outside of it
    sciplot::Plot2D plot;
    plot.drawPoints(points_inside[0], points_inside[1]).pointType(7).pointSize(0.5).lineColor("green");
    plot.drawPoints(points_outside[0], points_outside[1]).pointType(7).pointSize(0.5).lineColor("blue");
    plot.drawCurve(unit[0], unit[1]);
    plot.xlabel("x");
    plot.ylabel("y");
    plot.legend().hide();
    sciplot::Figure figure = {{plot}};
    sciplot::Canvas canvas = {{figure}};
    canvas.show();
}

/**
 * Plots the relative error of the Monte Carlo estimation of the value of pi.
 *
 * @param power The power of 2 used to calculate the maximum number of points in the Monte Carlo estimation.
 *
 * @return void
 *
 * @throws None
 */
void plot_error(int power){
    const double pi_actual = 4 * atan(1);
    //skips the case for n = 1 to reduce error a bit
    int n_to_skip = 1;
    //runs simulation num_runs times and averages the results to minimize stochastic effects
    int num_runs = 100;
    int length = power - n_to_skip + 1;

    //initializes vectors
    std::vector<double> relative_error(length);
    std::vector<double> logn(length);

    for (int run = 0; run < num_runs; run++){
        for (int i = 0; i < length; i++){
            //calculates the value of pi using n points (n = 2^(i+n_to_skip-1))
            int idx = i + n_to_skip - 1;
            double n = pow(2, idx);
            logn[i] = log2(n);
            std::vector<std::vector<double>> point_vector = points(n);
            std::vector<std::vector<std::vector<double>>> in_out = inside_outside(point_vector);
            std::vector<std::vector<double>> points_inside = in_out[0];
            std::vector<std::vector<double>> points_outside = in_out[1];

            int inside = points_inside[0].size();
            int outside = points_outside[0].size();

            double ratio = double(inside) / double(inside + outside);
            double pi = ratio * 4;
            relative_error[i] += fabs(pi - pi_actual) / pi_actual;
        }
    }
    //averages error and converts it to log scale for better visualisation
    for (int i = 0; i < length; i++){
        relative_error[i] /= num_runs;
        relative_error[i] = log2(relative_error[i]);
    }

    //quick calculation of the slope of the log-log plot to estimate order of 
    double slope = ((relative_error[length-1] - relative_error[0]) / (logn[length-1] - logn[0]));
    std::cout << "Slope: " << slope << std::endl;

    //plots the error vs n on a log2-log2 scale
    sciplot::Plot2D plot;
    plot.drawCurve(logn, relative_error).lineType(1).lineWidth(2).lineColor("black");
    plot.xlabel("log2(n)");
    plot.ylabel("log2(Error)");
    plot.legend().hide();
    sciplot::Figure figure = {{plot}};
    sciplot::Canvas canvas = {{figure}};
    canvas.show();
}

int main() {
    //plots points and errors for 1000 and [2, 2^20] points, respectively
    plot_points(1000);
    //plot_error(15);
}