#ifndef ICA_H
#define ICA_H

#include "Country.h"
#include <vector>
#include <functional>

class Imperialist_Competitive_Algorithm
{
private:
    int pop_size, dim, max_iter;
    std::function<double(const std::vector<double>&)> obj_func;
    std::vector<Country*> population, empires, colonies;
    std::vector<std::vector<double>> colours;
    std::vector<double> fitness_history, best_solution;
    double best_fitness;
    double tp;

    // colour generating - mostly for ploting
    std::vector<double> random_colour();
    // vector of colour generation
    std::vector<std::vector<double>> random_colours(int n);
    // country/empire fitness calculation
    void calculate_fitness();
    // setup of n first countries as emperors
    void create_empires();
    // setup of colonies to empires
    void create_colonies();
    // process of assimilation of small empires
    void assimilation(double beta);
    // movement of colony
    void revolution(double gamma);
    // revolt against emperor
    void mutiny();
    // war between empires - colonies shift between empires
    void imperial_war(double eta);
public:
    //setup of ICA
    Imperialist_Competitive_Algorithm(int pop_size, int dim, int max_iter, const std::function<double(const std::vector<double>&)>& obj_func);

    // algorithm run
    void optimize(double lb, double ub, double beta, double gamma, double eta);

    // algorithm clean run - without setup
    void alg_run(double beta, double gamma, double eta);

    // adding best empire from different proces
    void migrate_best(const std::vector<double>& elite_solution, const std::function<double(const std::vector<double>&)>& obj_func);

    //getters for fitness, solution coordinates
    double get_fitness();
    std::vector<double> get_best_solution();

    //setter for max iterations
    void set_max_iter(int max_iter);

    void check(int rank);

    ~Imperialist_Competitive_Algorithm();
};

#endif // ICA_H