#ifndef PICA_MP_H
#define PICA_MP_H

#include "ica.h"
#include <functional>
#include <vector>

class PICA_MP 
{
private:
    int rank, size;
    int dim;
    int pop_size;
    int max_iter;

    Imperialist_Competitive_Algorithm* ica;
    std::function<double(const std::vector<double>&)> obj_func;

    void print_results(double fitness, std::vector<double>& location);   
public:
    PICA_MP(int pop_size, int dim, int max_iter, const std::function<double(const std::vector<double>&)>& obj_func);

    void run(double lb, double ub, double beta, double gamma, double eta, int migration_cycles, int iterations_per_cycle);
    std::vector<double> get_best_solution() const;
    double get_best_fitness() const;

};

#endif