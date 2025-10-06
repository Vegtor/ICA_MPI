#include "pica_mp.h"
#include <mpi.h>
#include <iostream>
#include <cassert>

PICA_MP::PICA_MP(int pop_size, int dim, int max_iter,
    const std::function<double(const std::vector<double>&)>& obj_func)
    : pop_size(pop_size), dim(dim), max_iter(max_iter), obj_func(obj_func) 
{

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    ica = new Imperialist_Competitive_Algorithm(pop_size, dim, max_iter, obj_func);
}

void PICA_MP::print_results(double fitness, std::vector<double>& location)
{
    std::cout << "Best fitness: " << fitness << std::endl;
    std::cout << "Best solution location: [";
    for (size_t i = 0; i < location.size(); ++i)
    {
        std::cout << location[i];
        if (i != location.size() - 1)
            std::cout << ", ";
    }
    std::cout << "]" << std::endl;
}


void PICA_MP::run(double lb, double ub, double beta, double gamma, double eta, int migration_cycles, int iterations_per_cycle)
{
    ica->optimize(lb, ub, beta, gamma, eta);
    ica->set_max_iter(iterations_per_cycle);
    MPI_Barrier(MPI_COMM_WORLD);

    for (int cycle = 0; cycle < migration_cycles; ++cycle) 
    {
        int prev = (rank - 1 + size) % size;
        int next = (rank + 1) % size;
        int tag = 0;

        std::vector<double> send_solution = ica->get_best_solution();
        std::vector<double> recv_solution(dim);

        for (int r = 0; r < size; ++r)
        {
            if (rank == r) 
            {
                std::cerr << "[Cycle " << cycle << "] Rank " << rank
                    << ": sending to " << next << ", receiving from " << prev
                    << ", send_solution size = " << send_solution.size() << std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }

        int err = MPI_Sendrecv(send_solution.data(), dim, MPI_DOUBLE, next, tag, recv_solution.data(), dim, MPI_DOUBLE, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        ica->migrate_best(recv_solution, obj_func);
        ica->alg_run(beta, gamma, eta);
    }

    std::vector<double> local_best = ica->get_best_solution();
    double local_fitness = ica->get_fitness();

    std::vector<double> all_best_solutions;
    std::vector<double> all_fitnesses;

    if (rank == 0)
    {
        all_best_solutions.resize(size * dim);
        all_fitnesses.resize(size);
    }

    MPI_Gather(local_best.data(), dim, MPI_DOUBLE, all_best_solutions.data(), dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_fitness, 1, MPI_DOUBLE, all_fitnesses.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) 
    {
        int best_index = 0;
        double best_fitness = all_fitnesses[0];
        for (int i = 1; i < size; ++i) 
        {
            if (all_fitnesses[i] < best_fitness) 
            {
                best_fitness = all_fitnesses[i];
                best_index = i;
            }
        }

        std::vector<double> global_best(dim);
        std::copy(all_best_solutions.begin() + best_index * dim, all_best_solutions.begin() + (best_index + 1) * dim, global_best.begin());
        this->print_results(best_fitness, global_best);
    }
}

std::vector<double> PICA_MP::get_best_solution() const 
{
    return ica->get_best_solution();
}

double PICA_MP::get_best_fitness() const 
{
    return ica->get_fitness();
}