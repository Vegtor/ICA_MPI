#include "ica.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <random>

Imperialist_Competitive_Algorithm::Imperialist_Competitive_Algorithm(
    int pop_size, int dim, int max_iter,
    const std::function<double(const std::vector<double>&)>& obj_func)
    : pop_size(pop_size), dim(dim), max_iter(max_iter), obj_func(obj_func), best_fitness(INFINITY), tp(-1)
{
}

std::vector<double> Imperialist_Competitive_Algorithm::random_colour()
{
    static std::default_random_engine gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return { dist(gen), dist(gen), dist(gen) };
}

std::vector<std::vector<double>> Imperialist_Competitive_Algorithm::random_colours(int n)
{
    std::vector<std::vector<double>> result;
    for (int i = 0; i < n; ++i)
        result.push_back(random_colour());
    return result;
}

void Imperialist_Competitive_Algorithm::calculate_fitness()
{
    //different for cycle (auto with iteration)
    for (auto& c : population)
    {
        c->evaluate_fitness(obj_func);
        if (c->fitness < best_fitness)
        {
            best_fitness = c->fitness;
            best_solution = c->location;
        }
    }
}

void Imperialist_Competitive_Algorithm::create_empires()
{
    int n_empire = static_cast<int>(0.1 * pop_size);
    //change cycle for assign
    empires.assign(population.begin(), population.begin() + n_empire);
    colours = random_colours(n_empire);

    for (size_t i = 0; i < n_empire; ++i)
        empires[i]->colour = colours[i];

    double total_power = 0;
    //
    for (auto& e : empires)
        total_power += std::abs(e->fitness);
    this->tp = total_power;
}

void Imperialist_Competitive_Algorithm::create_colonies()
{
    // change for assign
    colonies.assign(population.begin() + empires.size(), population.end());
    std::vector<int> colonies_in_empires(empires.size());
    size_t assigned = 0;

    for (size_t i = 0; i < empires.size(); ++i)
    {
        colonies_in_empires[i] = static_cast<int>(std::floor(std::abs(empires[i]->fitness / this->tp) * colonies.size()));
        assigned += colonies_in_empires[i];
    }

    int diff = colonies.size() - assigned;
    for (size_t i = 0; diff > 0 && i < empires.size(); ++i)
    {
        int to_add = static_cast<int>(std::ceil(diff * std::abs(empires[i]->fitness / this->tp)));
        colonies_in_empires[i] += to_add;
        diff -= to_add;
    }

    //shuffle using std
    std::vector<int> indices(colonies.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(std::random_device{}()));

    int idx = 0;
    for (size_t i = 0; i < empires.size(); ++i)
    {
        for (size_t j = 0; j < colonies_in_empires[i] && idx < colonies.size(); ++j, ++idx)
        {
            // reafactoring
            Country* colony = colonies[indices[idx]];
            empires[i]->add_vassal(colony);
            colony->add_emperor(empires[i]);
            colony->colour = empires[i]->colour;
        }
    }
}

void Imperialist_Competitive_Algorithm::assimilation(double beta)
{
    for (auto& colony : colonies)
    {
        Country* emperor = colony->vassal_of_empire;
        double dist = 0;

        for (size_t i = 0; i < dim; ++i)
            dist += std::pow(emperor->location[i] - colony->location[i], 2);
        dist = std::sqrt(dist);


        if (dist != 0)
        {
            // if distance is 0, than shift is not needed (does it really make difference?)
            double shift = ((double)rand() / RAND_MAX) * beta * dist;
            for (size_t i = 0; i < dim; ++i)
                //directions are calculated every cycle, not ot itself
                colony->location[i] += shift * (emperor->location[i] - colony->location[i]) / dist;
        }
    }
}

void Imperialist_Competitive_Algorithm::revolution(double gamma)
{
    // just refactoring
    for (auto& colony : colonies)
    {
        for (size_t i = 0; i < dim; ++i)
            colony->location[i] += ((double)rand() / RAND_MAX) * 2 * gamma - gamma;
    }
}

void Imperialist_Competitive_Algorithm::mutiny()
{
    //
    for (Country* colony : colonies)
    {
        // using std function for minimum
        Country* nearest_imperialist = *std::min_element(empires.begin(), empires.end(), [&](Country* a, Country* b)
            {
                double country_a = 0;
                double country_b = 0;
                for (size_t i = 0; i < dim; ++i)
                {
                    country_a += std::pow(a->location[i] - colony->location[i], 2);
                    country_b += std::pow(b->location[i] - colony->location[i], 2);
                }
                return country_a < country_b;
            });

        if (colony->vassal_of_empire != nearest_imperialist)
        {
            std::vector<Country*>& vassals = colony->vassal_of_empire->vassals;
            auto idx = std::find(vassals.begin(), vassals.end(), colony);
            if (idx != vassals.end())
                vassals.erase(idx);
        }

        if (colony->fitness < nearest_imperialist->fitness)
        {
            colony->colour = nearest_imperialist->colour;
            colony->vassals = nearest_imperialist->vassals;
            colony->add_vassal(nearest_imperialist);
            nearest_imperialist->add_emperor(colony);
            colony->vassal_of_empire = nullptr;

            auto empire_idx = std::find(empires.begin(), empires.end(), nearest_imperialist);
            *empire_idx = colony;

            auto colony_idx = std::find(colonies.begin(), colonies.end(), colony);
            *colony_idx = nearest_imperialist;
        }
        else
        {
            colony->add_emperor(nearest_imperialist);
            colony->colour = nearest_imperialist->colour;
            nearest_imperialist->add_vassal(colony);
        }

    }
}

void Imperialist_Competitive_Algorithm::imperial_war(double eta)
{
    std::vector<double> total_power(empires.size());
    double max_power = -INFINITY;

    for (size_t i = 0; i < empires.size(); ++i)
    {
        double vassal_sum = 0;
        for (auto& v : empires[i]->vassals)
            vassal_sum += v->fitness;
        total_power[i] = empires[i]->fitness + eta * vassal_sum;
        if (total_power[i] > max_power)
            max_power = total_power[i];
    }

    std::vector<double> normalized_powers;
    double sum_norm_power = 0;
    for (auto power : total_power)
    {
        normalized_powers.push_back(power - max_power);
        sum_norm_power += power - max_power;
    }

    std::vector<double> D;
    for (auto norm : normalized_powers)
        D.push_back(norm / sum_norm_power - ((double)rand() / RAND_MAX));

    // using std functions
    int weakest_emp_idx = std::min_element(D.begin(), D.end()) - D.begin();
    int strongest_emp_idx = std::max_element(D.begin(), D.end()) - D.begin();

    if (!empires[weakest_emp_idx]->vassals.empty())
    {
        Country* weakest_vassal = empires[weakest_emp_idx]->weakest_vassal_removal();
        weakest_vassal->colour = empires[strongest_emp_idx]->colour;
        empires[strongest_emp_idx]->add_vassal(weakest_vassal);
        weakest_vassal->add_emperor(empires[strongest_emp_idx]);
    }
    else
    {
        empires[strongest_emp_idx]->add_vassal(empires[weakest_emp_idx]);
        empires[weakest_emp_idx]->add_emperor(empires[strongest_emp_idx]);
        colonies.push_back(empires[weakest_emp_idx]);
        empires.erase(empires.begin() + weakest_emp_idx);
    }
}

void Imperialist_Competitive_Algorithm::optimize(double lb, double ub, double beta, double gamma, double eta)
{
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    for (size_t i = 0; i < pop_size; ++i)
    {
        std::vector<double> pos(dim);
        for (int j = 0; j < dim; ++j)
        {
            double r = static_cast<double>(std::rand()) / (RAND_MAX + 1.0);
            pos[j] = lb + r * (ub - lb);
        }
        population.push_back(new Country(pos));
    }

    calculate_fitness();
    std::sort(population.begin(), population.end(), [](Country* a, Country* b)
        {
            return a->fitness < b->fitness;
        });
    create_empires();
    create_colonies();

    this->alg_run(beta, gamma, eta);
}

void Imperialist_Competitive_Algorithm::alg_run(double beta, double gamma, double eta)
{
    for (int i = 0; i < max_iter; ++i)
    {
        calculate_fitness();
        //fitness_history.push_back(best_fitness);
        assimilation(beta);
        revolution(gamma);
        mutiny();
        imperial_war(eta);
        if (empires.size() == 1)
            break;
    }
}

void Imperialist_Competitive_Algorithm::migrate_best(const std::vector<double>& elite_solution, const std::function<double(const std::vector<double>&)>& obj_func)
{
    auto worst = std::max_element(population.begin(), population.end(), [](Country* a, Country* b)
        {
            return a->fitness < b->fitness;
        });

    Country* worst_country = *worst; 
    worst_country->location = elite_solution;
    worst_country->evaluate_fitness(obj_func);
}

double Imperialist_Competitive_Algorithm::get_fitness()
{
    return this->best_fitness;
}

std::vector<double> Imperialist_Competitive_Algorithm::get_best_solution()
{
    return this->best_solution;
}

void Imperialist_Competitive_Algorithm::set_max_iter(int max_iter)
{
    this->max_iter = max_iter;
}

void Imperialist_Competitive_Algorithm::check(int rank)
{
    for (size_t i = 0; i < colonies.size(); ++i)
    {
        const Country* colony = colonies[i];
        if (!colony)
        {
            std::cerr << "Rank " << rank << "Error: colony at index " << i << " is nullptr!\n";
            break;
        }

        if (!colony->vassal_of_empire)
        {
            std::cerr << "Rank " << rank << "Error: colony at index " << i << " has null vassal_of_empire pointer!\n";
            break;
        }

        if (colony->location.size() != dim)
        {
            std::cerr << "Rank " << rank << "Error: colony at index " << i << " location vector size mismatch. Expected "
                << dim << ", got " << colony->location.size() << "\n";
            break;
        }
    }
}

Imperialist_Competitive_Algorithm::~Imperialist_Competitive_Algorithm()
{
    for (auto c : population) 
        delete c;
}
