#ifndef COUNTRY_H
#define COUNTRY_H

#include <vector>
#include <functional>
#include <limits>

class Country {
public:
    std::vector<double> location;
    double fitness;
    Country* vassal_of_empire;
    int index_in_list;
    double norm_imperialist_power;
    std::vector<Country*> vassals;
    std::vector<double> colour;

    // Constructor
    Country(const std::vector<double>& loc);

    // Destructor
    ~Country() = default;

    // Evaluate fitness using a provided objective function
    void evaluate_fitness(const std::function<double(const std::vector<double>&)>& objective_function);

    // Remove and return the weakest vassal
    Country* weakest_vassal_removal();

    // Add a vassal
    void add_vassal(Country* vassal);

    // Set emperor
    void add_emperor(Country* emperor);
};

#endif