#include "Country.h"

// Constructor
Country::Country(const std::vector<double>& loc)
    : location(loc),
    fitness(-std::numeric_limits<double>::infinity()),
    vassal_of_empire(nullptr),
    index_in_list(-1),
    norm_imperialist_power(0), 
    colour(3, 0.0) 
{
}

// Evaluate fitness
void Country::evaluate_fitness(const std::function<double(const std::vector<double>&)>& objective_function) 
{
    fitness = objective_function(location);
}

// Remove and return the weakest vassal
Country* Country::weakest_vassal_removal() 
{
    if (vassals.empty()) return nullptr;

    Country* weakest = vassals[0];
    size_t index = 0;

    for (size_t i = 1; i < vassals.size(); ++i) {
        if (vassals[i]->fitness < weakest->fitness) {
            weakest = vassals[i];
            index = i;
        }
    }

    vassals.erase(vassals.begin() + index);
    return weakest;
}

// Add a vassal
void Country::add_vassal(Country* vassal) 
{
    vassals.push_back(vassal);
}

// Set emperor
void Country::add_emperor(Country* emperor) 
{
    vassal_of_empire = emperor;
}
