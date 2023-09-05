#ifndef DICE_SET_H_
#define DICE_SET_H_

#include <unordered_map>
#include <vector>
#include "bytell_hash_map.hpp"

class DiceSet {
    public:
        void add(int element);
        void remove(int element);
        int get(int index);
        int size();
        //std::unordered_map<int, int> _map;
        ska::bytell_hash_map<int, int> _map;
    private:
        std::vector<int> _vector;
};

#endif // DICE_SET_H_
