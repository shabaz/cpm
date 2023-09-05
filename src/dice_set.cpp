#include <iostream>
#include "dice_set.h"

using namespace std;

void DiceSet::add(int element) {
    if (_map.count(element) == 0) {
        int index = _vector.size();
        _map[element] = index;
        _vector.push_back(element);
    }
}

void DiceSet::remove(int element) {
    if (_map.count(element) == 0)
        return;
    int index = _map[element];
    int lastElement = _vector.back();
    _vector[index] = lastElement;
    _map[lastElement] = index;
    _vector.pop_back();
    _map.erase(element);
}

int DiceSet::get(int index) {
    return _vector[index];
}

int DiceSet::size() {
    return _vector.size();
}
