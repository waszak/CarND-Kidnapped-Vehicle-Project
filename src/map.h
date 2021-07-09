/**
 * map.h
 *
 * Created on: Dec 12, 2016
 * Author: mufferm
 */

#ifndef MAP_H_
#define MAP_H_

#include <vector>
#include <map>
class Map
{
public:
    struct single_landmark_s
    {
        int id_i ; // Landmark ID
        float x_f; // Landmark x-position in the map (global coordinates)
        float y_f; // Landmark y-position in the map (global coordinates)
    };

    std::vector<single_landmark_s> landmark_list; // List of landmarks in the map
    std::map<int, single_landmark_s> getMap() const
    {
        std::map<int, single_landmark_s> map_;
        for(const auto & landmark: landmark_list)
        {
            map_.emplace(landmark.id_i, landmark);
        }
        return map_;
    }
};

#endif  // MAP_H_
