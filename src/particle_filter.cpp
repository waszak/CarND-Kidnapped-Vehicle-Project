/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <cassert>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[])
{
    /**
     * TODO: Set the number of particles. Initialize all particles to
     *   first position (based on estimates of x, y, theta and their uncertainties
     *   from GPS) and all weights to 1.
     * TODO: Add random Gaussian noise to each particle.
     * NOTE: Consult particle_filter.h for more information about this method
     *   (and others in this file).
     */
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    gen = std::mt19937_64(rd());

    num_particles = 100;
    particles.clear();
    for(int i = 0; i < num_particles; ++i)
    {
        Particle particle;
        particle.id = i;
        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
        particle.weight = 1;
        particles.push_back(particle);
    }
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate)
{
    /**
     * TODO: Add measurements to each particle and add random Gaussian noise.
     * NOTE: When adding noise you may find std::normal_distribution
     *   and std::default_random_engine useful.
     *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
     *  http://www.cplusplus.com/reference/random/default_random_engine/
     */

    std::normal_distribution<double> dist_x(0, std_pos[0]);
    std::normal_distribution<double> dist_y(0, std_pos[1]);
    std::normal_distribution<double> dist_theta(0, std_pos[2]);

    for(auto & particle: particles)
    {
        //add measurement
        if( fabs(yaw_rate) > ERROR)
        {
            double c = (velocity / yaw_rate);
            double d = delta_t * yaw_rate;
            particle.x +=  c * (sin(particle.theta + d)- sin(particle.theta));
            particle.y +=  c * (cos(particle.theta) - cos(particle.theta + d));
            particle.theta += d;
        }
        else
        {
            double c = delta_t * velocity;
            particle.x += c * cos(particle.theta);
            particle.y += c * sin(particle.theta);
        }

        //add noise
        particle.theta += dist_theta(gen);
        particle.x += dist_x(gen);
        particle.y += dist_y(gen);
    }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predictions,
                                     vector<LandmarkObs>& observations)
{
    /**
     * TODO: Find the predicted measurement that is closest to each
     *   observed measurement and assign the observed measurement to this
     *   particular landmark.
     * NOTE: this method will NOT be called by the grading code. But you will
     *   probably find it useful to implement this method and use it as a helper
     *   during the updateWeights phase.
     */


    for(auto & observation: observations)
    {
        double min_dist = std::numeric_limits<double>::max();
        for(auto const& predicted : predictions)
        {
            double distance = dist(predicted.x, predicted.y, observation.x, observation.y);
            if( distance < min_dist )
            {
                min_dist = distance;
                observation.id = predicted.id;
            }
        }
    }

}

vector<LandmarkObs> ParticleFilter::toMapCoordinates(const Particle& particle, const vector<LandmarkObs> &observations)
{
    vector<LandmarkObs> observations_map;
    for(auto const& obs: observations)
    {
        double x_map = cos(particle.theta) * obs.x - sin(particle.theta) * obs.y;
        x_map += particle.x;
        double y_map = sin(particle.theta) * obs.x + cos(particle.theta) * obs.y;
        y_map += particle.y;
        LandmarkObs observation = LandmarkObs{-1, x_map, y_map};
        observations_map.push_back(observation);
    }
    return observations_map;
}

vector<LandmarkObs> ParticleFilter::findLandmarksInRange(const Particle& particle, double range, const Map &map_landmarks)
{
    vector<LandmarkObs> predicted;
    for(auto const& landmark: map_landmarks.landmark_list)
    {
        if(fabs(particle.x - landmark.x_f) > range || fabs(particle.y - landmark.y_f) > range)
        {
            continue;
        }
        LandmarkObs obs = LandmarkObs{landmark.id_i, landmark.x_f, landmark.y_f};
        predicted.push_back(obs);
    }
    return predicted;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks)
{
    /**
     * TODO: Update the weights of each particle using a mult-variate Gaussian
     *   distribution. You can read more about this distribution here:
     *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
     * NOTE: The observations are given in the VEHICLE'S coordinate system.
     *   Your particles are located according to the MAP'S coordinate system.
     *   You will need to transform between the two systems. Keep in mind that
     *   this transformation requires both rotation AND translation (but no scaling).
     *   The following is a good resource for the theory:
     *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
     *   and the following is a good resource for the actual equation to implement
     *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
     */
    double sig_x = std_landmark[0];
    double sig_y = std_landmark[1];
    auto landmarkMap = map_landmarks.getMap();
    for(auto& particle: particles)
    {

        vector<LandmarkObs> mapped_observation = toMapCoordinates(particle, observations);
        vector<LandmarkObs> predicted = findLandmarksInRange(particle, sensor_range, map_landmarks);

        dataAssociation(predicted, mapped_observation);

        vector<double> sense_x;
        vector<double> sense_y;
        vector<int> associations;
        //reset weight
        particle.weight = 1;
        if (predicted.size() == 0)
        {
            particle.weight = 0;//0 or some EPS^len(mapped_observation)
            continue;
        }
        for(auto const& observation: mapped_observation)
        {
            //for debuging purposes detecting -nan
            assert(observation.id != -1);
            auto predicted = landmarkMap.at(observation.id);

            double gauss_norm = 1/(2.0 * M_PI * sig_x * sig_y);
            double diff_x = predicted.x_f - observation.x;
            double diff_y = predicted.y_f - observation.y;

            double power1 = (diff_x * diff_x) / (2.0 * sig_x * sig_x);
            double power2 = (diff_y * diff_y) / (2.0 * sig_y * sig_y);
            double power = power1 + power2;

            double weight = gauss_norm * exp(-power) ;
            particle.weight *= weight;


            sense_x.push_back(predicted.x_f);
            sense_y.push_back(predicted.y_f);
            associations.push_back(predicted.id_i);
        }
        setAssociations( particle, associations, sense_x, sense_y);
    }
}

void ParticleFilter::resample()
{
    /**
     * TODO: Resample particles with replacement with probability proportional
     *   to their weight.
     * NOTE: You may find std::discrete_distribution helpful here.
     *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
     */
    vector<double> weights;
    //double total_weight = 0;
    for( auto& particle: particles)
    {
        //total_weight += particle.weight;
        weights.push_back(particle.weight);
    }
    //assert(total_weight != 0);

    std::discrete_distribution<> d(weights.begin(), weights.end());
    vector<Particle> temp_particles;
    for(int i = 0; i < num_particles; ++i)
    {
        int idx = d(gen);
        temp_particles.push_back(particles[idx]);
    }
    particles = temp_particles;

}

void ParticleFilter::setAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y)
{
    // particle: the particle to which assign each listed association,
    //   and association's (x,y) world coordinates mapping
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
    vector<int> v = best.associations;
    std::stringstream ss;
    copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord)
{
    vector<double> v;

    if (coord == "X")
    {
        v = best.sense_x;
    }
    else
    {
        v = best.sense_y;
    }

    std::stringstream ss;
    copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
