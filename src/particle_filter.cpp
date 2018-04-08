/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    num_particles = 100;
    default_random_engine gen;
    
    // This line creates a normal (Gaussian) distribution for x
    normal_distribution<double> dist_x(x, std[0]);
    
    // TODO: Create normal distributions for y and theta
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    
    for (int i=0; i < num_particles; i++) {
        double sample_x, sample_y, sample_theta;
        
        sample_x = dist_x(gen);
        sample_y = dist_y(gen);
        sample_theta = dist_theta(gen);
        Particle sample_particle = {
            .id = i,
            .x =  sample_x,
            .y = sample_y,
            .theta = sample_theta,
            .weight = 1.0 / num_particles
        };
        
        particles.push_back(sample_particle);
        
    }
    
    
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    default_random_engine gen;
    
    for (int i = 0; i < num_particles; i++) {
        if(abs(yaw_rate) < 0.0001)
        {
            particles[i].x = particles[i].x + (delta_t*velocity*(cos(particles[i].theta)));
            particles[i].y = particles[i].y + (delta_t*velocity*(sin(particles[i].theta)));
        } else {
            double predicted_theta = particles[i].theta + yaw_rate * delta_t;
            particles[i].x += (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
            particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
            particles[i].theta = predicted_theta;
            
            
            normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
            normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
            normal_distribution<double> dist_t(particles[i].theta, std_pos[2]);
            
            particles[i].x = dist_x(gen);
            particles[i].y = dist_y(gen);
            particles[i].theta = dist_t(gen);
        }
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    
    // assoc vector holds the association of each item i
    // assoc_val vector holds the distance for each association
    double assoc[observations.size()];
    double assoc_val[observations.size()];
    for (int x = 0; x < observations.size(); x++) {
        assoc[x] = 0;
        assoc_val[x] = INFINITY;
        observations[x].id = 0;
    }
    
    for (int j = 0; j < observations.size(); j++) {
        for (int i = 0; i < predicted.size(); i++) {
            double distance = sqrt(((predicted[i].x - observations[j].x) * (predicted[i].x - observations[j].x)) + ((predicted[i].y - observations[j].y) * (predicted[i].y - observations[j].y)));
            if (distance <= assoc_val[j]) {
                assoc[j] = i;
                assoc_val[j] = distance;
                observations[j].id = i;
            }
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
        std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    
    double s_x = sqrt(std_landmark[0]);
    double s_y = sqrt(std_landmark[1]);
    double C = (1 / (2 * M_PI*s_x * s_y));
    
    for (int i = 0; i < num_particles; i++) {
        vector<LandmarkObs> transform_obs;
        
        for (int j = 0; j < observations.size(); j++) {
            LandmarkObs landmark;
            landmark.x = (observations[j].x * cos(particles[i].theta)) - (observations[j].y * sin( particles[i].theta)) + particles[i].x;
            landmark.y = (observations[j].x * sin(particles[i].theta)) + (observations[j].y * cos( particles[i].theta)) + particles[i].y;
            transform_obs.push_back(landmark);
        }
        
        //getting landmarks in range
        vector<LandmarkObs> predicted;
        int id = 0;
        for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            double distance = sqrt(((map_landmarks.landmark_list[j].x_f - particles[i].x) * (map_landmarks.landmark_list[j].x_f - particles[i].x)) + ((map_landmarks.landmark_list[j].y_f - particles[i].y) * (map_landmarks.landmark_list[j].y_f - particles[i].y)));
            if (distance <= sensor_range) {
                LandmarkObs landmark;
                landmark.id = id;
                landmark.x = map_landmarks.landmark_list[j].x_f;
                landmark.y = map_landmarks.landmark_list[j].y_f;
                predicted.push_back(landmark);
                id++;
            }
        }
        
        double weight = 1.0;
        
        if (predicted.size() > 0) {
            dataAssociation(predicted,transform_obs);
            for (int j = 0; j < transform_obs.size(); j++) {
                weight = weight * (C*exp(-((((predicted[transform_obs[j].id].x - transform_obs[j].x)*(predicted[transform_obs[j].id].x - transform_obs[j].x))) +
                                         (((predicted[transform_obs[j].id].y - transform_obs[j].y)*(predicted[transform_obs[j].id].y - transform_obs[j].y)))) / (2 * s_y * s_y)));
            }
            particles[i].weight = weight;
        } else {
            particles[i].weight = 1 / num_particles;
        }
    }
    
    double sum = 0;
    for (int i = 0; i < num_particles; i++) {
        sum+= particles[i].weight;
    }
    // normalizing weights
    if (sum > 0 ) {
        for (int i = 0; i < num_particles; i++) {
            particles[i].weight = particles[i].weight / sum;
        }
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    vector<double> weights;
    for (int i = 0; i < num_particles; i++) {
        weights.push_back(particles[i].weight);
    }
    
    std::discrete_distribution<> d(weights.begin(),weights.end());
    for (int i = 0; i < num_particles; i++) {
        int x = d(gen);
        particles[i] = particles[x];
    }
    
    double sum = 0;
    for (int i = 0; i < num_particles; i++) {
        sum+= particles[i].weight;
    }
    // normalizing weights
    if (sum > 0 ) {
        for (int i = 0; i < num_particles; i++) {
            particles[i].weight = particles[i].weight / sum;
        }
    }
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    
    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
