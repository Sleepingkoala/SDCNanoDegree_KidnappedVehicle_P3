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

///* declare a random engine for generating random numbers
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
      ///* set particle number and random noise gaussian distributions
       num_particles = 101;
       
       normal_distribution<double> dist_x_noise(0, std[0]);
       normal_distribution<double> dist_y_noise(0, std[1]);
       normal_distribution<double> dist_theta_noise(0, std[2]);
       
       for (int i = 0; i< num_particles; i++){
           Particle p;
           ///* initialize all particles to the first position and all weights to 1
           p.id = i;
           p.x = x;
           p.y = y;
           p.theta = theta;
           p.weight = 1.0;
           
           ///* add random gaussian noise to each particle
           p.x += dist_x_noise(gen);
           p.y += dist_y_noise(gen);
           p.theta += dist_theta_noise(gen);
          
           ///* save each particle into particles set
           particles.push_back(p);
       }
       is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
        
        ///* set random noise gaussian distributions
        normal_distribution<double> dist_x_noise(0, std_pos[0]);
        normal_distribution<double> dist_y_noise(0, std_pos[1]);
        normal_distribution<double> dist_theta_noise(0, std_pos[2]);
        
        ///* add measurent(velocity, yawrate) to each particle
        for(int i = 0; i< num_particles; i++){
        
        if(fabs(yaw_rate)> 0.0001){
            particles[i].x += velocity/ yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
            particles[i].y += velocity/ yaw_rate * (-cos(particles[i].theta + yaw_rate * delta_t) + cos(particles[i].theta));
            particles[i].theta +=  yaw_rate* delta_t;
        }
        else{
            particles[i].x +=  velocity * delta_t *cos(particles[i].theta);
            particles[i].y +=  velocity * delta_t * sin(particles[i].theta);
        }
        ///* add random gaussian noises to each particle
        particles[i].x += dist_x_noise(gen);
        particles[i].y += dist_y_noise(gen);
        particles[i].theta += dist_theta_noise(gen);        
        }         	        
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
        
        ///* for each landmark observation to search all the predicted landmarks to find the closest one 
        ///* then has the closest predicted landmark id and data association is bonded.
        for (int i = 0; i < observations.size(); i++){
            
            LandmarkObs obs = observations[i];
            
            double min_distance = numeric_limits<double>::max(); // set initial distance minimum to possible maximum
            int map_id = -1; // set initial map landmark id from map placeholder to observation association
            
            for (int j = 0; j< predicted.size(); j++){
                LandmarkObs pred = predicted[j];
                
                double distance = dist(pred.x, pred.y, obs.x, obs.y);
                if (distance < min_distance){
                    min_distance = distance;
                    map_id = pred.id;
                }
            }
            
            observations[i].id = map_id;
        }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
	
	///* transformation and associates for each particle
	for (int i = 0; i< num_particles; i++){
	
	    double p_x = particles[i].x;
	    double p_y = particles[i].y;
	    double p_theta = particles[i].theta;
	    
	    vector<LandmarkObs> predictions;	    
	    for (int j = 0; j < map_landmarks.landmark_list.size(); j++){
	    
	        float landmark_x = map_landmarks.landmark_list[j].x_f;
	        float landmark_y = map_landmarks.landmark_list[j].y_f;
	        float landmark_id = map_landmarks.landmark_list[j].id_i;
	        
	        if(fabs(p_x - landmark_x) <=sensor_range && fabs(p_y - landmark_y) <= sensor_range){	        
	            predictions.push_back(LandmarkObs{landmark_id, landmark_x, landmark_y});
	        }
	    }
	   
	   
	   ///**** 1. convert the observations from vehicle coordinate system to the map coordinate system
	   vector<LandmarkObs> transformed_obs;
	   for(int j = 0; j< observations.size(); j++){	       
	       double transformed_x = cos(p_theta) * observations[j].x - sin(p_theta)*observations[j].y + p_x;
	       double transformed_y = sin(p_theta) * observations[j].x + cos(p_theta)* observations[j].y + p_y;
	       transformed_obs.push_back(LandmarkObs{observations[j].id, transformed_x, transformed_y});
	   }
	   
	   ///****2. associate each transformed observation with a predicted landmark id
	   dataAssociation(predictions, transformed_obs);
	   
	   particles[i].weight = 1;  // re-initialization to 1 for the following final weight product 
	   for (int j = 0; j< transformed_obs.size(); j++){
	       double observation_x, observation_y, prediction_x, prediction_y;
	       observation_x = transformed_obs[j].x;
	       observation_y = transformed_obs[j].y;
	       
	       int associated_prediction_id = transformed_obs[j].id;
	       for (int k = 0; k < predictions.size(); k++){
	           if(predictions[k].id == associated_prediction_id){
	               prediction_x = predictions[k].x;
	               prediction_y = predictions[k].y;
	           }
	       }
	       
	       ///**** 3. calculate the final weight of this particle
	       double std_x = std_landmark[0];
	       double std_y = std_landmark[1];
	       double observation_weight = (1/(2*M_PI* std_x* std_y)) * 
	                                                             exp(-(pow(observation_x - prediction_x, 2)/(2*pow(std_x,2)) + pow(observation_y-prediction_y,2)/(2*pow(std_y,2))));
	       
	       particles[i].weight *= observation_weight;	    
	   }	   
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
        vector<Particle> new_particles;
        vector<double> weights;
        for (int i = 0; i< num_particles; i++){
            weights.push_back(particles[i].weight);
        }
        
        ///* use the resampling wheel for resampling
        uniform_int_distribution<int> dist_uni_int(0, num_particles -1);
        auto index = dist_uni_int(gen);
        
        double max_weight = *max_element(weights.begin(), weights.end());
        
        uniform_real_distribution<double> dist_uni_real(0.0, max_weight);
        
        double beta = 0.0;
        
        for (int i = 0; i< num_particles; i++){
            beta += dist_uni_real(gen) * 2.0;
            while(beta > weights[index]){
                beta -= weights[index];
                index = (index + 1) % num_particles;
            }
            new_particles.push_back(particles[index]);
        }
        particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

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
