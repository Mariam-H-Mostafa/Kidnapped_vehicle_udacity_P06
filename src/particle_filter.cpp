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
#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  std::default_random_engine gen;
  num_particles = 100;  // TODO: Set the number of particles
 
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  double sample_x, sample_y, sample_theta;

  for(int i = 0;i<num_particles;i++)
  {
   
   Particle particle;
   sample_x = dist_x(gen);
   sample_y = dist_y(gen);
   sample_theta = dist_theta(gen);

   particle.id = i;
   particle.x = sample_x;
   particle.y = sample_y;
   particle.theta = sample_theta;
   particle.weight = 1.0;
   particles.push_back(particle);
   weights.push_back(1.0);
  }
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    std::default_random_engine gen;
   
   
  for(int i =0;i<num_particles;i++)
  {
    if(fabs(yaw_rate) < 0.0001)
    {
    
      particles[i].x =  particles[i].x + velocity*delta_t*cos(particles[i].theta);
      particles[i].y =  particles[i].y + velocity*delta_t*sin(particles[i].theta);
      particles[i].theta = particles[i].theta;
    }
    else
    {
     particles[i].x =  particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta+(yaw_rate*delta_t))-sin(particles[i].theta));
     particles[i].y =  particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta)-cos(particles[i].theta+(yaw_rate*delta_t)));
     particles[i].theta = particles[i].theta + yaw_rate*delta_t;
      
    }
   
   normal_distribution<double> dist_xf(particles[i].x, std_pos[0]);
   normal_distribution<double> dist_yf(particles[i].y, std_pos[1]);
   normal_distribution<double> dist_thetaf(particles[i].theta, std_pos[2]);
   particles[i].x =dist_xf(gen);
   particles[i].y = dist_yf(gen);
   particles[i].theta =dist_thetaf(gen);
  }  

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
 
  double diff;
  LandmarkObs Obs;
  LandmarkObs pr;
  double nearest = std::numeric_limits<double>::max();
  int index;
  for(unsigned int i =0;i<observations.size();i++)
  {
    Obs = observations[i];
    for(unsigned int j =0;j<predicted.size();j++)
    {
      pr = predicted[j];
      
      diff = dist(Obs.x,Obs.y,pr.x,pr.y);
     
      if(diff<nearest)
      {
        nearest = diff;
        index= predicted[j].id;
      }
    }
    observations[i].id = index;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
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
  vector<LandmarkObs> trans_observations;
  vector<LandmarkObs> predictions;
  double gauss_norm;
  double exponent;
  double mu_x;
  double mu_y;
  double obs_x;
  double obs_y;
  double distance;
  int idx;
  vector<double> sense_x;
  vector<double> sense_y;
  long double weight_perObs;
  vector<int> associations;
  
  double W_S = 0.0;
  gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
  
  for(int i = 0;i<num_particles;i++)
  {
    trans_observations.clear();
    for(unsigned int j=0;j<observations.size();j++)
    { 
        LandmarkObs trans_obs;
        trans_obs.x = particles[i].x + (cos(particles[i].theta))*(observations[j].x) - (sin(particles[i].theta))*(observations[j].y);

        trans_obs.y = particles[i].y + (sin(particles[i].theta))*(observations[j].x) + (cos(particles[i].theta))*(observations[j].y);
        trans_obs.id = observations[j].id;
        trans_observations.push_back(trans_obs);

    }
    predictions.clear();
    for(unsigned int l = 0; l<map_landmarks.landmark_list.size(); l++)
    {

      double lm_x = map_landmarks.landmark_list[l].x_f;
      double lm_y = map_landmarks.landmark_list[l].y_f;
      int lm_id = map_landmarks.landmark_list[l].id_i;
      distance = dist(lm_x,lm_y,particles[i].x,particles[i].y);
      if(distance<=sensor_range)
      {
        predictions.push_back(LandmarkObs{ lm_id, lm_x, lm_y });
        associations.push_back(lm_id);   
        sense_x.push_back(lm_x);
        sense_y.push_back(lm_y);
      }
     
    }
    SetAssociations(particles[i], associations, sense_x, sense_y);
    dataAssociation(predictions, trans_observations);
    
    particles[i].weight =1.0;
    for(unsigned int k =0;k<trans_observations.size();k++){   
      int id_search = trans_observations[k].id;
      obs_x = trans_observations[k].x;
      obs_y = trans_observations[k].y;
      auto itr = std::find_if(predictions.begin(),predictions.end(),[id_search](const LandmarkObs& landmark)
                            {return landmark.id==id_search;});
      if(itr!=predictions.cend()){
            idx = std::distance(predictions.begin(),itr);
       }
           mu_x = predictions[idx].x;
           mu_y = predictions[idx].y; 

          
          exponent = (pow(obs_x - mu_x, 2) / (2 * pow(std_landmark[0], 2)))
               + (pow(obs_y - mu_y, 2) / (2 * pow(std_landmark[1], 2)));

          weight_perObs =gauss_norm*exp(-exponent);
          if(weight_perObs>0)
          {
            particles[i].weight*=weight_perObs;
          }       
}
 W_S+=particles[i].weight;

}
  for(unsigned int s = 0;s<particles.size();s++)
  {
    particles[s].weight = particles[s].weight/W_S;
    weights[s] = particles[s].weight;

  }
}
void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
    std::vector<Particle> particles_resampled;
    std::random_device rd;
    std::mt19937 gen(rd());
  

    std::discrete_distribution<> d(weights.begin(),weights.end());

    for(int n=0; n<num_particles; ++n) {
        particles_resampled.push_back(particles[d(gen)]);
    }
  particles = particles_resampled;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations.clear();
  particle.sense_x.clear();
  particle.sense_y.clear();
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}