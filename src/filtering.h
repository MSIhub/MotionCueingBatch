// Implement Filtering.h

#ifndef FILTERING_H
#define FILTERING_H

#include <iostream>
#include <sstream>// ostringstream for logging
#include <chrono>
#include <iomanip>//setprecision
#include <fstream>//reading files
#include <string>
#include <iterator>
#include <map>
#include <cstdint>//int header
#include<math.h>
#include <ctime>


//MACROS: CUEING
# define M_PI           3.14159265358979323846  /* pi */
#define KERNEL_LENGTH 4	
constexpr double SP7_ZERO_POSE[6] = { 0.0, 0.0, 0.401, 0.0, 0.0, 0.0 };

namespace filtering
{
    static double initial_time = 0.0;
    #pragma region STRUCTURE
    struct SP7Pose
    {
       double x, y, z, roll, pitch, yaw;
    };

	struct SP7Vel
	{
		double vx, vy, vz, vroll, vpitch, vyaw;
	};


	struct CueData
	{
		double t_prev = 0.0;
		double velocity_prev = 0.0;
		double position_prev = 0.0;
		double acc_fltrd_scaled_prev = 0.0;

		double t = 0.0;
		double acc_fltrd = 0.0;
		double acc_fltrd_scaled = 0.0;
		double velocity = 0.0;
		double position = 0.0;

		double Input_Buff[KERNEL_LENGTH];
		int circ_buff_idx = KERNEL_LENGTH - 1;
	};


	struct CueDataVel
	{
		double t_prev = 0.0;
		double velocity_fltr_scaled_prev = 0.0;
		double position_prev = 0.0;

		double t = 0.0;
		double velocity_fltr = 0.0;
		double velocity_fltr_scaled = 0.0;
		double position = 0.0;

		double Input_Buff[KERNEL_LENGTH];
		int circ_buff_idx = KERNEL_LENGTH - 1;
	};

#pragma endregion


	//Filter Kernel Calculations
	void calc_kernel_high_pass(int filter_length, double cut_off, double* fltr_kernel_arr_, bool);
	void calc_kernel_low_pass(int filter_length, double cut_off, double* fltr_kernel_arr_, bool isLogging);


	//FUNCTION PROTOTYPES 
	void cue_translational_channel(double sig_acc_input, double sig_time, double* kernel, double scale_factor, int data_index, CueData* cue_data, double* out_pos_, double* out_vel_, double* out_t_);
	void cue_rotational_channel(double sig_vel_input, double sig_time, double* kernel, double scale_factor, int data_index, CueDataVel* cue_data, double* out_pos_, double* out_vel_, double* out_t_);
	void cue_tilt_coordination_channel(double sig_acc_input, double sig_time, double* kernel, double scale_factor, int  data_index, CueData* cue_data, double* out_ang_, double* out_t_);
	double Intergration_Trapezoidal(double input_curr, double input_prev, double output_prev, double t_prev, double t_curr);
	double Convolve_rt(double* h, int h_size, double x_in, double* x, int* circ_index);

}

#endif
