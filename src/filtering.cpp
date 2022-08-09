#include "filtering.h"

void filtering::calc_kernel_high_pass(int filter_length, double cut_off, double* fltr_kernel_arr_, bool isLogging)
{
	double threshold = 0.0;
	//Calculate the first low-pass filter kernel 
	for (int i = 0; i < filter_length; i++)
	{
		threshold = (i - (double)(filter_length / 2));
		if (threshold == 0)
		{
			fltr_kernel_arr_[i] = 2 * M_PI * cut_off;
		}
		if (threshold != 0)
		{
			fltr_kernel_arr_[i] = sin(2 * M_PI * cut_off * threshold) / threshold;
			fltr_kernel_arr_[i] *= (0.42 - 0.5 * cos(2 * M_PI * i / filter_length) + 0.08 * cos(4 * M_PI * i / filter_length));
		}
	}

	//Change low-pass filter to a high-pass filter using spectral inversion
	for (int i = 0; i < filter_length; i++)
	{
		fltr_kernel_arr_[i] = -fltr_kernel_arr_[i];
	}
	fltr_kernel_arr_[filter_length / 2] += 1;

	if (isLogging)
	{
		std::string log_filename = "high_pass_kernel_";
		log_filename += std::to_string(std::time(nullptr));
		log_filename += ".dat";

		std::fstream output_sig_rt_kernel_fptr;
		output_sig_rt_kernel_fptr.open("log/" + log_filename, std::fstream::in | std::fstream::out | std::fstream::app);

		if (output_sig_rt_kernel_fptr.is_open())
		{
			for (int i = 0; i < filter_length; i++)
			{
				output_sig_rt_kernel_fptr << "\n";
				output_sig_rt_kernel_fptr << fltr_kernel_arr_[i];//column vector
			}
			output_sig_rt_kernel_fptr.close();
		}

	}

}

void filtering::calc_kernel_low_pass(int filter_length, double cut_off, double* fltr_kernel_arr_, bool isLogging)
{
	double threshold = 0.0;
	//Calculate the first low-pass filter kernel 
	for (int i = 0; i < filter_length; i++)
	{
		threshold = (i - (double)(filter_length / 2));
		if (threshold == 0)
		{
			fltr_kernel_arr_[i] = 2 * M_PI * cut_off;
		}
		if (threshold != 0)
		{
			fltr_kernel_arr_[i] = sin(2 * M_PI * cut_off * threshold) / threshold;
			fltr_kernel_arr_[i] *= (0.42 - 0.5 * cos(2 * M_PI * i / filter_length) + 0.08 * cos(4 * M_PI * i / filter_length));
		}
	}

	if (isLogging)
	{
		std::string log_filename = "low_pass_kernel_";
		log_filename += std::to_string(std::time(nullptr));
		log_filename += ".dat";

		std::fstream output_sig_rt_kernel_fptr;
		output_sig_rt_kernel_fptr.open("log" + log_filename, std::fstream::in | std::fstream::out | std::fstream::app);

		if (output_sig_rt_kernel_fptr.is_open())
		{
			for (int i = 0; i < filter_length; i++)
			{
				output_sig_rt_kernel_fptr << "\n";
				output_sig_rt_kernel_fptr << fltr_kernel_arr_[i];//column vector
			}
			output_sig_rt_kernel_fptr.close();
		}

	}

}

void filtering::cue_translational_channel(double sig_acc_input, double sig_time, double* kernel, double scale_factor, int data_index, filtering::CueData* cue_data, double* out_pos_, double* out_vel_, double* out_t_)
{
	/* filter => scale => integrate => integrate => update prev*/
	//Convolving with designed kernel will filter the signal
	//Offsetting should not affect the integration value as the value will be modified  based on that

	cue_data->t = sig_time - initial_time;
	cue_data->acc_fltrd = Convolve_rt(&kernel[0], KERNEL_LENGTH, sig_acc_input, cue_data->Input_Buff, &(cue_data->circ_buff_idx));
	cue_data->acc_fltrd_scaled = scale_factor * cue_data->acc_fltrd; //scale 
	cue_data->velocity = Intergration_Trapezoidal(cue_data->acc_fltrd_scaled, cue_data->acc_fltrd_scaled_prev, cue_data->velocity_prev, cue_data->t_prev, cue_data->t); //Integration
	cue_data->position = Intergration_Trapezoidal(cue_data->velocity, cue_data->velocity_prev, cue_data->position_prev, cue_data->t_prev, cue_data->t);

	//Updating the output values
	*(out_pos_) = cue_data->position + SP7_ZERO_POSE[data_index];//Offseting with respect to zeropose of SP7
	*(out_vel_) = cue_data->velocity;
	*(out_t_) = cue_data->t;
	//Updating the previous to current
	cue_data->t_prev = cue_data->t;
	cue_data->acc_fltrd_scaled_prev = cue_data->acc_fltrd_scaled;
	cue_data->velocity_prev = cue_data->velocity;
	cue_data->position_prev = cue_data->position;
}

void filtering::cue_rotational_channel(double sig_vel_input, double sig_time, double* kernel, double scale_factor, int data_index, filtering::CueDataVel* cue_data, double* out_pos_, double* out_vel_, double* out_t_)
{
	/* filter => scale => integrate => integrate => update prev*/
	//Convolving with designed kernel will filter the signal
	//Offsetting should not affect the integration value as the value will be modified  based on that

	cue_data->t = sig_time - initial_time;
	cue_data->velocity_fltr = Convolve_rt(&kernel[0], KERNEL_LENGTH, sig_vel_input, cue_data->Input_Buff, &(cue_data->circ_buff_idx));
	cue_data->velocity_fltr_scaled = scale_factor * cue_data->velocity_fltr; //scale 
	cue_data->position = Intergration_Trapezoidal(cue_data->velocity_fltr_scaled, cue_data->velocity_fltr_scaled_prev, cue_data->position_prev, cue_data->t_prev, cue_data->t);

	//Updating the output values
	*(out_pos_) = cue_data->position + SP7_ZERO_POSE[data_index];//Offsetting with respect to zero pose of SP7
	*(out_vel_) = cue_data->velocity_fltr_scaled;
	*(out_t_) = cue_data->t;
	//Updating the previous to current
	cue_data->t_prev = cue_data->t;
	cue_data->velocity_fltr_scaled_prev = cue_data->velocity_fltr_scaled;
	cue_data->position_prev = cue_data->position;
}

void filtering::cue_tilt_coordination_channel(double sig_acc_input, double sig_time, double* kernel, double scale_factor, int data_index, filtering::CueData* cue_data, double* out_ang_, double* out_t_)
{
	// LPF => Tilt cooridnation =>Integration => rate limit [x and y axis only -> not the up axis] 
	/*
	Max tilt angle	      = 5 deg
	Max tilt rate		  = 5 deg/s
	Max tilt acceleration = 8 deg/s^2
	Ref - 10.1177/0037549716675955*/

	double rate_limit_factor = 1;
	double g = 9.81;
	cue_data->t = sig_time - initial_time;
	cue_data->acc_fltrd = filtering::Convolve_rt(&kernel[0], KERNEL_LENGTH, sig_acc_input, cue_data->Input_Buff, &(cue_data->circ_buff_idx));
	cue_data->acc_fltrd_scaled = scale_factor * cue_data->acc_fltrd; //scale 
	// Tilt scaling factor
	double K = 5;//form factor [3-6]
	double theta_max = 5 * (M_PI / 180);
	double Acc_max = g * theta_max;
	//double y_tilt_ref = Acc_max * tanh(cue_data->acc_fltrd_scaled / K * Acc_max);
	// Tilt coordination --> Linerzation [Reid and Nahon]
	cue_data->velocity = cue_data->acc_fltrd_scaled / g;
	if (data_index == 0)
	{
		cue_data->velocity = -cue_data->velocity;
	}
	// Integration 
	cue_data->position = filtering::Intergration_Trapezoidal(cue_data->velocity, cue_data->velocity_prev, cue_data->position_prev, cue_data->t_prev, cue_data->t);
	//Rate Limit [1-5 deg]
	cue_data->position = rate_limit_factor * cue_data->position;
	//Updating the output values
	*(out_ang_) = cue_data->position;
	*(out_t_) = cue_data->t;
	//Updating the previous to current
	cue_data->t_prev = cue_data->t;
	cue_data->acc_fltrd_scaled_prev = cue_data->acc_fltrd_scaled;
	cue_data->velocity_prev = cue_data->velocity;
	cue_data->position_prev = cue_data->position;
}


#pragma region Helper functions RT

double filtering::Convolve_rt(double* h, int h_size, double x_in, double* x, int* circ_index)
{
	/*
	Real time convolution for one input signal and filter kernel
	Convolution with MAC(multiply and accumulate) and circulation buffering (shifting)
	executes only when run first and then treats as a global variable
	*/
	//circ_index = h_size - 1;

	x[*(circ_index)] = x_in;

	//Convolution
	double y_out = 0.0;
	for (int k = 0; k < h_size; k++)
		y_out += h[k] * x[(k + *(circ_index)) % h_size]; //MAC

	//Signal shift circularly through  array x in time-reversed order
	*(circ_index) += h_size - 1;
	*(circ_index) %= h_size;
	return y_out;
}

double filtering::Intergration_Trapezoidal(double input_curr, double input_prev, double output_prev, double t_prev, double t_curr)
{
	if (input_curr == 0 && input_prev == 0)
		return 0.0;
	return output_prev + ((t_curr - t_prev) * ((input_curr + input_prev) / 2)); // Trapezoidal integral
}
