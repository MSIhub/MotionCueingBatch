#include "mcaFilter.h"



/**
 *
 *
 */
Matrix3d McaFilter::angular_to_gimbal_matrix(double a1, double a2, double a3)
{
  Matrix3d gimbal;
  gimbal << 1, (sin(a1) * sin(a2)) / cos(a2), -(cos(a1) * sin(a2)) / cos(a2), 0,
      cos(a1), sin(a1), 0, -sin(a1) / cos(a2), cos(a1) / cos(a2);

  return gimbal;
}

/**
 *
 */
Matrix3d McaFilter::get_R_process(double roll, double pitch, double yaw)
{
  Matrix3d R_p;

  R_p << cos(pitch) * cos(roll),
      sin(yaw) * sin(roll) + cos(yaw) * cos(roll) * sin(pitch),
      cos(roll) * sin(yaw) * sin(pitch) - cos(yaw) * sin(roll), -sin(pitch),
      cos(yaw) * cos(pitch), cos(pitch) * sin(yaw), cos(pitch) * sin(roll),
      cos(yaw) * sin(pitch) * sin(roll) - cos(roll) * sin(yaw),
      cos(yaw) * cos(roll) + sin(yaw) * sin(pitch) * sin(roll);

  return R_p;
}

/**
 *
 */
Matrix3d McaFilter::get_T_process(double roll, double pitch, double yaw)
{
  Matrix3d T_p;
  T_p << 1, (cos(yaw) * sin(pitch)) / cos(pitch),
      (sin(yaw) * sin(pitch)) / cos(pitch), 0, -sin(yaw), cos(yaw), 0,
      cos(yaw) / cos(pitch), sin(yaw) / cos(pitch);

  return T_p;
}

/**
 * \param   f_L is low pass filtered trasnlation data
 *  \return is the acceleration of platform in platform frame
 */
Vector3d McaFilter::tilt_cord(Vector3d f_L)
{
  return {f_L[1] / g_acc, -f_L[0] /g_acc, 0.0};
}


/**
 * \param   f_L is low pass filtered trasnlation data
 *  \return is the acceleration of platform in platform frame
 */

void McaFilter::loadParams()
{
	std::string filename = "src/param.yaml";
	std::ifstream param;


	param.open(filename);
	if (!param.is_open())
	{
		std::cout << "file: " << filename << " could not open!" << std::endl;
		return;
	}

	while (param)
	{
		std::string key;
		double value;
		std::getline(param, key, ':');
		param >> value;
		param.get(); // catch empty line
		if (!param)
			return;
		paramMap[key] = value; //paramMap.insert(std::pair<std::string, float>(key, value));
	}
	param.close();
	return;
}

void McaFilter::printParams()
{
	std::map<std::string, double>::iterator itr;
	for (itr = paramMap.begin(); itr != paramMap.end(); ++itr)
	{
		std::cout << itr->first << ": " << itr->second << std::endl;
	}
}

void McaFilter::filtering(float data[NUM_DATA])
{
	t = data[0]; // TODO We need some origin to start and convert t to relative time
	double f_ggx = data[1]; // Specific force along x - it really is ax
	double f_ggy = data[2]; // Specific force along y - it really is ay
	double f_ggz = data[3]; // Specific force along z - it really is az

	double roll = data[4] * pi_180;
	double pitch = data[5] * pi_180;
	double yaw = data[6] * pi_180;

	double vroll = data[7] * pi_180;
	double vpitch = data[8] * pi_180;
	double vyaw = data[9] * pi_180;

	/*Insert filter logic here*/
	if (init_run)
	{
		//cal kernel
		calculateKernels();
		init_run = false;
	}


	//------------- Cueing function -------------------------//
	filtering::SP7Pose* pose = new filtering::SP7Pose{ 0,0,0,0,0,0 };// initialization to avoid sending garbage value
	filtering::SP7Vel* velocity = new filtering::SP7Vel{ 0,0,0,0,0,0 };
	double timestamp = 0.0;

	//Translational channel
	cue_translational_channel(f_ggx, t, &high_pass_kernel[0], paramMap["k_ax"], 0, c_ax, &(pose->x), &(velocity->vx), &timestamp);
	cue_translational_channel(f_ggy, t, &high_pass_kernel[0], paramMap["k_ay"], 1, c_ay, &(pose->y), &(velocity->vy), &timestamp);
	cue_translational_channel(f_ggz, t, &high_pass_kernel[0], paramMap["k_az"], 2, c_az, &(pose->z), &(velocity->vz), &timestamp);


	// Tilt coordination channel
	double tilt_x = 0.0;
	cue_tilt_coordination_channel(f_ggx, t, &low_pass_kernel[0], paramMap["k_ax"], 0, c_tcx, &tilt_x, &timestamp);
	double tilt_y = 0.0;
	cue_tilt_coordination_channel(f_ggy, t, &low_pass_kernel[0], paramMap["k_ay"], 1, c_tcy, &tilt_y, &timestamp);

	// Rotational channel
	cue_rotational_channel(vroll, t, &high_pass_kernel[0], paramMap["k_vroll"], 3, c_vroll, &(pose->roll), &(velocity->vroll), &timestamp);
	pose->roll += tilt_x; // adding tilt effect
	cue_rotational_channel(vpitch, t, &high_pass_kernel[0], paramMap["k_vpitch"], 4, c_vpitch, &(pose->pitch), &(velocity->vpitch), &timestamp);
	pose->pitch += tilt_y;// adding tilt effect
	cue_rotational_channel(vyaw, t, &high_pass_kernel[0], paramMap["k_vyaw"], 5, c_vyaw, &(pose->yaw), &(velocity->vyaw), &timestamp);

	//Pack the output
	pos[0] = pose->x;
	pos[1] = pose->y;
	pos[2] = pose->z;
	a[0] = pose->roll;
	a[1] = pose->pitch;
	a[2] = pose->yaw;

	vel[0] = velocity->vx;
	vel[1] = velocity->vy;
	vel[2] = velocity->vz;
	theta_dot_h[0] = velocity->vroll;
	theta_dot_h[1] = velocity->vpitch;
	theta_dot_h[2] = velocity->vyaw;

	delete pose, velocity;
	

}

void McaFilter::getData(double data[NUM_DATA])
{
	data[0] = pos[0];
	data[1] = pos[1];
	data[2] = pos[2];
	data[3] = a[0];
	data[4] = a[1];
	data[5] = a[2];

	data[6] = vel[0];
	data[7] = vel[1];
	data[8] = vel[2];
	data[9] = theta_dot_h[0];
	data[10] = theta_dot_h[1];
	data[11] = theta_dot_h[2];
}

void McaFilter::logFilteredData(double data[NUM_DATA], double timestamp)
{
	
		//Preparing the data stream
		std::stringstream ss;
		ss.precision(8);// max to micro meter
		ss << std::fixed << timestamp << delimiter;
		for (int i = 0; i < NUM_DATA; i++)
			ss << data[i] << delimiter;
		
		if (log_fptr.is_open())
		{
			log_fptr << "\n";
			log_fptr << ss.str();
		}
}

void McaFilter::reset()
{
	t = 0.;
	tprev = 0.;
	a = Vector3d::Zero();
	aprev = Vector3d::Zero();
	pos = Vector3d::Zero();
	vel = Vector3d::Zero();
	theta_dot_h = Vector3d::Zero();

}

void McaFilter::calculateKernels()
{
	//Kernel calculation beforehand
	filtering::calc_kernel_high_pass(KERNEL_LENGTH, paramMap["hp_ax"], &high_pass_kernel[0], log_data);
	filtering::calc_kernel_low_pass(KERNEL_LENGTH, paramMap["lp_ax"], &low_pass_kernel[0], log_data);
}