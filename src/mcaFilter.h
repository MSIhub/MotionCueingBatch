#ifndef MCA_FILTER_H
#define MCA_FILTER_H

#define _USE_MATH_DEFINES

#include <math.h>
#include <map>
#include "filtering.h"
#include <Eigen/Dense>
using namespace Eigen;

#define SOFTSATURATION
#define NUM_DATA 12


// Imported from common to be adjusted
#ifdef M_PI
static constexpr double pi = M_PI;
#else
#error "undefined M_PI"
#endif






using Vector6d = Eigen::Matrix<double, 6, 1>;

constexpr double pi_180{pi / 180.};
constexpr double pi_180_inv{180. / pi};
constexpr double deg2rad(double angleInDegree)
{
  return angleInDegree * pi_180;
}
constexpr double rad2deg(double angleRad) { return angleRad * pi_180_inv; }

//

inline double linear_scaling(double x, double f) { return x / f; };

template <typename T> T clamp(const T &n, const T &lower, const T &upper)
{
  return std::max(lower, std::min(n, upper));
}

class McaFilter {
private:
    std::map<std::string, double> paramMap;

    double t{};     // actual sample time
    double tprev{}; // previous sample time
    double a_g{};
    double b_g{};

    Vector3d aprev{ 0, 0, 0 };
    Vector3d a{ 0, 0, 0 };
    Vector3d pos{ 0, 0, 0.005 };
    Vector3d vel{ 0, 0, 0 };
    Vector3d theta_dot_h{ 0, 0, 0 };

	double high_pass_kernel[KERNEL_LENGTH];
	double low_pass_kernel[KERNEL_LENGTH];
    bool init_run;
    

    //CueData
	filtering::CueData* c_ax;
	filtering::CueData* c_ay;
	filtering::CueData* c_az;
	filtering::CueData* c_tcx;
	filtering::CueData* c_tcy;
	filtering::CueDataVel* c_vroll;
	filtering::CueDataVel* c_vpitch;
	filtering::CueDataVel* c_vyaw;

    //logging
    bool log_data;
	std::string delimiter = "\t"; //tab limited text file with 8 point precision
    std::fstream log_fptr;
    std::string log_filename = "log_data_batch_";

public:
    void loadParams();
    void printParams();
    McaFilter():init_run(true), log_data(true)
    {
        loadParams();
		// Init CueData
		c_ax = new filtering::CueData{};
		c_ay = new filtering::CueData{};
		c_az = new filtering::CueData{};
		c_tcx = new filtering::CueData{};
		c_tcy = new filtering::CueData{};
		c_vroll = new filtering::CueDataVel{};
		c_vpitch = new filtering::CueDataVel{};
		c_vyaw = new filtering::CueDataVel{};
        //
        filtering::initial_time = 0.0;

    #pragma region Log data
		log_filename += std::to_string(std::time(nullptr));
		log_filename += ".log";
		log_fptr.open(log_filename, std::fstream::out | std::fstream::app);
    #pragma endregion

    }
    McaFilter(McaFilter &c) = delete; // deleting copy constructor
		//Hard saturation
	bool saturate_min(double& x, double threshold, double& x_dot);
	bool saturate_max(double& x, double threshold, double& x_dot);

    Matrix3d angular_to_gimbal_matrix(double a1, double a2, double a3);
    Matrix3d get_R_process(double roll, double pitch, double yaw);
    Matrix3d get_T_process(double roll, double pitch, double yaw);
    Vector3d tilt_cord(Vector3d f_L);

    void calculateKernels();
    void filtering(float data[NUM_DATA]);
    void getData(double data[NUM_DATA]);
    void reset();
    void logFilteredData(double data[NUM_DATA], double timestamp);

    ~McaFilter() { log_fptr.close(); }
};




namespace
{
    // static constexpr double fs {47.092}; // Should be 60

    static constexpr double g_acc{ 9.8 }; //< Gravity acceleration
    // scaling factor for linear scaling
    static double scale_x{ 50. };
    static double scale_y{ 50. };
    static double scale_z{ 100. };
    static double scale_wx{ 2. };
    static double scale_wy{ 2. };
    static double scale_wz{ 2. };
    // Saturation acceleration limits

    static double Sat_ax_Min{ -15. };
    static double Sat_ay_Min{ -15. };
    static double Sat_az_Min{ -15. };
    static double Sat_ax_Max{ 15. };
    static double Sat_ay_Max{ 15. };
    static double Sat_az_Max{ 15. };

    // Saturation angular velocity limits
    static double Sat_wx_Min{ -0.4 };
    static double Sat_wy_Min{ -0.4 };
    static double Sat_wz_Min{ -0.4 };
    static double Sat_wx_Max{ 0.4 };
    static double Sat_wy_Max{ 0.4 };
    static double Sat_wz_Max{ 0.4 };

    // Complex way to do zeroPose = [0 0 0.401 0 0 0];
    static const Vector6d
        zeroPose((Vector6d() << 0.0, 0.0, 0.401, 0.0, 0.0, 0.0).finished());

    // Boundary translation limits

    static double xMinLimit{ -0.05 };
    static double yMinLimit{ -0.05 };
    static double zMinLimit{ 0.396 };
    static double xMaxLimit{ 0.05 };
    static double yMaxLimit{ 0.05 };
    static double zMaxLimit{ 0.451 };

    // Boundary rotational limits
    static double thetaxMinLimit{ deg2rad(-5.) };
    static double thetayMinLimit{ deg2rad(-5.) };
    static double thetazMinLimit{ deg2rad(-180.) };
    static double thetaxMaxLimit{ deg2rad(5.) };
    static double thetayMaxLimit{ deg2rad(5.) };
    static double thetazMaxLimit{ deg2rad(180.) };

    // Mid point of Z axis
    static double offset{ (zMaxLimit + zMinLimit) / 2. };

} // namespace McaParam

#endif
