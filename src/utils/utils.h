#ifndef FRENET_OPTIMAL_TRAJECTORY_UTILS_H
#define FRENET_OPTIMAL_TRAJECTORY_UTILS_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <msgpack.hpp>
#include <nlohmann/json.hpp>
#include <random>
#include <tuple>
#include <vector>

using namespace std;
using json = nlohmann::json;

namespace utils {

inline double deg2rad(double deg) { return deg * M_PI / 180.0; }
inline double rad2deg(double rad) { return rad * 180.0 / M_PI; }
double wrap_angle(double angle);  // (-pi, pi]

typedef array<vector<double>, 2> WayPoints;  // [x_vec, y_vec]

struct Pose {
  double x;    // [m]
  double y;    // [m]
  double yaw;  // [rad]

  // 重载了加法减法运算符 +/-，使得两个 Pose 对象可以直接相加减。
  Pose operator+(const Pose& p) const {
    return Pose{x + p.x, y + p.y, yaw + p.yaw};
  }

  Pose operator-(const Pose& p) const {
    return Pose{x - p.x, y - p.y, yaw - p.yaw};
  }

  friend std::ostream& operator<<(std::ostream& os, const Pose& p) {
    os << " Pose(" << p.x << ", " << p.y << "," << utils::rad2deg(p.yaw)
       << "). ";
    return os;
  }

  MSGPACK_DEFINE(x, y, yaw);
};
struct Twist {
  double vx;
  double vy;
  double yaw_rate;

  // 重载了输出流运算符 <<，使得可以直接使用 std::cout 或其他输出流对象来打印 Twist 对象的内容。输出格式为：
  // Twist(vx值, vy值, yaw_rate的度数值).
  friend std::ostream& operator<<(std::ostream& os, const Twist& t) {
    os << " Twist(" << t.vx << ", " << t.vy << "," << utils::rad2deg(t.yaw_rate)
       << "). ";
    return os;
  }

  // 高效的二进制序列化
  // Twist 对象就可以方便地被转换为二进制格式进行存储或传输，并且可以轻松地从二进制格式恢复
  MSGPACK_DEFINE(vx, vy, yaw_rate);
};
struct Accel {
  double ax;
  double ay;
  double yaw_accel;

  MSGPACK_DEFINE(ax, ay, yaw_accel);
};

enum class CoordinateType {
  CARTESIAN = 0,  //
  FRENET = 1
};

enum LonMotionMode {
  UNDEFINED = 0,  //
  Following = 1,  //
  VelocityKeeping = 2
};

inline double norm(double x, double y) { return sqrt(pow(x, 2) + pow(y, 2)); }

inline void as_unit_vector(tuple<double, double>& vec) {
  double magnitude = norm(get<0>(vec), get<1>(vec));
  if (magnitude > 0) {
    get<0>(vec) = get<0>(vec) / magnitude;
    get<1>(vec) = get<1>(vec) / magnitude;
  }
}

inline double dot(const tuple<double, double>& vec1,
                  const tuple<double, double>& vec2) {
  return get<0>(vec1) * get<0>(vec2) + get<1>(vec1) * get<1>(vec2);
}

inline double cross_product(const tuple<double, double>& vec1,
                            const tuple<double, double>& vec2) {
  return get<0>(vec1) * get<1>(vec2) - get<1>(vec1) * get<0>(vec2);
}

bool LoadJsonFile(string scene_path, json* j);

double genGaussianNoise(double mean, double std_dev);

}  // namespace utils

#endif  // FRENET_OPTIMAL_TRAJECTORY_UTILS_H
