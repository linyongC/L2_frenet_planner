#include "utils.h"

namespace utils {

bool LoadJsonFile(string scene_path, json* j) {
  ifstream f(scene_path);
  if (!f.is_open()) {
    cout << "Failed to open file: " << scene_path << endl;
    return false;
  }

  std::string s((std::istreambuf_iterator<char>(f)),
                std::istreambuf_iterator<char>());

  if (s.empty()) {
    cout << "File is empty: " << scene_path << endl;
    return false;
  }

  try {
    *j = json::parse(s);
  } catch (const std::exception& e) {
    cout << "Json parsing error for file: " << scene_path << endl;
    return false;
  }

  f.close();
  return true;
}

double wrap_angle(double angle) {
  while (angle > M_PI) {
    angle -= 2 * M_PI;
  }
  while (angle <= -M_PI) {
    angle += 2 * M_PI;
  }
  return angle;
}

double genGaussianNoise(double mean, double std_dev) {
  static std::random_device rd; // 提供随机数种子
  static std::mt19937 gen(rd()); // Mersenne Twister随机数生成器
  std::normal_distribution<double> d(mean, std_dev);
  return d(gen); // 生成符合该正态分布的随机数
}

}  // namespace utils