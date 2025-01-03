#include <gflags/gflags.h>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <utility>

#include "FrenetOptimalTrajectory/FrenetOptimalTrajectory.h"
#include "FrenetOptimalTrajectory/py_cpp_struct.h"
#include "FrenetPath.h"
#include "Lane/Lane.h"
#include "py_cpp_struct.h"
#include "utils/coordinate_utils.h"
#include "utils/data_log.h"
#include "utils/debug.h"
#include "utils/utils.h"

using namespace std;
using namespace utils;
using json = nlohmann::json;

const std::string REPO_DIR = "/home/sheng/Projects/l2_frenet_planner/";

DEFINE_string(scene_path, REPO_DIR + "config/scenes/one_lane_slow_down.json",
              "Path to scene config file");
DEFINE_string(hyper_path, REPO_DIR + "config/hyperparameters.json",
              "Path to hyperparameter config file");
DEFINE_bool(store_data, false, "turn on flag to store running data.");
DEFINE_string(data_path, REPO_DIR + "build/data.bin",
              "Path to store running data.");
DEFINE_bool(local_planning, true, "turn on flag to enable local planning.");

double get_duration_ms(
    std::chrono::time_point<std::chrono::high_resolution_clock> end,
    std::chrono::time_point<std::chrono::high_resolution_clock> start) {
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  return static_cast<double>(duration.count()) / 1000.0;
}

double round_to_tenth(double x) { return std::round(x * 10.0) / 10.0; }

bool InitFrenetHyperParameters() {
  json hyper_j;
  if (!LoadJsonFile(FLAGS_hyper_path, &hyper_j)) {
    return false;
  }
  FrenetHyperparameters::getInstance().Init(hyper_j);
  return true;
}

void UpdateFrenetCoordinates(const Car& car, const utils::WayPoints& wp,
                             FrenetInitialConditions* fot_ic) {
  fot_ic->yaw_c = car.getPose().yaw;

  Car car_f;
  utils::ToFrenet(car, wp, &car_f);  // 将car_c投影到local_wp上面获取car_f信息
  fot_ic->s = car_f.getPose().x;
  fot_ic->s_d = car_f.getTwist().vx;
  fot_ic->s_dd = car_f.getAccel().ax;
  fot_ic->d = car_f.getPose().y;
  fot_ic->d_d = car_f.getTwist().vy;
  fot_ic->d_dd = car_f.getAccel().ay;
}

// void InitFrenetInitialConditions(const Car& car, const json& scene_j,
//                                  const utils::WayPoints& wp,
//                                  FrenetInitialConditions* fot_ic) {
//   UpdateFrenetCoordinates(car, wp, fot_ic);
// }

void InitWaypoints(const json& scene_j, WayPoints* wp) {
  for (const auto& point : scene_j["wp"]) {
    (*wp)[0].push_back(point[0]);
    (*wp)[1].push_back(point[1]);
  }
}

void InitObstacles(const json& scene_j, const vector<Lane>& lanes,
                   vector<Obstacle>* const obstacles) {
  utils::WayPoints ref_wp;
  InitWaypoints(scene_j, &ref_wp);

  const auto& fot_hp = FrenetHyperparameters::getConstInstance();
  // const double dt = fot_hp.dt;
  vector<Obstacle> obstacles_f;
  for (const auto& ob_j : scene_j["obs"]) {
    Pose ob_pose = {ob_j["pose"][0], ob_j["pose"][1], ob_j["pose"][2]};
    Obstacle ob(ob_pose, ob_j["length"], ob_j["width"],
                fot_hp.obstacle_clearance);
    if (ob_j.contains("speed_profile")) {
      std::map<double, double> spd_profile;
      for (const auto& spd_j : ob_j["speed_profile"]) {
        spd_profile[spd_j[0]] = spd_j[1];
      }
      ob.setSpeedLookupTable(spd_profile);
      // ob.predictPoses(0.0, fot_hp.maxt, dt);
      ob.setTwist({spd_profile.begin()->second, 0.0, 0.0});
    }
    obstacles_f.push_back(std::move(ob));
  }

  // convert to Cartesian coordinates
  // note: the obstacle initial pose defined in scene in Frenet frame is w.r.t.
  // ref_wp
  for (const auto& ob_f : obstacles_f) {
    std::unique_ptr<Obstacle> ob_c = nullptr;
    utils::ToCartesian(ob_f, ref_wp, ob_c);
    obstacles->push_back(std::move(*ob_c));
  }

  // set lane id for each obstacle as well.
  for (auto& ob : *obstacles) {
    ob.clearLaneIds();
    for (std::size_t i = 0; i < lanes.size(); ++i) {
      if (point_in_lane(lanes[i], ob.getPose().x, ob.getPose().y)) {
        ob.addLaneId(i);
      }
    }
    auto& lane_ids = ob.getLaneIds();
    vector<int> vec_lane_ids(lane_ids.begin(), lane_ids.end());
    cout << "ob lane_ids: " << utils::vector_to_str(vec_lane_ids) << endl;
  }
}

void UpdateNextPlanningStateLocal(const FrenetPath* best_frenet_path_local,
                                  const WayPoints& wp_local,
                                  Car* next_planning_state_local) {
  const auto& fot_hp = FrenetHyperparameters::getConstInstance();

  // update ego car to next state
  double next_s = best_frenet_path_local->s[1];
  double next_d = best_frenet_path_local->d[1];
  double next_s_d = best_frenet_path_local->s_d[1];
  double next_d_d = best_frenet_path_local->d_d[1];
  double next_s_dd = best_frenet_path_local->s_dd[1];
  double next_d_dd = best_frenet_path_local->d_dd[1];
  double next_yaw_f = utils::wrap_angle(std::atan2(next_d_d, next_s_d));
  double next_yaw_d_f = (next_s_d * next_d_dd - next_d_d * next_s_dd) /
                        (next_s_d * next_s_d + next_d_d * next_d_d);
  Pose pose_c;
  Twist twist_c;
  Accel accel_c;
  ToCartesian({next_s, next_d, next_yaw_f}, {next_s_d, next_d_d, next_yaw_d_f},
              {next_s_dd, next_d_dd, 0.0}, wp_local, &pose_c, &twist_c,
              &accel_c);

  // PATCH: The yaw and yaw rate will be unexpected if path too short
  if (std::fabs(best_frenet_path_local->s.back() -
                best_frenet_path_local->s.front()) < 10.0) {
    pose_c.yaw = utils::wrap_angle(best_frenet_path_local->yaw[1]);
    twist_c.yaw_rate = utils::wrap_angle(best_frenet_path_local->yaw[2] -
                                         best_frenet_path_local->yaw[0]) /
                       2.0 / fot_hp.dt;
    twist_c.vx = (best_frenet_path_local->s[2] - best_frenet_path_local->s[0]) /
                 2.0 / fot_hp.dt;
    twist_c.vy = 0.0;
    accel_c.ax =
        (best_frenet_path_local->s_d[2] - best_frenet_path_local->s_d[0]) /
        2.0 / fot_hp.dt;
    accel_c.ay = 0.0;
    accel_c.yaw_accel = utils::wrap_angle(best_frenet_path_local->yaw[2] -
                                          2.0 * best_frenet_path_local->yaw[1] +
                                          best_frenet_path_local->yaw[0]) /
                        fot_hp.dt / fot_hp.dt;
  }

  next_planning_state_local->setPose(pose_c);
  next_planning_state_local->setTwist(twist_c);
  next_planning_state_local->setAccel(accel_c);
  next_planning_state_local->setTargetLaneId(best_frenet_path_local->lane_id);
}

Car InitEgoCar(const json& scene_j, const vector<Lane>& lanes) {
  Pose ego_car_pose = {scene_j["pose"][0], scene_j["pose"][1],
                       scene_j["pose"][2]};
  Twist ego_car_twist = {scene_j["vel"][0], scene_j["vel"][1],
                         scene_j["vel"][2]};
  Car ego_car(ego_car_pose, ego_car_twist, {0, 0, 0});

  // initialize lane_id for ego car
  int lane_id = -1;
  for (std::size_t i = 0; i < lanes.size(); ++i) {
    if (point_in_lane(lanes[i], ego_car_pose.x, ego_car_pose.y)) {
      lane_id = i; // 判断ego所在车道
      break;
    }
  }
  ego_car.setCurLaneId(lane_id);
  ego_car.setTargetLaneId(lane_id);

  // align ego heading w.r.t lane heading 将车辆heading设置为lane heading
  if (lane_id != -1) {
    const auto& wp = lanes[lane_id].GetWayPoints();
    double lane_heading =
        utils::wrap_angle(std::atan2(wp[1][1] - wp[1][0], wp[0][1] - wp[0][0]));
    ego_car.setPose({ego_car_pose.x, ego_car_pose.y, lane_heading});
  }

  return ego_car;
}

void InitLanes(const json& scene_j, vector<Lane>* lanes) {
  utils::WayPoints ref_wp;
  InitWaypoints(scene_j, &ref_wp); // 获取scene_j中的waypoints并存入ref_wp

  const double lane_width = scene_j["lane_width"];
  const int num_lanes_left = scene_j["num_lanes_left"];
  const int num_lanes_right = scene_j["num_lanes_right"];
  const int total_lanes = num_lanes_left + num_lanes_right + 1;

  for (int i = -num_lanes_left; i <= num_lanes_right; ++i) {
    utils::WayPoints lane_wp;
    if (i == 0) {
      lane_wp = ref_wp; // 当前ego车道
    } else {
      utils::ShiftWaypoints(ref_wp, -i * lane_width, // 先使用三次样条插值，在沿yaw垂直方向shift，得到lane_wp
                            &lane_wp);  // positive shift to the left
    }
    int lane_id = i + num_lanes_left; // 从左往右依次为0，1，2，3，，，
    Lane lane(lane_id, lane_wp, lane_width); // 实例化lane并添加到vector中
    lanes->push_back(lane);
  }

  // lane associations
  for (std::size_t i = 0; i < lanes->size(); ++i) {
    if (i > 0) {
      (*lanes)[i].SetLeftLane(&(*lanes)[i - 1]);
    }
    if (i < lanes->size() - 1) {
      (*lanes)[i].SetRightLane(&(*lanes)[i + 1]); // 此处应该是SetRightLane
    }
  }

  // calculate lane boundaries
  utils::WayPoints left_lane_bound;
  utils::WayPoints right_lane_bound;
  for (int i = -num_lanes_left; i <= num_lanes_right; ++i) { // 计算每条lane的边界
    if (i == -num_lanes_left) { 
      utils::ShiftWaypoints(ref_wp, (-i + 0.5) * lane_width, &left_lane_bound);
    }
    utils::ShiftWaypoints(ref_wp, (-i - 0.5) * lane_width, &right_lane_bound);

    Lane& lane = (*lanes)[i + num_lanes_left];
    lane.SetLeftLaneBoundary(left_lane_bound);
    lane.SetRightLaneBoundary(right_lane_bound);

    left_lane_bound = right_lane_bound;
  }
}

/**
 * @brief Estimate the change of pose based on speed and yaw rate
 *
 * @param speed_meas
 * @param yaw_rate_meas
 * @return Pose
 */
Pose EstimateChangeOfPose(double speed_meas, double yaw_rate_meas) {
  double delta_x = 0.0;
  double delta_y = 0.0;
  double delta_yaw = 0.0;

  const auto& fot_hp = FrenetHyperparameters::getConstInstance();
  constexpr double kSimuStep = 0.002;  // [s]
  const int kNumSteps = static_cast<int>(fot_hp.dt / kSimuStep);
  for (int i = 0; i < kNumSteps; ++i) {
    delta_x += speed_meas * std::cos(delta_yaw) * kSimuStep;
    delta_y += speed_meas * std::sin(delta_yaw) * kSimuStep;
    delta_yaw += yaw_rate_meas * kSimuStep;
  }
  return {delta_x, delta_y, delta_yaw};
}

/**
 * @brief Estimate the change of pose by add-on noise
 * TEMPORARY FIX.
 * @return Pose
 */
Pose EstimateChangeOfPoseV2(const Car& next_planning_state_local,
                            double* speed_meas, double* yaw_rate_meas) {
  // TODO: tmp fix until control and vehicle model added to better related
  // speed/yaw_rate sensor measurements to the estimated state change
  double delta_x_noise = 0.0;
  double delta_y_noise = 0.0;
  double delta_yaw_noise = 0.0;

  const auto& fot_hp = FrenetHyperparameters::getConstInstance();
  double speed_noise = utils::genGaussianNoise(fot_hp.sensor_speed_offset,
                                               fot_hp.sensor_speed_noise_std);
  double yaw_rate_noise = utils::genGaussianNoise(
      fot_hp.sensor_yaw_rate_offset, fot_hp.sensor_yaw_rate_noise_std);

  constexpr double kSimuStep = 0.002;  // [s]
  const int kNumSteps = static_cast<int>(fot_hp.dt / kSimuStep); // 0.1/0.002 = 50
  for (int i = 0; i < kNumSteps; ++i) {
    delta_x_noise += speed_noise * std::cos(delta_yaw_noise) * kSimuStep;
    delta_y_noise +=
        (next_planning_state_local.getTwist().vx * std::sin(delta_yaw_noise) +
         speed_noise * std::sin(delta_yaw_noise)) *
        kSimuStep;
    delta_yaw_noise += yaw_rate_noise * kSimuStep;
  }

  // or directly
  // double delta_x_noise = utils::genGaussianNoise(0.1, 0.1);
  // double delta_y_noise = utils::genGaussianNoise(0.04, 0.04);

  double delta_x = next_planning_state_local.getPose().x + delta_x_noise;
  double delta_y = next_planning_state_local.getPose().y + delta_y_noise;
  double delta_yaw = next_planning_state_local.getPose().yaw + delta_yaw_noise;

  *speed_meas = next_planning_state_local.getTwist().vx + speed_noise;
  *yaw_rate_meas =
      next_planning_state_local.getTwist().yaw_rate + yaw_rate_noise;

  return {delta_x, delta_y, delta_yaw};
}

// void UpdateSensorMeasurements(const Car& ego_car, const Car& next_ego_car,
//                               double* speed_meas, double* yaw_rate_meas) {
//   const auto& fot_hp = FrenetHyperparameters::getConstInstance();
//   const double TimeStep = fot_hp.dt;
//   double position_change =
//       utils::norm(next_ego_car.getPose().x - ego_car.getPose().x,
//                   next_ego_car.getPose().y - ego_car.getPose().y);
//   double yaw_change =
//       utils::wrap_angle(next_ego_car.getPose().yaw - ego_car.getPose().yaw);

//   *speed_meas = position_change / TimeStep;
//   *yaw_rate_meas = yaw_change / TimeStep;

//   // add offset and noise
//   *speed_meas += utils::genGaussianNoise(
//       fot_hp.sensor_speed_offset, fot_hp.sensor_speed_noise_std);  // [m/s]
//   *yaw_rate_meas +=
//       utils::genGaussianNoise(fot_hp.sensor_yaw_rate_offset,
//                               fot_hp.sensor_yaw_rate_noise_std);  // [rad/s]
// }

int main(int argc, char** argv) {
  // 调用gflags库来解析命令行参数。
  // 这里true的参数表示在解析后将命令行参数从argv中移除，使得后续代码不会看到这些已经处理的参数。
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if (!InitFrenetHyperParameters()) { // FLAGS_hyper_path来自命令行参数，初始化FrenetHyperParameters
    return 1;
  }
  json scene_j;
  if (!LoadJsonFile(FLAGS_scene_path, &scene_j)) { // FLAGS_scene_path来自命令行参数的scene_path
    cout << "Fail to load simulation scene: " << FLAGS_scene_path << endl;
    return 1;
  }

  vector<Lane> lanes; // 定义类型是Lane类的vector，未初始化
  InitLanes(scene_j, &lanes); // 初始化车道：三次样条插值、平移或取多个车道中心线、计算车道边界等

  Car ego_car = InitEgoCar(scene_j, lanes); // 初始化ego

  vector<Obstacle> obstacles;
  InitObstacles(scene_j, lanes, &obstacles); // 初始化bostacles，增加laneid，转化为cartesian坐标

  double target_speed = scene_j["target_speed"];

  const auto& fot_hp = FrenetHyperparameters::getConstInstance();
  const double TimeStep = fot_hp.dt; // 0.1s
  int sim_loop = 200;
  double total_runtime = 0.0;  // [ms]
  double timestamp = 0.0;      // [s], simulation timestamp
  int i = 0;
  std::vector<DataFrame> data_frames; // 使用msgpack::pack可将data_frames序列化存储
  bool reach_goal = false;
  std::vector<FrenetPath>
      best_frenet_paths_local;  // store the best frenet path for each lane
  std::vector<std::vector<FrenetPath>> frenet_paths_local_all;

  Car planning_init_point_local = ego_car; // init规划起点为当前位置
  if (FLAGS_local_planning) {
    planning_init_point_local.setPose({0, 0, 0}); // 使能local_planning时将当前位置置为0
  }

  double yaw_rate_meas = 0.0;
  double speed_meas = 0.0;
  Pose pose_change_est = {0.0, 0.0, 0.0};
  Car planning_init_point_wrt_last_frame; // 基于last frame的规划起点

  unordered_map<int, WayPoints> wp_lanes_local; // 哈希表（无序映射），它的键是int类型，值是WayPoints类型
  std::vector<Obstacle> obstacles_local;

  for (; i < sim_loop; ++i) {
    auto start = std::chrono::high_resolution_clock::now(); //获取当前时间点的时间戳

    // Loop each lane here and may initialize fot_ic for each lane
    best_frenet_paths_local.clear();
    frenet_paths_local_all.clear();
    frenet_paths_local_all.resize(lanes.size());
    for (const auto& lane : lanes) {
      WayPoints wp = lane.GetWayPoints();

      // break if near goal
      if (utils::norm(ego_car.getPose().x - wp[0].back(),
                      ego_car.getPose().y - wp[1].back()) < 10.0) {
        reach_goal = true;
        break;
      }

      // prediction on obstacles
      for (auto& ob : obstacles) { // 通过auto&声明的ob是对vector中元素的引用，更改会影响原obstacles
        std::unique_ptr<Obstacle> ob_f = nullptr;
        int ob_laneid = 0;
        if (!ob.getLaneIds().empty()) {
          ob_laneid = *(ob.getLaneIds().begin());
        }
        utils::ToFrenet(ob, lanes[ob_laneid].GetWayPoints(), ob_f); // 障碍物信息转换到障碍物所在lane的frenet坐标系
        ob_f->predictPoses(timestamp, fot_hp.maxt, TimeStep); // 预测轨迹（5s），匀速直线运动
        std::unique_ptr<Obstacle> ob_c = nullptr;
        utils::ToCartesian(*ob_f, lanes[ob_laneid].GetWayPoints(), ob_c); // 障碍物信息转换到cartesian坐标系
        ob = std::move(*ob_c);
      }

      // convert obstacles to local coordinate w.r.t. ego car
      obstacles_local.clear();
      WayPoints wp_local;
      if (FLAGS_local_planning) {
        ToLocal(obstacles, ego_car.getPose(), &obstacles_local);// 障碍物转换到车辆坐标系下
        ToLocal(wp, ego_car.getPose(), &wp_local);// wp转换到车辆坐标系下
      } else {
        obstacles_local = obstacles;
        wp_local = wp;
      }
      wp_lanes_local[lane.GetLaneId()] = wp_local;

      // run frenet optimal trajectory 简称fot
      FrenetInitialConditions fot_ic(wp_local, obstacles_local,
                                     lane.GetLaneId());
      fot_ic.target_speed = target_speed;
      fot_ic.lane_width = lane.GetLaneWidth();

      // update Frenet coordinate of ego car
      UpdateFrenetCoordinates(planning_init_point_local, wp_local, &fot_ic); // 获取fot_ic的s/d等信息，将car_c投影到local_wp上面获取car_f信息

      // local planning w.r.t. ego car，fot_ic是初始化的fot，涵盖投影点的frenet信息
      FrenetOptimalTrajectory fot = FrenetOptimalTrajectory(fot_ic, fot_hp); // 横纵向轨迹规划，求代价
      FrenetPath* best_frenet_path_per_lane = fot.getBestPath();
      if (!best_frenet_path_per_lane || best_frenet_path_per_lane->x.empty()) {
        cerr << "Fail to find a feasible path at timestamp: " << timestamp
             << ", at lane: " << lane.GetLaneId() << endl;
      } else {
        // update cost for each frenet path based on lane
        best_frenet_path_per_lane->lane_id = lane.GetLaneId();
        best_frenet_path_per_lane->c_lane_change =            // 对当前lane的最有轨迹增加变道代价，变道时非目标lane的轨迹的代价增加
            std::abs(best_frenet_path_per_lane->lane_id -
                     planning_init_point_local.getTargetLaneId());
        best_frenet_path_per_lane->cf +=
            fot_hp.klane * best_frenet_path_per_lane->c_lane_change;

        best_frenet_paths_local.push_back(*best_frenet_path_per_lane);
      }

      // store frenet paths at each d for each lane
      map<double, FrenetPath*> d_to_best_path_local; // 键是double类型的值，而值是指向FrenetPath对象的指针
      for (auto fp : fot.getFrenetPaths()) {
        double d = round_to_tenth(fp->d.back()); // 终端的d保留一位小数
        if (d_to_best_path_local.count(d) == 0) { // 没有和终端d对应的映射
          d_to_best_path_local[d] = fp;
        } else {
          if (fp->cf < d_to_best_path_local[d]->cf) { // 如果有，找到最小代价的轨迹，存到终端d处
            d_to_best_path_local[d] = fp;
          }
        }
      }

      for (const auto& d_fp : d_to_best_path_local) {
        d_fp.second->lane_id = lane.GetLaneId(); // 存储当前车道每一个终点d对应的最有路径
        frenet_paths_local_all[lane.GetLaneId()].push_back(*(d_fp.second));
      }
    }

    if (reach_goal) {
      break;
    }

    if (best_frenet_paths_local.empty()) {
      cerr << "Fail to find a feasible path at timestamp: " << timestamp
           << " for all lanes. Terminate!" << endl;
      break;
    }

    // choose from best trajectory along each lane based on cost
    FrenetPath* best_frenet_path_local = nullptr;
    for (auto& fp : best_frenet_paths_local) {
      if (!best_frenet_path_local) {
        best_frenet_path_local = &fp;
      } else {
        if (fp.cf < best_frenet_path_local->cf) {
          best_frenet_path_local = &fp;
        }
      }
    }

    auto plan_end = std::chrono::high_resolution_clock::now();
    double plan_duration = get_duration_ms(plan_end, start);

    // convert (x,y,yaw) of paths to global for debug purposes
    std::vector<FrenetPath> best_frenet_paths_global = best_frenet_paths_local;
    if (FLAGS_local_planning) {
      for (auto& fp : best_frenet_paths_global) {   // 将最优路径转换到global坐标系下，用于debug
        std::size_t fp_size = fp.x.size();
        for (std::size_t i = 0; i < fp_size; ++i) {
          Pose pose_l({fp.x[i], fp.y[i], fp.yaw[i]});
          Pose pose_g;
          ToGlobal(pose_l, ego_car.getPose(), &pose_g);
          fp.x[i] = pose_g.x;
          fp.y[i] = pose_g.y;
          fp.yaw[i] = pose_g.yaw;
        }
      }
    }

    FrenetPath* best_frenet_path_g = nullptr;
    for (auto& fp : best_frenet_paths_global) {
      if (!best_frenet_path_g) {
        best_frenet_path_g = &fp;
      } else {
        if (fp.cf < best_frenet_path_g->cf) {
          best_frenet_path_g = &fp;
        }
      }
    }

    // save current frame data
    if (FLAGS_store_data) {
      DataFrame df;
      df.timestamp = timestamp;
      df.ego_car = ego_car;
      df.best_frenet_path = *best_frenet_path_g;
      df.lanes = lanes;
      df.obstacles = obstacles;
      df.frenet_paths = best_frenet_paths_global;
      df.obstacles_local = obstacles_local;
      df.wp_lanes_local = wp_lanes_local;
      df.planning_init_point_local = planning_init_point_local;
      df.best_frenet_path_local = *best_frenet_path_local;
      df.frenet_paths_local = best_frenet_paths_local;
      df.frenet_paths_local_all = frenet_paths_local_all;
      df.speed_meas = speed_meas;
      df.yaw_rate_meas = yaw_rate_meas;
      df.pose_change_est = pose_change_est;
      df.planning_init_point_wrt_last_frame =
          planning_init_point_wrt_last_frame.getPose();
      data_frames.push_back(std::move(df));
    }

    // update
    timestamp += TimeStep;
    // update obstacle to next state
    for (Obstacle& ob : obstacles) { // 更新obstacles的状态为下一帧
      Pose ob_pose_next = ob.getPredictPoseAtTimestamp(timestamp);
      ob.setPose(ob_pose_next);
    }

    // update next planning state w.r.t. local frame
    Car next_planning_state_local = planning_init_point_local;
    UpdateNextPlanningStateLocal( // 根据最优轨迹更新下一个时刻ego的状态（在当前时刻ego的local坐标系下）
        best_frenet_path_local, wp_lanes_local[best_frenet_path_local->lane_id],
        &next_planning_state_local);

    // next frame
    Car next_ego_car =
        next_planning_state_local;  // set to next_planning_state_global

    // std::cout << "planning_init_point_local: "
    //           << plannig_init_point_local.getPose()
    //           << planning_init_point_local.getTwist() << std::endl;
    // std::cout << "next_planning_local: " <<
    // next_planning_state_local.getPose()
    //           << next_planning_state_local.getTwist() << std::endl;
    // std::cout << "ego_car: " << ego_car.getPose() << ego_car.getTwist()
    //           << std::endl;

    if (FLAGS_local_planning) { // 将下一个时刻ego的状态转换到global坐标系下，只针对pose，twist和accel还是针对当前车辆坐标系
      ToGlobal(next_planning_state_local, ego_car.getPose(), &next_ego_car);
    }

    // std::cout << "next_ego_car: " << next_ego_car.getPose()
    //           << next_ego_car.getTwist() << std::endl;

    // TODO: remove, update this way not accurate
    // Comment (sz): not accurate due to the way we update sensor measurements
    // UpdateSensorMeasurements(ego_car, next_ego_car, &speed_meas,
    //                          &yaw_rate_meas);
    // // update planning init point
    // pose_change_est = EstimateChangeOfPose(
    //     speed_meas, yaw_rate_meas);  // estimation of state change

    ego_car = next_ego_car;   //更新ego_car的信息为下一个step，pose/twist/accel

    // TEMP FIX: use add-on noise for now. // 得到local坐标系下的姿态偏移量！！！！！
    pose_change_est = EstimateChangeOfPoseV2(next_planning_state_local, // 待理解详细的估计过程 24.11.04
                                             &speed_meas, &yaw_rate_meas);

    planning_init_point_wrt_last_frame = next_planning_state_local;
    // TODO: update Initial Planning Point in local frame
    planning_init_point_local = planning_init_point_wrt_last_frame;
    if (FLAGS_local_planning) {
      planning_init_point_local.setPose( // 下一个周期的规划起点：下一周期车辆在本周期local的位置 - 姿态偏移量
          planning_init_point_wrt_last_frame.getPose() - pose_change_est); // 这个差值就是noise，对应实际的误差
    }

    auto end = std::chrono::high_resolution_clock::now();
    double cycle_duration = get_duration_ms(end, start);
    total_runtime += cycle_duration;

    cout << "\n#" << i << ", simtime: " << timestamp
         << "[s]. Plan: " << plan_duration << "[ms], Cycle: " << cycle_duration
         << "[ms]. x=" << ego_car.getPose().x << ", y=" << ego_car.getPose().y
         << ", yaw=" << utils::rad2deg(ego_car.getPose().yaw)
         << "[deg]. vx=" << ego_car.getTwist().vx
         << ", vy=" << ego_car.getTwist().vy
         << ", w=" << utils::rad2deg(ego_car.getTwist().yaw_rate)
         << "[deg/s]. ax=" << ego_car.getAccel().ax << endl;
    cout << "speed_meas: " << speed_meas
         << ", yaw_rate_meas: " << utils::rad2deg(yaw_rate_meas) << endl;
    cout << "Pose change estimate,  x: " << pose_change_est.x
         << ", y: " << pose_change_est.y
         << ", yaw: " << utils::rad2deg(pose_change_est.yaw) << endl;
  }
  cout << "Total runtime: " << total_runtime << " [ms] for # " << i
       << " iterations." << endl;

  if (FLAGS_store_data) {
    if (save_data(FLAGS_data_path, data_frames)) {
      cout << "Data logging at: " << FLAGS_data_path << endl;
    }
  }
  return 0;
}