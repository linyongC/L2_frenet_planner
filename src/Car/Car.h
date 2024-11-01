#ifndef FRENETOPTIMALTRAJECTORY_CAR_H
#define FRENETOPTIMALTRAJECTORY_CAR_H

#include "utils/geometry.h"
#include "utils/utils.h"

#include <memory>
#include <msgpack.hpp>
#include <tuple>
#include <vector>

using namespace std;
using namespace utils;

// Lincoln MKZ configuration
const double VEHICLE_LENGTH = 4.93;
const double VEHICLE_WIDTH = 1.86;

class Car {
 public:
  Car() : Car({0, 0, 0}) {}
  Car(Pose pose_) : Car(pose_, {0, 0, 0}) {}
  Car(Pose pose_, Twist twist_) : Car(pose_, twist_, {0, 0, 0}) {}
  Car(Pose pose_, Twist twist_, Accel accel_)
      : Car(pose_, twist_, accel_, VEHICLE_LENGTH, VEHICLE_WIDTH) {}
  Car(Pose pose_, Twist twist_, Accel accel_, double length_, double width_)
      : pose(pose_),
        twist(twist_),
        accel(accel_),
        length(length_),
        width(width_) {}

  Box getBox() { return utils::pose_to_box(pose, length, width); }

  void setPose(Pose p) { pose = p; }
  const Pose &getPose() const { return pose; }
  Pose *mutablePose() { return &pose; }

  void setTwist(Twist t) { twist = t; }
  const Twist &getTwist() const { return twist; }
  Twist *mutableTwist() { return &twist; }

  void setAccel(Accel a) { accel = a; }
  const Accel &getAccel() const { return accel; }
  Accel *mutableAccel() { return &accel; }

  void setCurLaneId(int lane_id) { cur_lane_id = lane_id; }
  int getCurLaneId() { return cur_lane_id; }
  void setTargetLaneId(int lane_id) { target_lane_id = lane_id; }
  int getTargetLaneId() { return target_lane_id; }

  //  private:  // sheng: set public to use msgpack
  Pose pose;    // pose w.r.t global frame w.r.t表示with respect to
  Twist twist;  // velocity w.r.t vehicle frame
  Accel accel;  // acceleration w. = .t vehicle frame
  double length;
  double width;

  int cur_lane_id = -1;
  int target_lane_id = -1;

  MSGPACK_DEFINE(pose, twist, accel, length, width);
};

#endif  // FRENETOPTIMALTRAJECTORY_CAR_H
