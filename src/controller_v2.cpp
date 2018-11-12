
// This controller is for ICRA'19

#include <ros/ros.h>
#include <tf/transform_listener.h>
#include <std_srvs/Empty.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/TransformStamped.h>
#include <geometry_msgs/Vector3.h>
#include <kalman/State.h>
#include <cmath>

#include <vector>
#include <fstream>
#include <string>
#include <sstream>

#include "pid.hpp"
#include "MinimumJerkTrajectory.h"
#include "MinimumSnapTrajectory.h"
#include "CircleTrajectory.h"
#include "SineTrajectory.h"
#include "CrazyflieModel.h"

// For GaoModel
#include "Gao.h"
#include "modelLoader.h"
#include "qpUpdate.h"

enum Task_Type
{
  TakeOff = 0,
  GoTo = 1,
  GoHome = 2,
  Land = 3,
  ByeBye = 4,
  Emergency = 5
};

struct task
{
  // type of task
  Task_Type taskType;

  // parameters, different tasks have different parameters, some don't have parameters
  double parameters[5];
};

// 1 :
// 2 : snap
// 3 : anyVelWeight model
#define TRJ 2

double get(
  const ros::NodeHandle& n,
  const std::string& name) {
  double value;
  n.getParam(name, value);
  return value;
}

class Controller
{
public:

  Controller(
    const std::string& worldFrame,
    const std::string& frame,
    const ros::NodeHandle& n)
    : m_worldFrame(worldFrame)
    , m_frame(frame)
    , m_pubNav()
    , m_listener()

    , m_pidX(
        get(n, "PIDs/X/kp"),
        get(n, "PIDs/X/kd"),
        get(n, "PIDs/X/ki"),
        get(n, "PIDs/X/minOutput"),
        get(n, "PIDs/X/maxOutput"),
        get(n, "PIDs/X/integratorMin"),
        get(n, "PIDs/X/integratorMax"),
        "x")
    , m_pidY(
        get(n, "PIDs/Y/kp"),
        get(n, "PIDs/Y/kd"),
        get(n, "PIDs/Y/ki"),
        get(n, "PIDs/Y/minOutput"),
        get(n, "PIDs/Y/maxOutput"),
        get(n, "PIDs/Y/integratorMin"),
        get(n, "PIDs/Y/integratorMax"),
        "y")
    , m_pidZ(
        get(n, "PIDs/Z/kp"),
        get(n, "PIDs/Z/kd"),
        get(n, "PIDs/Z/ki"),
        get(n, "PIDs/Z/minOutput"),
        get(n, "PIDs/Z/maxOutput"),
        get(n, "PIDs/Z/integratorMin"),
        get(n, "PIDs/Z/integratorMax"),
        "z")
    , m_pidYaw(
        get(n, "PIDs/Yaw/kp"),
        get(n, "PIDs/Yaw/kd"),
        get(n, "PIDs/Yaw/ki"),
        get(n, "PIDs/Yaw/minOutput"),
        get(n, "PIDs/Yaw/maxOutput"),
        get(n, "PIDs/Yaw/integratorMin"),
        get(n, "PIDs/Yaw/integratorMax"),
        "yaw")

    , m_state(Idle)
    , m_goal()
    , m_subscribeGoal()

    , m_serviceTakeoff()
    , m_serviceLand()
    , m_serviceGoHome()

    , m_thrust(0)
    , m_startZ(0)
    , X(0.0)
    , Y(0.0)
    , Z(0.0)
    , Vx(0.0)
    , Vy(0.0)
    , Vz(0.0)
    , rawX(0.0)
    , rawY(0.0)
    , rawZ(0.0)
    , rawVx(0.0)
    , rawVy(0.0)
    , rawVz(0.0)
    , startX(0.0)
    , startY(0.0)
    , startZ(0.0)

    , newGoalDeltaX(0.0)
    , newGoalDeltaY(0.0)
    , newGoalDeltaZ(0.0)

    , sineGoalAxis('z')
    , sineGoalOmega(1.0)
    , sineGoalTimes(1)

    , m_subscribeGaosGoal()
    , m_subscribeSineGoal()
    , m_subscribeCircleGoal()

    , gravityCompensation(42500)
    , trajectoryStartTime(ros::Time::now().toSec())
    , trajectoryTrackTime(0.0)
    , hoverStateCounter(0)

    , m_pidXTrack(
        get(n, "PIDs/XTrack/kp"),
        get(n, "PIDs/XTrack/kd"),
        get(n, "PIDs/XTrack/ki"),
        get(n, "PIDs/XTrack/minOutput"),
        get(n, "PIDs/XTrack/maxOutput"),
        get(n, "PIDs/XTrack/integratorMin"),
        get(n, "PIDs/XTrack/integratorMax"),
        "z_track")
    , m_pidYTrack(
        get(n, "PIDs/YTrack/kp"),
        get(n, "PIDs/YTrack/kd"),
        get(n, "PIDs/YTrack/ki"),
        get(n, "PIDs/YTrack/minOutput"),
        get(n, "PIDs/YTrack/maxOutput"),
        get(n, "PIDs/YTrack/integratorMin"),
        get(n, "PIDs/YTrack/integratorMax"),
        "z_track")
    , m_pidZTrack(
        get(n, "PIDs/ZTrack/kp"),
        get(n, "PIDs/ZTrack/kd"),
        get(n, "PIDs/ZTrack/ki"),
        get(n, "PIDs/ZTrack/minOutput"),
        get(n, "PIDs/ZTrack/maxOutput"),
        get(n, "PIDs/ZTrack/integratorMin"),
        get(n, "PIDs/ZTrack/integratorMax"),
        "z_track")

    , m_autoState(Not_Ready_To_Track)
    , goalShallNotChange(false)

    , newTrajectoryAwaits(false)
    , newSineTrajectoryAwaits(false)
    , newCircleTrajectoryAwaits(false)

    , goHomeState(noWannaGoHome)
    , trj()
#if TRJ == 2
    , GaoModel("/home/sun/Downloads/rmuAnyVelWeightRDUCorrect_7_500_317.json")
    , trj_time()
    , weight(get(n, "Aggs"))
#elif TRJ == 3
    , GaoModel("/home/sun/Downloads/rmuAnyVelWeightRDUCorrect_7_500_317.json")
    , weight(get(n, "Aggs"))
#endif
    , trj_sine()
    , tasks()
    , trajectoryType(Linear)
    , sineLog()
    , logSequence(get(n, "logBase"))
    , sineLogSequence(get(n, "sineLogBase"))
    , sineXLogSequence(get(n, "sineXLogBase"))
    , sineYLogSequence(get(n, "sineYLogBase"))
    , sineZLogSequence(get(n, "sineZLogBase"))
    , averageVelocity(get(n, "SnapVel"))
  {
    ros::NodeHandle nh;
    m_listener.waitForTransform(m_worldFrame, m_frame, ros::Time(0), ros::Duration(10.0));
    m_pubNav = nh.advertise<geometry_msgs::Twist>("cmd_vel", 1);
    m_subscribeGoal = nh.subscribe("goal", 1, &Controller::goalChanged, this);

    m_serviceTakeoff = nh.advertiseService("takeoff", &Controller::takeoff, this);
    m_serviceLand = nh.advertiseService("land", &Controller::land, this);
    m_serviceGoHome = nh.advertiseService("goHome", &Controller::goHome, this);

    m_filteredState = nh.subscribe("/filteredState", 1 , &Controller::filteredStateReceived, this);
    m_rawState = nh.subscribe("/rawState", 1 , &Controller::rawStateReceived, this);

    m_subscribeGaosGoal = nh.subscribe("gaosGoal", 1, &Controller::gaosNewGoalReceived, this);
    m_subscribeSineGoal = nh.subscribe("sineGoal", 1, &Controller::newSineGoalReceived, this);
    m_subscribeCircleGoal = nh.subscribe("circleGoal", 1, &Controller::newCircleGoalReceived, this);

    // m_subscribeLogSequnceBase = nh.subscribe("logSequenceBase", 1, &Controller::updateLogSequenceBase, this);
  }

  void run(double frequency)
  {
    ros::NodeHandle node;
    ros::Timer timer = node.createTimer(ros::Duration(1.0 / frequency), &Controller::iteration, this);
    ros::spin();
  }

private:

  void filteredStateReceived(const kalman::State& msg)
  {
    this -> X = msg.position.x;
    this -> Y = msg.position.y;
    this -> Z = msg.position.z;

    this -> Vx = msg.velocity.x;
    this -> Vy = msg.velocity.y;
    this -> Vz = msg.velocity.z;
  }

  void rawStateReceived(const kalman::State& msg)
  {
    this -> rawX = msg.position.x;
    this -> rawY = msg.position.y;
    this -> rawZ = msg.position.z;

    this -> rawVx = msg.velocity.x;
    this -> rawVy = msg.velocity.y;
    this -> rawVz = msg.velocity.z;
  }

  void goalChanged(
    const geometry_msgs::PoseStamped::ConstPtr& msg)
  {
    if (goalShallNotChange == false)
    {
      m_goal = *msg;
      goalShallNotChange = true;

      startX = m_goal.pose.position.x;
      startY = m_goal.pose.position.y;
      startZ = m_goal.pose.position.z;

      ROS_ERROR("startXYZ set: (%.2f, %.2f, %.2f)", startX, startY, startZ);
    }
  }

  void gaosNewGoalReceived(const geometry_msgs::Vector3& msg)
  {
    newGoalDeltaX = msg.x;
    newGoalDeltaY = msg.y;
    newGoalDeltaZ = msg.z;

    newTrajectoryAwaits = true;

    ROS_ERROR("New goal(delta) received: (%.2f, %.2f, %.2f)", newGoalDeltaX, newGoalDeltaY, newGoalDeltaZ);
  }

  void newSineGoalReceived(const geometry_msgs::Vector3& msg)
  {
    if (floor(msg.x) == 1)
      sineGoalAxis = 'x';
    else if (floor(msg.x) == 2)
      sineGoalAxis = 'y';
    else if (floor(msg.x) == 3)
      sineGoalAxis = 'z';
    else
      sineGoalAxis = 'z';

    sineGoalOmega = msg.y;
    sineGoalTimes = msg.z;
    newSineTrajectoryAwaits = true;

    ROS_ERROR("New sine trajectory: Axis: %c, omega: %f, time: %f", sineGoalAxis, sineGoalOmega, sineGoalTimes);

  }

  void newCircleGoalReceived(const geometry_msgs::Vector3& msg)
  {

  }

  bool takeoff(
    std_srvs::Empty::Request& req,
    std_srvs::Empty::Response& res)
  {
    ROS_INFO("Takeoff requested!");
    m_state = TakingOff;

    tf::StampedTransform transform;
    m_listener.lookupTransform(m_worldFrame, m_frame, ros::Time(0), transform);
    m_startZ = transform.getOrigin().z();
    m_startZ = Z;

    return true;
  }

  bool land(
    std_srvs::Empty::Request& req,
    std_srvs::Empty::Response& res)
  {
    ROS_INFO("Landing requested!");
    m_state = Landing;
    m_pidZ.reset(); //changed by GAO
    return true;
  }

  bool goHome(
    std_srvs::Empty::Request& req,
    std_srvs::Empty::Response& res)
  {
    ROS_INFO("Go home requested!");

    goHomeState = wannaGoHome;
    newTrajectoryAwaits = true;

    return true;
  }

  void getTransform(
    const std::string& sourceFrame,
    const std::string& targetFrame,
    tf::StampedTransform& result)
  {
    m_listener.lookupTransform(sourceFrame, targetFrame, ros::Time(0), result);
  }

  void pidReset()
  {
    m_pidX.reset();
    m_pidY.reset();
    m_pidZ.reset();
    m_pidYaw.reset();

    m_pidXTrack.reset();
    m_pidYTrack.reset();
    m_pidZTrack.reset();
  }

  void iteration(const ros::TimerEvent& e)
  {
    float dt = e.current_real.toSec() - e.last_real.toSec();

    switch (m_state)
    {
    case TakingOff:
    {
      tf::StampedTransform transform;
      m_listener.lookupTransform(m_worldFrame, m_frame, ros::Time(0), transform);
      if (this->Z > m_startZ + 0.05 || m_thrust > 50000)
      {
        gravityCompensation = m_thrust - 1500;
        ROS_ERROR("gravity compensation: %f", gravityCompensation);
        pidReset();
        m_state = Automatic;
        m_thrust = 0;
      }
      else
      {
        if (m_thrust < gravityCompensation)
          m_thrust = gravityCompensation;

        m_thrust += 1000 * dt;
        geometry_msgs::Twist msg;
        msg.linear.z = m_thrust;
        m_pubNav.publish(msg);
      }

    }
    break;
    case Landing:
    {
      m_goal.pose.position.z = m_startZ + 0.05;
      tf::StampedTransform transform;
      m_listener.lookupTransform(m_worldFrame, m_frame, ros::Time(0), transform);
      if (Z <= m_startZ + 0.05) {
        m_state = Idle;
        geometry_msgs::Twist msg;
        m_pubNav.publish(msg);

        m_goal.pose.position.z = startZ;
      }
    }
    // intentional fall-thru
    case Automatic:
    {
      tf::StampedTransform transform;
      m_listener.lookupTransform(m_worldFrame, m_frame, ros::Time(0), transform);

      geometry_msgs::PoseStamped targetWorld;
      targetWorld.header.stamp = transform.stamp_;
      targetWorld.header.frame_id = m_worldFrame;
      targetWorld.pose = m_goal.pose;

      geometry_msgs::PoseStamped targetDrone;
      m_listener.transformPose(m_frame, targetWorld, targetDrone);

      tfScalar roll, pitch, yaw;
      tf::Matrix3x3(
        tf::Quaternion(
          targetDrone.pose.orientation.x,
          targetDrone.pose.orientation.y,
          targetDrone.pose.orientation.z,
          targetDrone.pose.orientation.w
        )).getRPY(roll, pitch, yaw);

      geometry_msgs::Twist msg;
      targetDrone.pose.position.x = m_goal.pose.position.x - X;
      targetDrone.pose.position.y = m_goal.pose.position.y - Y;
      targetDrone.pose.position.z = m_goal.pose.position.z - Z;

      // ROS_WARN("goal x: %f, y: %f, z: %f", m_goal.pose.position.x, m_goal.pose.position.y, m_goal.pose.position.z);

      // add hover mode
      double positionError = sqrt(targetDrone.pose.position.x * targetDrone.pose.position.x
                                  + targetDrone.pose.position.y * targetDrone.pose.position.y
                                  + targetDrone.pose.position.z * targetDrone.pose.position.z);

      double velocityError = sqrt(Vx * Vx + Vy * Vy + Vz * Vz);

      if (m_autoState != Tracking && positionError < 0.05 && velocityError < 0.6)
      {
        hoverStateCounter++;
      }
      else
      {
        hoverStateCounter = 0;

        if (m_autoState == Ready_To_Track) // not ready to track anymore
          m_autoState = Not_Ready_To_Track;
      }

      if (hoverStateCounter >= 20)
      {
        // gravity compenstation reset
        gravityCompensation += m_pidZ.getIntegral() * m_pidZ.ki();
        m_pidZ.setIntegral(0.0);

        // jump into Ready_To_Track state
        if (m_autoState == Not_Ready_To_Track)
          m_autoState = Ready_To_Track;
        //else if (m_autoState == Recovering_From_Tracking)
        //{
        //m_autoState = Ready_To_Track;
        //ROS_ERROR("Recovered from tracking!");
        //}
      }

      if ( m_autoState == Recovering_From_Tracking)
      {
        sineLog << ros::Time::now().toSec() - trajectoryStartTime << ", " << m_goal.pose.position.x << ", " << m_goal.pose.position.y << ", " << m_goal.pose.position.z << ", " << X << ", " << Y << ", " << Z  << ", " << "0.0" << ", " << "0.0" << ", " << "0.0" << ", " << rawVx << ", " << rawVy << ", " << rawVz << ", " << Vx << ", " << Vy << ", " << Vz << ", " << msg.linear.z << ", " << msg.linear.x << ", " << msg.linear.y << ", " << "0.0" << ", " << "0.0" << ", " << "0.0" << std::endl;
        // ROS_WARN("velE: %f, xe: %f, ye: %f, ze: %f", velocityError, abs(targetDrone.pose.position.x), abs(targetDrone.pose.position.y), abs(targetDrone.pose.position.z));
      }

      if ( m_autoState == Recovering_From_Tracking && velocityError < 0.04 && abs(targetDrone.pose.position.x) < 0.02 && abs(targetDrone.pose.position.y) < 0.02 && abs(targetDrone.pose.position.z) < 0.02 )
      {
        m_autoState = Ready_To_Track;
        ROS_ERROR("Recovered from tracking!");
        sineLog.close();
      }

      //Finite State machines
      if (m_autoState == Not_Ready_To_Track)
      {
        msg.linear.x = m_pidX.update(targetDrone.pose.position.x);
        msg.linear.y = m_pidY.update(targetDrone.pose.position.y);
        msg.linear.z = gravityCompensation + m_pidZ.update(targetDrone.pose.position.z, 0 - Vz);
      }
      else if (m_autoState == Recovering_From_Tracking)
      {
        double targetPx = m_goal.pose.position.x, targetPy = m_goal.pose.position.y, targetPz = m_goal.pose.position.z;

        double targetVx = 0.0, targetVy = 0.0, targetVz = 0.0, targetAx = 0.0, targetAy = 0.0, targetAz = 0.0;

        //double ddrx = m_pidXTrack.update((targetPx - X), (targetVx - Vx));
        //double ddry = m_pidYTrack.update((targetPy - Y), (targetVy - Vy));

        //msg.linear.z =  Thrust2PWM(PWM2Thrust(gravityCompensation) +
        //                           PWM2Thrust(gravityCompensation) / 9.81 * (targetAz + m_pidZTrack.update(targetPz - Z, targetVz - Vz)));

        //msg.linear.x = (ddrx * cos(yaw) + ddry * sin(yaw)) / 9.81;
        //msg.linear.y = (ddrx * sin(yaw) - ddry * cos(yaw)) / 9.81;

        msg.linear.x = m_pidX.update(targetDrone.pose.position.x);
        msg.linear.y = m_pidY.update(targetDrone.pose.position.y);
        msg.linear.z = gravityCompensation + m_pidZ.update(targetDrone.pose.position.z, 0 - Vz);


        //ROS_WARN("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f",
        //         ros::Time::now().toSec() - trajectoryStartTime,
        //         targetPx, targetPy, targetPz,
        //         X, Y, Z,
        //         targetVx, targetVy, targetVz,
        //         rawVx, rawVy, rawVz,
        //         Vx, Vy, Vz,
        //         msg.linear.z, msg.linear.x, msg.linear.y,
        //         targetAx, targetAy, targetAz);

      }
      else if (m_autoState == Ready_To_Track)
      {
        msg.linear.x = m_pidX.update(targetDrone.pose.position.x);
        msg.linear.y = m_pidY.update(targetDrone.pose.position.y);
        msg.linear.z = gravityCompensation + m_pidZ.update(targetDrone.pose.position.z, 0 - Vz);

        if (goHomeState == haveArrivedAtHome)
        {
          // land
          ROS_INFO("Start landing!");
          m_state = Landing;
          m_pidZ.reset();
          goHomeState = noWannaGoHome;
        }

        if (newTrajectoryAwaits == true)
        {
          std::ostringstream fileName;
          fileName << "/home/sun/Desktop/icra19/log/test/" << "test_" << logSequence << ".csv";
          sineLog.open(fileName.str());

#if TRJ == 2
          double wayPointsX[2] = {X, m_goal.pose.position.x + newGoalDeltaX};
          double wayPointsY[2] = {Y, m_goal.pose.position.y + newGoalDeltaY};
          double wayPointsZ[2] = {Z, m_goal.pose.position.z + newGoalDeltaZ};

          if (goHomeState == wannaGoHome)
          {
            wayPointsX[1] = startX;
            wayPointsY[1] = startY;
            wayPointsZ[1] = startZ;
            goHomeState = onTheWayHome;
          }

          //ROS_WARN("X, goal, delta: %f, %f, %f", X, m_goal.pose.position.x, newGoalDeltaX);
          //ROS_WARN("Y, goal, delta: %f, %f, %f", Y, m_goal.pose.position.y, newGoalDeltaY);
          //ROS_WARN("Z, goal, delta: %f, %f, %f", Z, m_goal.pose.position.z, newGoalDeltaZ);
          

          // Code below are usd to generate a suitable time
          double xin[7] = {-(wayPointsX[1] - wayPointsX[0]), -(wayPointsY[1] - wayPointsY[0]), -(wayPointsZ[1] - wayPointsZ[0]), 0.0, 0.0, 0.0, log10(weight)};
          VX vec = GaoModel.evalQP(xin);
          MX acce = acceGen(vec);
          trj_time.generateTrajectory(vec, acce);
          double distanceToNextGoal = sqrt(pow(wayPointsX[1] - wayPointsX[0],2) + pow(wayPointsY[1]-wayPointsY[0],2) + pow(wayPointsZ[1]-wayPointsZ[0],2));
          double snapVelocity = distanceToNextGoal / trj_time.getDuration();
          ROS_ERROR("snap vel: %f", snapVelocity);
          sineLog << "# snap vel: " << snapVelocity << endl;
          sineLog << "# gt " << wayPointsX[1] - wayPointsX[0] << " " << wayPointsY[1] - wayPointsY[0] << " " << wayPointsZ[1] - wayPointsZ[0] << endl;

          trj.generateTrajectory(wayPointsX, wayPointsY, wayPointsZ, 2, snapVelocity);
#elif TRJ == 3
          double xin[7] = { -(m_goal.pose.position.x + newGoalDeltaX - X),
                            -(m_goal.pose.position.y + newGoalDeltaY - Y),
                            -(m_goal.pose.position.z + newGoalDeltaZ - Z),
                            Vx, Vy, Vz, log10(weight)
                          }; // initial point

          if (goHomeState == wannaGoHome)
          {
            xin[0] = -(startX - X);
            xin[1] = -(startY - Y);
            xin[2] = -(startZ - Z);
            goHomeState = onTheWayHome;
          }

          VX vec = GaoModel.evalQP(xin);
          MX acce = acceGen(vec);

          trj.generateTrajectory(vec, acce);

          double refY = Y;
          double refX = X;
          double refZ = Z;

          trj.setPositionOffset(refX - xin[0], refY - xin[1], refZ - xin[2]);
          sineLog << "# Aggr: " << weight << endl;
          sineLog << "# gt " << -xin[0] << " " << -xin[1] << " " << -xin[2] << endl;

#endif
          trajectoryType = Linear;

          hoverStateCounter = 0;
          m_autoState = Tracking;
          newTrajectoryAwaits = false;
          ROS_ERROR("new linear trajectory generated!");
          trajectoryStartTime = ros::Time::now().toSec();
        }
        else if (newSineTrajectoryAwaits == true)
        {
          std::ostringstream fileName;
          if (sineGoalAxis == 'x')
          {
            fileName << "/home/sun/Desktop/icra19/log/sinex/" << "sinex_" << sineXLogSequence << ".csv";
            sineXLogSequence++;
          }
          else if (sineGoalAxis == 'y')
          {
            fileName << "/home/sun/Desktop/icra19/log/siney/" << "siney_" << sineYLogSequence << ".csv";
            sineYLogSequence++;
          }
          else 
          {
            fileName << "/home/sun/Desktop/icra19/log/sinez/" << "sinez_" << sineZLogSequence << ".csv";
            sineZLogSequence++;
          }

          // fileName << "/home/sun/Desktop/icra19/log/" << "sine_" << sineLogSequence << ".csv";
          
          sineLog.open(fileName.str());

          // std::cout << "Start logging in: " << fileName.str() << std::endl;
          // sineLog << "New Sine Trj, omega: " << sineGoalOmega << std::endl;

          trj_sine.generateTrajectory(m_goal.pose.position.x, m_goal.pose.position.y, m_goal.pose.position.z,
                                      sineGoalAxis, sineGoalOmega, sineGoalTimes);

          sineLog << "# Omega: " << trj_sine.omega() << " Axis: " << trj_sine.axis() << endl;

          trajectoryType = Sine;

          hoverStateCounter = 0;
          m_autoState = Tracking;
          newSineTrajectoryAwaits = false;
          ROS_ERROR("new sine trajectory generated!");
          trajectoryStartTime = ros::Time::now().toSec();
        }
        else if (newCircleTrajectoryAwaits == true)
        {
          ROS_ERROR("new circular trajectory generated!");

        }
      }
      else if (m_autoState == Tracking)
      {
        if (0 && (abs(targetDrone.pose.position.x) >= 0.06 || abs(targetDrone.pose.position.y) >= 0.06 || abs(targetDrone.pose.position.z) >= 0.06))
        {
          ROS_ERROR("position way off!, planning new one");
        }

        trajectoryTrackTime = ros::Time::now().toSec() - trajectoryStartTime;
        double targetPx, targetPy, targetPz, targetVx, targetVy, targetVz, targetAx, targetAy, targetAz, targetPhi, targetTheta, targetPsi;

        if (trajectoryType == Sine)
        {
          trj_sine.sampleTrajectory(trajectoryTrackTime, targetPx, targetPy, targetPz, targetVx, targetVy, targetVz, targetAx, targetAy, targetAz, targetPhi, targetTheta, targetPsi);
          ROS_WARN("%f, %f, %f, %f, %f, %f, %f", trajectoryTrackTime, targetPx, targetPy, targetPz, X, Y, Z);

        }
        else if (trajectoryType == Circular)
          int placeHolder = 1;
        else
          trj.sampleTrajectory(trajectoryTrackTime, targetPx, targetPy, targetPz, targetVx, targetVy, targetVz, targetAx, targetAy, targetAz, targetPhi, targetTheta, targetPsi);

        m_goal.pose.position.x = targetPx;
        m_goal.pose.position.y = targetPy;
        m_goal.pose.position.z = targetPz;

        double ddrx = 57.3 * targetAx + m_pidXTrack.update((targetPx - X), (targetVx - Vx));
        double ddry = 57.3 * targetAy + m_pidYTrack.update((targetPy - Y), (targetVy - Vy));


        msg.linear.z =  Thrust2PWM(PWM2Thrust(gravityCompensation) +
                                   PWM2Thrust(gravityCompensation) / 9.81 * (targetAz + m_pidZTrack.update(targetPz - Z, targetVz - Vz)));

        msg.linear.x = (ddrx * cos(yaw) + ddry * sin(yaw)) / 9.81 * cos(targetPhi) * cos(targetTheta);
        msg.linear.y = (ddrx * sin(yaw) - ddry * cos(yaw)) / 9.81 * cos(targetPhi) * cos(targetTheta);

        msg.linear.z = msg.linear.z / cos(targetPhi) / cos(targetTheta);

        if (msg.linear.x < -30)
          msg.linear.x = -30;
        else if (msg.linear.x > 30)
          msg.linear.x = 30;

        if (msg.linear.y < -30)
          msg.linear.y = -30;
        else if (msg.linear.y > 30)
          msg.linear.y = 30;

        // if (trajectoryType == Sine)
        if (1)
          sineLog << trajectoryTrackTime << ", " << targetPx << ", " << targetPy << ", " << targetPz << ", " << X << ", " << Y << ", " << Z  << ", " << targetVx << ", " << targetVy << ", " << targetVz << ", " << rawVx << ", " << rawVy << ", " << rawVz << ", " << Vx << ", " << Vy << ", " << Vz << ", " << msg.linear.z << ", " << msg.linear.x << ", " << msg.linear.y << ", " << targetAx << ", " << targetAy << ", " << targetAz << std::endl;

        // ROS_WARN("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f", trajectoryTrackTime, targetPx, targetPy, targetPz, X, Y, Z, targetVx, targetVy, targetVz, rawVx, rawVy, rawVz, Vx, Vy, Vz, msg.linear.z, msg.linear.x, msg.linear.y, targetAx, targetAy, targetAz);

        hoverStateCounter = 0;

        double endPx, endPy, endPz;
        if (trajectoryType == Sine)
          trj_sine.getEndPosition(endPx, endPy, endPz);
        else if (trajectoryType == Circular)
          int placeHolder = 0;
        else
          trj.getEndPosition(endPx, endPy, endPz);

        double finalPositionError = sqrt((endPx - X) * (endPx - X) + (endPy - Y) * (endPy - Y) + (endPz - Z) * (endPz - Z));

        if (trajectoryType == Linear && (finalPositionError < 0.05 || trajectoryTrackTime >= trj.getEndTime()))
        {
          // sineLog.close();
          logSequence++;

          ROS_ERROR("Trajectory tracking complete!");

          m_goal.pose.position.x = endPx;
          m_goal.pose.position.y = endPy;
          m_goal.pose.position.z = endPz;

          m_autoState = Not_Ready_To_Track;
          m_autoState = Recovering_From_Tracking;

          m_pidXTrack.reset();
          m_pidYTrack.reset();
          m_pidZTrack.reset();

          if (goHomeState == onTheWayHome)
            goHomeState = haveArrivedAtHome;
        }

        if (trajectoryType == Sine && finalPositionError < 0.15 && trajectoryTrackTime >= trj_sine.getEndTime())
        {
          // sineLog.close();

          ROS_ERROR("Sine trajectory tracking complete!");

          m_goal.pose.position.x = endPx;
          m_goal.pose.position.y = endPy;
          m_goal.pose.position.z = endPz;

          m_autoState = Not_Ready_To_Track;
          m_autoState = Recovering_From_Tracking;

          m_pidXTrack.reset();
          m_pidYTrack.reset();
          m_pidZTrack.reset();

        }

      }
      else
      {
        ROS_ERROR("Automatic sub state error!");
        exit(0);
      }

      if (msg.linear.z > 65000)
      {
        msg.linear.z = 65000;
      }

      if (msg.linear.z < 10000)
      {
        msg.linear.z = 10000;
      }

      msg.angular.z = m_pidYaw.update(yaw);
      m_pubNav.publish(msg);
    }
    break;

    case Idle:
    {
      geometry_msgs::Twist msg;
      m_pubNav.publish(msg);
      m_autoState = Not_Ready_To_Track;

      m_pidX.reset();
      m_pidY.reset();
      m_pidZ.reset();
      m_pidYaw.reset();
    }
    break;
    }
  }

private:

  enum State
  {
    Idle = 0,
    Automatic = 1,
    TakingOff = 2,
    Landing = 3,
  };

  enum AutoState
  {
    Tracking = 0,
    Recovering_From_Tracking = 1,
    Ready_To_Track = 2,
    Not_Ready_To_Track = 3,
  };

  enum GoHomeState
  {
    noWannaGoHome = 0, // default value
    wannaGoHome = 1,
    onTheWayHome = 2,
    haveArrivedAtHome = 3,
  };

  enum TrajectoryType
  {
    Linear = 0,
    Sine = 1,
    Circular = 2,
  };

private:
  std::string m_worldFrame;
  std::string m_frame;
  ros::Publisher m_pubNav;
  tf::TransformListener m_listener;
  PID m_pidX;
  PID m_pidY;
  PID m_pidZ;
  PID m_pidYaw;
  State m_state;
  geometry_msgs::PoseStamped m_goal;
  ros::Subscriber m_subscribeGoal;
  ros::ServiceServer m_serviceTakeoff;
  ros::ServiceServer m_serviceLand;
  ros::ServiceServer m_serviceGoHome;
  float m_thrust;
  float m_startZ;

  ros::Subscriber m_filteredState;
  ros::Subscriber m_rawState;

  ros::Subscriber m_subscribeGaosGoal;
  ros::Subscriber m_subscribeSineGoal;
  ros::Subscriber m_subscribeCircleGoal;

  double X;
  double Y;
  double Z;

  double Vx;
  double Vy;
  double Vz;

  double rawX;
  double rawY;
  double rawZ;

  double rawVx;
  double rawVy;
  double rawVz;

  // start point's coordinates in XYZ
  double startX;
  double startY;
  double startZ;

  // for linear goal
  double newGoalDeltaX;
  double newGoalDeltaY;
  double newGoalDeltaZ;

  // for sine trajectory goal
  char sineGoalAxis;
  double sineGoalOmega;
  double sineGoalTimes;

  float gravityCompensation;

  double trajectoryStartTime;
  double trajectoryTrackTime;

  int hoverStateCounter;

  PID m_pidXTrack;
  PID m_pidYTrack;
  PID m_pidZTrack;

  //CircleTrajectory trj3;
#if TRJ == 2
  TDMST trj;
  mlpModelQP GaoModel;
  Gao trj_time;
  double weight;
#elif TRJ == 3
  mlpModelQP GaoModel;
  Gao trj;
  double weight;
#endif

  SineTrajectory trj_sine;

  AutoState m_autoState;
  bool goalShallNotChange;

  bool newTrajectoryAwaits;
  bool newSineTrajectoryAwaits;
  bool newCircleTrajectoryAwaits;

  GoHomeState goHomeState;

  std::vector<task> tasks;
  TrajectoryType trajectoryType;

  std::ofstream sineLog;

  unsigned int logSequence;
  unsigned int sineLogSequence;
  unsigned int sineXLogSequence;
  unsigned int sineYLogSequence;
  unsigned int sineZLogSequence;

  double averageVelocity;

};

int main(int argc, char **argv)
{
  ros::init(argc, argv, "controller");

  // Read parameters
  ros::NodeHandle n("~");
  std::string worldFrame;
  n.param<std::string>("worldFrame", worldFrame, "/world");
  std::string frame;
  n.getParam("frame", frame);
  double frequency;
  n.param("frequency", frequency, 100.0);

  Controller controller(worldFrame, frame, n);
  controller.run(frequency);

  return 0;
}
