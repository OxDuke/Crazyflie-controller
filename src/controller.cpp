#include <ros/ros.h>
#include <tf/transform_listener.h>
#include <std_srvs/Empty.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/TransformStamped.h>
#include <geometry_msgs/Vector3.h>
#include <cmath>

#include <fstream>
#include <string>

#include "pid.hpp"
#include "MinimumJerkTrajectory.h"
#include "MinimumSnapTrajectory.h"
#include "CircleTrajectory.h"
#include "Gao.h"
#include "CrazyflieModel.h"
#include "modelLoader.h"
#include "qpUpdate.h"

// 1 : anyVel model
// 2 : snap
// 3 : anyVelWeight model
#define TRJ 3

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
    , m_serivceGoHome()
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
    , m_subscribeGaosGoal()
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
    , trj()
#if TRJ == 1
    , GaoModel("/home/sun/Downloads/anyVel_6_200_317.json")
#elif TRJ == 2
    , snapVelocity(get(n, "SnapVel"))
#elif TRJ == 3
    , GaoModel("/home/sun/Downloads/rmuAnyVelWeightRDUCorrect_7_500_317.json")
    , weight(get(n, "Aggs"))
#endif
    , logFile()
    , logSequence(get(n, "logBase"))
  {
    ros::NodeHandle nh;
    m_listener.waitForTransform(m_worldFrame, m_frame, ros::Time(0), ros::Duration(10.0));
    m_pubNav = nh.advertise<geometry_msgs::Twist>("cmd_vel", 1);
    m_subscribeGoal = nh.subscribe("goal", 1, &Controller::goalChanged, this);
    m_serviceTakeoff = nh.advertiseService("takeoff", &Controller::takeoff, this);
    m_serviceLand = nh.advertiseService("land", &Controller::land, this);
    m_serivceGoHome = nh.advertiseService("goHome", &Controller::goHome, this);

    m_Pos = nh.subscribe("/Pos", 1, &Controller::positionReceived, this);
    m_Vel = nh.subscribe("/Vel", 1, &Controller::velocityReceived, this);
    m_rawPos = nh.subscribe("/RawPos", 1, &Controller::rawPositionReceived, this);
    m_rawVel = nh.subscribe("/RawVel", 1, &Controller::rawVelocityReceived, this);

    m_subscribeGaosGoal = nh.subscribe("gaosGoal", 1, & Controller::gaosNewGoalReceived, this);
  }

  void run(double frequency)
  {
    ros::NodeHandle node;
    ros::Timer timer = node.createTimer(ros::Duration(1.0 / frequency), &Controller::iteration, this);
    ros::spin();
  }

private:
  void positionReceived(const geometry_msgs::Vector3& msg)
  {
    this->X = msg.x;
    this->Y = msg.y;
    this->Z = msg.z;
  }

  void velocityReceived(const geometry_msgs::Vector3& msg)
  {
    this->Vx = msg.x;
    this->Vy = msg.y;
    this->Vz = msg.z;
  }

  void rawPositionReceived(const geometry_msgs::Vector3& msg)
  {
    this->rawX = msg.x;
    this->rawY = msg.y;
    this->rawZ = msg.z;
  }

  void rawVelocityReceived(const geometry_msgs::Vector3& msg)
  {
    this->rawVx = msg.x;
    this->rawVy = msg.y;
    this->rawVz = msg.z;
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

      //ROS_ERROR("startXYZ set: (%.2f, %.2f, %.2f)", startX, startY, startZ);
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
    newGoalDeltaX = startX - m_goal.pose.position.x;
    newGoalDeltaY = startY - m_goal.pose.position.y;
    newGoalDeltaZ = startZ - m_goal.pose.position.z;

    newTrajectoryAwaits = true;

    ROS_INFO("Go home requested!");
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

      // add hover mode
      double positionError = sqrt(targetDrone.pose.position.x * targetDrone.pose.position.x
                                  + targetDrone.pose.position.y * targetDrone.pose.position.y
                                  + targetDrone.pose.position.z * targetDrone.pose.position.z);

      double velocityError = sqrt(Vx * Vx + Vy * Vy + Vz * Vz);

      if (positionError < 0.08 && velocityError < 0.10 && m_autoState != Tracking)
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

      //if ( m_autoState == Recovering_From_Tracking)
      //{
      //  ROS_ERROR("Rec: p_e: %f, %f, %f, v_e: %f", targetDrone.pose.position.x, targetDrone.pose.position.y, targetDrone.pose.position.z, velocityError);
      //}

      if (m_autoState == Recovering_From_Tracking && velocityError < 0.10 && abs(targetDrone.pose.position.x) < 0.015 && abs(targetDrone.pose.position.y) < 0.015 && abs(targetDrone.pose.position.z) < 0.015 )
      {
        m_autoState = Ready_To_Track;
        ROS_ERROR("Recovered from tracking!");

        logFile.close();
        logSequence++;
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
        
        logFile << ros::Time::now().toSec() - trajectoryStartTime << ", " << targetPx << ", " << targetPy << ", " << targetPz << ", " << X << ", " << Y << ", " << Z  << ", " << targetVx << ", " << targetVy << ", " << targetVz << ", " << rawVx << ", " << rawVy << ", " << rawVz << ", " << Vx << ", " << Vy << ", " << Vz << ", " << msg.linear.z << ", " << msg.linear.x << ", " << msg.linear.y << ", " << targetAx << ", " << targetAy << ", " << targetAz << std::endl;


      }
      else if (m_autoState == Ready_To_Track)
      {
        msg.linear.x = m_pidX.update(targetDrone.pose.position.x);
        msg.linear.y = m_pidY.update(targetDrone.pose.position.y);
        msg.linear.z = gravityCompensation + m_pidZ.update(targetDrone.pose.position.z, 0 - Vz);
        if (newTrajectoryAwaits == true)
        {
          std::ostringstream fileName;
          fileName << "/home/sun/Desktop/iros18/" << "log_" << logSequence << ".csv";
          logFile.open(fileName.str());

#if TRJ == 2
          double wayPointsX[2] = {X, m_goal.pose.position.x + newGoalDeltaX};
          double wayPointsY[2] = {Y, m_goal.pose.position.y + newGoalDeltaY};
          double wayPointsZ[2] = {Z, m_goal.pose.position.z + newGoalDeltaZ};
          trj.generateTrajectory(wayPointsX, wayPointsY, wayPointsZ, 2, snapVelocity);

#elif TRJ == 1
          if (newGoalDeltaX > 300 && newGoalDeltaY > 300 && newGoalDeltaZ > 300)
          {
            trj.generateTrajectory("/home/sun/Downloads/threeTraj/goto_2_1000.txt");
            trj.setPositionOffset(3.0 + X, 3.0 + Y, 1.0 + Z);
          }
          else if (newGoalDeltaX > 255 && newGoalDeltaY > 255 && newGoalDeltaZ > 255)
          {
            trj.generateTrajectory("/home/sun/Downloads/threeTraj/goback_2_1000.txt");
            trj.setPositionOffset(-3.0 + X, -3.0 + Y, -1.0 + Z);
          }
          else
          {
            double xin[6] = { -newGoalDeltaX, -newGoalDeltaY, -newGoalDeltaZ, Vx, Vy, Vz}; // initial point
            VX vec = GaoModel.evalQP(xin);

            trj.generateTrajectory(vec);

            double refY = Y;
            double refX = X;
            double refZ = Z;

            trj.setPositionOffset(refX - xin[0], refY - xin[1], refZ - xin[2]);
          }
#elif TRJ == 3
          double xin[7] = { -(m_goal.pose.position.x + newGoalDeltaX - X),
                            -(m_goal.pose.position.y + newGoalDeltaY - Y),
                            -(m_goal.pose.position.z + newGoalDeltaZ - Z),
                            Vx, Vy, Vz, log10(weight)
                          }; // initial point
          VX vec = GaoModel.evalQP(xin);
          MX acce = acceGen(vec);

          trj.generateTrajectory(vec, acce);

          double refY = Y;
          double refX = X;
          double refZ = Z;

          trj.setPositionOffset(refX - xin[0], refY - xin[1], refZ - xin[2]);
#endif
          hoverStateCounter = 0;
          m_autoState = Tracking;
          newTrajectoryAwaits = false;
          ROS_ERROR("new trajectory generated!");
          trajectoryStartTime = ros::Time::now().toSec();
        }
      }
      else if (m_autoState == Tracking)
      {
#if TRJ == 1
        //online trajectory change
        if (0 && newTrajectoryAwaits == true)
        {
          ROS_ERROR("new trajectory pops up!, track new one");

          double xin[6] = { -newGoalDeltaX, -newGoalDeltaY, -newGoalDeltaZ, Vx, Vy, Vz}; // initial point
          VX vec = GaoModel.evalQP(xin);

          trj.generateTrajectory(vec);

          double refY = Y;
          double refX = X;
          double refZ = Z;

          trj.setPositionOffset(refX - xin[0], refY - xin[1], refZ - xin[2]);
          trj.setTimeOffset(ros::Time::now().toSec() - trajectoryStartTime);
          newTrajectoryAwaits = false;
          //trajectoryStartTime = ros::Time::now().toSec();
        }
#elif TRJ == 3
        //online trajectory change
        if (0 && newTrajectoryAwaits == true)
        {
          ROS_ERROR("new trajectory pops up!, track new one");

          double xin[7] = { -newGoalDeltaX, -newGoalDeltaY, -newGoalDeltaZ, Vx, Vy, Vz, log10(weight)}; // initial point
          VX vec = GaoModel.evalQP(xin);
          MX acce = acceGen(vec);

          trj.generateTrajectory(vec, acce);

          double refY = Y;
          double refX = X;
          double refZ = Z;

          trj.setPositionOffset(refX - xin[0], refY - xin[1], refZ - xin[2]);
          trj.setTimeOffset(ros::Time::now().toSec() - trajectoryStartTime);
          newTrajectoryAwaits = false;
          //trajectoryStartTime = ros::Time::now().toSec();
        }
#endif
        if (0 && (abs(targetDrone.pose.position.x) >= 0.06 || abs(targetDrone.pose.position.y) >= 0.06 || abs(targetDrone.pose.position.z) >= 0.06))
        {
          ROS_ERROR("position way off!, planning new one");

#if TRJ == 3
          double endPx, endPy, endPz;
          trj.getEndPosition(endPx, endPy, endPz);
          double xin[7] = { -(endPx - X), -(endPy - Y), -(endPz - Z), Vx, Vy, Vz, log10(weight)}; // initial point
          VX vec = GaoModel.evalQP(xin);
          MX acce = acceGen(vec);

          trj.generateTrajectory(vec, acce);

          double refX = X, refY = Y, refZ = Z;

          trj.setPositionOffset(refX - xin[0], refY - xin[1], refZ - xin[2]);
          trj.setTimeOffset(ros::Time::now().toSec() - trajectoryStartTime);
#endif
        }

        trajectoryTrackTime = ros::Time::now().toSec() - trajectoryStartTime;
        double targetPx, targetPy, targetPz, targetVx, targetVy, targetVz, targetAx, targetAy, targetAz, targetPhi, targetTheta, targetPsi;
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

        // ROS_WARN("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f", trajectoryTrackTime, targetPx, targetPy, targetPz, X, Y, Z, targetVx, targetVy, targetVz, rawVx, rawVy, rawVz, Vx, Vy, Vz, msg.linear.z, msg.linear.x, msg.linear.y, targetAx, targetAy, targetAz);
        logFile << trajectoryTrackTime << ", " << targetPx << ", " << targetPy << ", " << targetPz << ", " << X << ", " << Y << ", " << Z  << ", " << targetVx << ", " << targetVy << ", " << targetVz << ", " << rawVx << ", " << rawVy << ", " << rawVz << ", " << Vx << ", " << Vy << ", " << Vz << ", " << msg.linear.z << ", " << msg.linear.x << ", " << msg.linear.y << ", " << targetAx << ", " << targetAy << ", " << targetAz << std::endl;
        
        hoverStateCounter = 0;

        double endPx, endPy, endPz;
        trj.getEndPosition(endPx, endPy, endPz);

        double finalPositionError = sqrt((endPx - X) * (endPx - X) + (endPy - Y) * (endPy - Y) + (endPz - Z) * (endPz - Z));

        if (finalPositionError < 0.05 || trajectoryTrackTime >= trj.getEndTime())
        {
          double trajectoryEndTime = ros::Time::now().toSec();

          ROS_ERROR("Trajectory track complete! Took %fs", trajectoryEndTime - trajectoryStartTime);
          
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
  ros::ServiceServer m_serivceGoHome;
  float m_thrust;
  float m_startZ;
  ros::Subscriber m_Pos;
  ros::Subscriber m_Vel;
  ros::Subscriber m_rawPos;
  ros::Subscriber m_rawVel;

  ros::Subscriber m_subscribeGaosGoal;

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

  double newGoalDeltaX;
  double newGoalDeltaY;
  double newGoalDeltaZ;

  float gravityCompensation;

  double trajectoryStartTime;
  double trajectoryTrackTime;

  int hoverStateCounter;

  PID m_pidXTrack;
  PID m_pidYTrack;
  PID m_pidZTrack;

  //MinimumJerkTrajectory trjX;
  //MinimumJerkTrajectory trjY;
  //MinimumJerkTrajectory trjZ;

  //MinimumSnapTrajectory trj2X;
  //MinimumSnapTrajectory trj2Y;
  //MinimumSnapTrajectory trj2Z;

  //CircleTrajectory trj3;
#if TRJ == 2
  TDMST trj;
  double snapVelocity;
#elif TRJ == 1
  mlpModelQP GaoModel;
  Gao trj;
#elif TRJ == 3
  mlpModelQP GaoModel;
  Gao trj;
  double weight;
#endif

  AutoState m_autoState;
  bool goalShallNotChange;
  bool newTrajectoryAwaits;

  std::ofstream logFile;
  int logSequence;

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
