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
#include "Gao.h"
#include "CrazyflieModel.h"
#include "modelLoader.h"
#include "qpMPCUpdate.h"

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
    , dynamicGoalPx(0.0)
    , dynamicGoalPy(0.0)
    , dynamicGoalPz(0.0)
    , dynamicGoalVx(0.0)
    , dynamicGoalVy(0.0)
    , dynamicGoalVz(0.0)
    , lastKeyGoalPx(0.0)
    , lastKeyGoalPy(0.0)
    , lastKeyGoalPz(0.0)
    , m_subscribeDynamicTracking()
    , gravityCompensation(42000)
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
    , m_flyingState(Hovering)
    , goalShallNotChange(false)
    , newDynamicGoalAwaits(false)
    , trj()
    ,mpchorizon(get(n, "mpchorizon"))
    , GaoModel("/home/sun/crazyflie_ws/src/crazyflie_ros/crazyflie_controller/config/pensol_6_100_200_316.json")
    , m_dynamicGoalPos()
    , dynamicTrackingFlag(false)
    , lastTimeNewDynamicGoal(0.0)
    , lastTImeNewTrajectory(0.0)
    , m_dynamicGoalVel()
  {
    ros::NodeHandle nh;
    m_listener.waitForTransform(m_worldFrame, m_frame, ros::Time(0), ros::Duration(10.0));
    m_pubNav = nh.advertise<geometry_msgs::Twist>("cmd_vel", 1);
    m_subscribeGoal = nh.subscribe("goal", 1, &Controller::goalChanged, this);
    m_serviceTakeoff = nh.advertiseService("takeoff", &Controller::takeoff, this);
    m_serviceLand = nh.advertiseService("land", &Controller::land, this);
    m_Pos = nh.subscribe("/Pos", 1, &Controller::positionReceived, this);
    m_Vel = nh.subscribe("/Vel", 1, &Controller::velocityReceived, this);
    m_rawPos = nh.subscribe("/RawPos", 1, &Controller::rawPositionReceived, this);
    m_rawVel = nh.subscribe("/RawVel", 1, &Controller::rawVelocityReceived, this);

    m_dynamicGoalPos = nh.subscribe("/goalPos", 1, &Controller::dynamicGoalPositionReceived, this);
    m_dynamicGoalVel = nh.subscribe("/goalVel", 1, &Controller::dynamicGoalVelocityReceived, this);

    m_subscribeDynamicTracking = nh.subscribe("dynamic_tracking", 1, & Controller::dynamicTrackingSwitch, this);
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

      //init dynamic goal
      dynamicGoalPx = startX;
      dynamicGoalPy = startY;
      dynamicGoalPz = startZ;

      lastKeyGoalPx = startX;
      lastKeyGoalPy = startY;
      lastKeyGoalPz = startZ;
    }
  }

  // no explicit use of the message, every time a message is received, the switch is flipped
  void dynamicTrackingSwitch(const geometry_msgs::Vector3& msg)
  {
    dynamicTrackingFlag = ! dynamicTrackingFlag;

    if (dynamicTrackingFlag == false)
    {
      ROS_ERROR("Dynamic tracking is off!");
    }
    else
    {
      ROS_ERROR("Dynamic tracking is on!");
    }
  }

  void dynamicGoalPositionReceived(const geometry_msgs::Vector3& msg)
  {
    this->dynamicGoalPx = msg.x;
    this->dynamicGoalPy = msg.y;
    this->dynamicGoalPz = msg.z + 1.0;

    newDynamicGoalAwaits = true;

    //if(abs(this->dynamicGoalPx - lastKeyGoalPx) > 0.4 || abs(this->dynamicGoalPy - lastKeyGoalPy) > 0.4 || abs(this->dynamicGoalPz - lastKeyGoalPz) > 0.4)
    //{
      //newDynamicGoalAwaits = true;
      //ROS_WARN("key goal updated!");

      //lastKeyGoalPx = this->dynamicGoalPx;
      //lastKeyGoalPy = this->dynamicGoalPy;
      //lastKeyGoalPz = this->dynamicGoalPz;
    //}

    //ROS_WARN("new dynamic goal position received:(%f, %f, %f)", msg.x, msg.y, msg.z);
  }

  void dynamicGoalVelocityReceived(const geometry_msgs::Vector3& msg)
  {
    this -> dynamicGoalVx = msg.x;
    this -> dynamicGoalVy = msg.y;
    this -> dynamicGoalVz = msg.z;
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
        m_state = Flying;
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
    case Flying:
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

      if (positionError < 0.08 && velocityError < 0.10 && m_flyingState == Hovering)
      {
        hoverStateCounter++;
      }
      else
      {
        hoverStateCounter = 0;
      }

      if (hoverStateCounter >= 20)
      {
        // gravity compenstation reset
        gravityCompensation += m_pidZ.getIntegral() * m_pidZ.ki();
        m_pidZ.setIntegral(0.0);
        hoverStateCounter = 0;
      }

      // new FSM
      if (m_flyingState == Hovering) // hovering
      {
        //ROS_INFO("Hovering! Goal: %f, %f, %f, %f, %f, %f", m_goal.pose.position.x, m_goal.pose.position.y, m_goal.pose.position.z, targetDrone.pose.position.x, targetDrone.pose.position.y, targetDrone.pose.position.z);

        msg.linear.x = m_pidX.update(targetDrone.pose.position.x, 0 - Vx);
        msg.linear.y = m_pidY.update(targetDrone.pose.position.y, 0 - Vy);
        msg.linear.z = gravityCompensation + m_pidZ.update(targetDrone.pose.position.z, 0 - Vz);

        if (dynamicTrackingFlag == true)
        {
          m_flyingState = Following;
        }

      }
      else if (m_flyingState == Following)
      {
        // check if it is neccessary to generate a new trajectory
        double distanceToGoal = (dynamicGoalPx - X) * (dynamicGoalPx - X) + (dynamicGoalPy - Y) * (dynamicGoalPy - Y) + (dynamicGoalPz - Z) * (dynamicGoalPz - Z);
        distanceToGoal = sqrt(distanceToGoal);

        if (distanceToGoal < 0.50) // PID is enough for small movements
        {
          ROS_INFO("small F");

          msg.linear.x = m_pidX.update(dynamicGoalPx - X, dynamicGoalVx - Vx);
          msg.linear.y = m_pidY.update(dynamicGoalPy - Y, dynamicGoalVy - Vy);
          msg.linear.z = gravityCompensation + m_pidZ.update(dynamicGoalPz - Z, dynamicGoalVz - Vz);

          m_goal.pose.position.x = dynamicGoalPx;
          m_goal.pose.position.y = dynamicGoalPy;
          m_goal.pose.position.z = dynamicGoalPz;

        }
        else
        {
          ROS_INFO("big F");

          // plan new trajectory
          if (newDynamicGoalAwaits == true && ros::Time::now().toSec() - lastTImeNewTrajectory > 0.4)
          {
            // generate a new trajectory
            double xin[6] = { -(this->dynamicGoalPx - X), -(this->dynamicGoalPy - Y), -(this->dynamicGoalPz - Z), Vx, Vy, Vz}; // initial point
            VX vec = GaoModel.evalQP(xin);
            MX acce = qpMPC::acceGen(vec);
            trj.generateTrajectory(vec, acce, mpchorizon);
            trj.setPositionOffset(X - xin[0], Y - xin[1], Z - xin[2]);

            newDynamicGoalAwaits = false;
            trajectoryStartTime = ros::Time::now().toSec();
            lastTImeNewTrajectory = trajectoryStartTime;
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

          //ROS_WARN("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f",
                   //trajectoryTrackTime, targetPx, targetPy, targetPz, X, Y, Z, targetVx, targetVy, targetVz, rawVx, rawVy,
                   //rawVz, Vx, Vy, Vz, msg.linear.z, msg.linear.x, msg.linear.y, targetAx, targetAy, targetAz);

        }

        if (dynamicTrackingFlag == false)
        {
          m_flyingState = Hovering;

          m_goal.pose.position.x = m_goal.pose.position.x;
          m_goal.pose.position.y = m_goal.pose.position.y;
          m_goal.pose.position.z = m_goal.pose.position.z;

          m_pidXTrack.reset();
          m_pidYTrack.reset();
          m_pidZTrack.reset();
        }
      }
      else
      {
        ROS_ERROR("Flying sub state error!");
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
      m_flyingState = Hovering;
    }
    break;
    }
  }

private:

  enum State
  {
    Idle = 0,
    Flying = 1,
    TakingOff = 2,
    Landing = 3,
  };

  enum FlyingState
  {
    Hovering = 0,
    Following = 1,
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
  float m_thrust;
  float m_startZ;
  ros::Subscriber m_Pos;
  ros::Subscriber m_Vel;
  ros::Subscriber m_rawPos;
  ros::Subscriber m_rawVel;

  ros::Subscriber m_subscribeDynamicTracking;
  ros::Subscriber m_dynamicGoalPos;
  ros::Subscriber m_dynamicGoalVel;

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

  double dynamicGoalPx;
  double dynamicGoalPy;
  double dynamicGoalPz;

  double dynamicGoalVx;
  double dynamicGoalVy;
  double dynamicGoalVz;

  double lastKeyGoalPx;
  double lastKeyGoalPy;
  double lastKeyGoalPz;

  float gravityCompensation;

  double trajectoryStartTime;
  double trajectoryTrackTime;

  int hoverStateCounter;

  PID m_pidXTrack;
  PID m_pidYTrack;
  PID m_pidZTrack;

  mlpModelMPCQP GaoModel;
  Gao trj;
  double mpchorizon;

  FlyingState m_flyingState;
  bool goalShallNotChange;
  bool newDynamicGoalAwaits;
  bool dynamicTrackingFlag;

  double lastTimeNewDynamicGoal;
  double lastTImeNewTrajectory;

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
