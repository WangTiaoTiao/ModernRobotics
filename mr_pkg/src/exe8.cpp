#include "mr_h.h"
#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/joint_state.hpp"

using namespace std::chrono_literals;

class JointStatePublisher : public rclcpp::Node
{
public:
  using JointState = sensor_msgs::msg::JointState;

  JointStatePublisher() : Node("exe8_node")
  {
    js_pub_ = this->create_publisher<JointState>("/joint_states", 10);
    
    js_msg_ = std::make_shared<JointState>();
    js_msg_->name = {"shoulder_1_joint", "shoulder_2_joint", "elbow_joint", 
                     "wrist_1_joint", "wrist_2_joint", "wrist_3_joint"};
    js_msg_->header.stamp = this->now();
    js_msg_->position = {0, 0, 0, 0, 0, 0};

    UR5INFO();

    timer_ = this->create_wall_timer(10ms, std::bind(&JointStatePublisher::timer_cb, this));

    RCLCPP_INFO(this->get_logger(), "Create exe8_node !");
  }

  ~JointStatePublisher() {
    RCLCPP_INFO(this->get_logger(), "Destory exe8_node !");
  }

private:

  void step(const VectorXf& dd_thetalist) {
    float dt = 0.01;
    d_thetalist_ = d_thetalist_ + dt * dd_thetalist;
    thetalist_ = thetalist_ + dt * d_thetalist_;
  }

  void timer_cb() {
  
    auto dd_thetalist = ForwardDynamics(thetalist_, d_thetalist_, taulist_, g_, f_tip, Mlist_, Glist_, Slist_);
    
    step(dd_thetalist);

    js_msg_->header.stamp = this->now();
    for (int i = 0; i < 6; ++i) js_msg_->position[i] = deg2rad(rad2deg(thetalist_[i]));

    js_pub_->publish(*js_msg_);
  }

  void UR5INFO() {
    float l1 = 0.425, l2 = 0.39225, w1 = 0.13585, w2 = 0.093, w3 = 0.0823, h1 = 0.089159, h2 = 0.09465;

    Slist_.reserve(6);
    VectorXf s(6);
    s << 0, 0, 1, 0, 0, 0;
    Slist_.emplace_back(s);
    s << 0, 1, 0, -h1, 0, 0;
    Slist_.emplace_back(s);
    s << 0, 1, 0, -h1, 0, l1;
    Slist_.emplace_back(s);
    s << 0, 1, 0, -h1, 0, l1 + l2;
    Slist_.emplace_back(s);
    s << 0, 0, -1, -w1, l1 + l2, 0;
    Slist_.emplace_back(s);
    s << 0, 1, 0, h2 - h1, 0, l1 + l2;
    Slist_.emplace_back(s);

    Mlist_.reserve(7);
    Mlist_.emplace_back(Matrix4f{{1, 0, 0,  0}, {0, 1, 0,  0}, { 0,  0, 1,  h1}, {0, 0, 0, 1}});
    Mlist_.emplace_back(Matrix4f{{0, 0, 1,  0}, {0, 1, 0, w1}, {-1,  0, 0,   0}, {0, 0, 0, 1}});
    Mlist_.emplace_back(Matrix4f{{1, 0, 0,  0}, {0, 1, 0,  0}, { 0,  0, 1,  l1}, {0, 0, 0, 1}});
    Mlist_.emplace_back(Matrix4f{{0, 0, 1,  0}, {0, 1, 0,  0}, {-1,  0, 0,  l2}, {0, 0, 0, 1}});
    Mlist_.emplace_back(Matrix4f{{1, 0, 0,  0}, {0, 1, 0, w2}, { 0,  0, 1,   0}, {0, 0, 0, 1}});
    Mlist_.emplace_back(Matrix4f{{1, 0, 0,  0}, {0, 1, 0,  0}, { 0,  0, 1,  h2}, {0, 0, 0, 1}});
    Mlist_.emplace_back(Matrix4f{{1, 0, 0,  0}, {0, 0, 1, w3}, { 0, -1, 0,   0}, {0, 0, 0, 1}});

    Glist_.reserve(6);
    Glist_.emplace_back(getG_cylinder(3.7,    0.06, 0.15));
    Glist_.emplace_back(getG_cylinder(8.393,  0.06, 0.56));
    Glist_.emplace_back(getG_cylinder(2.275,  0.06, 0.50));
    Glist_.emplace_back(getG_cylinder(1.219,  0.6,  0.12));
    Glist_.emplace_back(getG_cylinder(1.219,  0.6,  0.12));
    Glist_.emplace_back(getG_cylinder(0.1879, 0.6,  0.12));

    thetalist_ = Eigen::VectorXf(6);
    thetalist_ << 0, 0, 0, 0, 0, 0;
  
    d_thetalist_ = Eigen::VectorXf(6);
    d_thetalist_ << 0, 0, 0, 0, 0, 0;
  
    taulist_ = Eigen::VectorXf(6);
    taulist_ << 0, 0, 0, 0, 0, 0;

    f_tip = Eigen::VectorXf(6);
    f_tip << 0, 0, 0, 0, 0, 0;

    g_ << 0, 0, -9.81;
  }

  MatrixXf getG_cylinder(const float &m, const float &r, const float &l) const {
    MatrixXf g = MatrixXf::Zero(6, 6);
    float ixy = m * (3 * r * r + l * l) / 12;
    float iz = m * r * r * 0.5;
    g.diagonal() << ixy, ixy, iz, m, m, m;
    return g;
  }

  std::vector<Eigen::VectorXf> Slist_;
  std::vector<Eigen::Matrix4f> Mlist_;
  std::vector<Eigen::MatrixXf> Glist_;
  Eigen::VectorXf thetalist_;
  Eigen::VectorXf d_thetalist_;
  Eigen::VectorXf taulist_;
  Eigen::Vector3f g_;
  Eigen::VectorXf f_tip; 
  
  rclcpp::Publisher<JointState>::SharedPtr js_pub_;
  rclcpp::TimerBase::SharedPtr timer_;
  std::shared_ptr<JointState> js_msg_;
};

int main(int argc, char ** argv)
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<JointStatePublisher>());
  rclcpp::shutdown();
  return 0;
}
